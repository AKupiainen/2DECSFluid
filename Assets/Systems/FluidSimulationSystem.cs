using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Transforms;
using UnityEngine; // Added for Debug.Log

[UpdateInGroup(typeof(SimulationSystemGroup))]
public partial class FluidSimulationSystem : SystemBase
{
    private EntityQuery fluidParticleQuery;
    
    protected override void OnCreate()
    {
        fluidParticleQuery = GetEntityQuery(typeof(FluidParticle), typeof(LocalTransform));
    }

    // Improved smoothing kernels with better stability
    [BurstCompile]
    private static class SmoothingKernels
    {
        // Poly6 kernel for density
        public static float Poly6(float2 r, float h)
        {
            float r2 = math.lengthsq(r);
            float h2 = h * h;
            
            if (r2 >= h2) return 0;
            
            float coef = 4.0f / (math.PI * math.pow(h, 8));
            return coef * math.pow(h2 - r2, 3);
        }
        
        // Spiky kernel gradient for pressure - strengthened for better repulsion
        public static float2 SpikyGradient(float2 r, float h)
        {
            float d = math.length(r);
            
            if (d > h || d < 1e-6f) return float2.zero;
            
            // Increased coefficient for stronger repulsion
            float coef = -45.0f / (math.PI * math.pow(h, 5));
            
            // Enhanced gradient strength at close distances
            float factor;
            if (d < 0.2f * h)
            {
                // Much stronger at very close distances
                factor = coef * math.pow(h - d, 2) / math.max(d, 0.01f * h) * 3.0f;
            }
            else
            {
                factor = coef * math.pow(h - d, 2) / d;
            }
            
            return r * factor;
        }
        
        // Viscosity Laplacian for viscosity forces
        public static float ViscosityLaplacian(float2 r, float h)
        {
            float d = math.length(r);
            
            if (d >= h) return 0;
            
            float coef = 40.0f / (math.PI * math.pow(h, 5));
            return coef * (h - d);
        }
    }

    [BurstCompile]
    private struct CalculateDensityJob : IJobParallelFor
    {
        [ReadOnly] public NativeArray<LocalTransform> Positions;
        [ReadOnly] public NativeArray<FluidParticle> ParticlesRead;
        [WriteOnly] public NativeArray<FluidParticle> ParticlesWrite;
        [ReadOnly] public FluidSimulationParameters Params;
        
        public void Execute(int index)
        {
            FluidParticle particle = ParticlesRead[index];
            float2 position = Positions[index].Position.xy;
            
            float density = 0;
            
            // Count very close particles for strong repulsion
            int closeParticleCount = 0;
            
            // Sum density contributions from all particles
            for (int i = 0; i < Positions.Length; i++)
            {
                if (i == index) continue;
                
                float2 otherPosition = Positions[i].Position.xy;
                float2 r = position - otherPosition;
                float d = math.length(r);
                
                // Track particles that are too close
                if (d < 0.25f * Params.SmoothingRadius)
                {
                    closeParticleCount++;
                }
                
                density += ParticlesRead[i].Mass * SmoothingKernels.Poly6(r, Params.SmoothingRadius);
            }
            
            // Add self-contribution to density
            density += particle.Mass * SmoothingKernels.Poly6(float2.zero, Params.SmoothingRadius);
            
            // Add extra density for very close particles to trigger stronger repulsion
            if (closeParticleCount > 0)
            {
                density += closeParticleCount * particle.Mass * 2.0f;
            }
            
            // Ensure reasonable density
            density = math.max(density, Params.TargetDensity * 0.1f);
            
            // Calculate pressure using NON-LINEAR pressure response for stronger repulsion
            float pressureDiff = density - Params.TargetDensity;
            
            // More aggressive pressure response at high densities
            float pressure;
            if (pressureDiff > 0)
            {
                // Use a quadratic pressure response for stronger repulsion
                pressure = Params.PressureCoefficient * pressureDiff * pressureDiff * 1.5f;
            }
            else
            {
                pressure = 0;
            }
            
            // Create a new particle with updated values
            FluidParticle updatedParticle = particle;
            updatedParticle.Density = density;
            updatedParticle.Pressure = pressure;
            
            ParticlesWrite[index] = updatedParticle;
        }
    }

    [BurstCompile]
    private struct CalculateForceJob : IJobParallelFor
    {
        [ReadOnly] public NativeArray<LocalTransform> Positions;
        [ReadOnly] public NativeArray<FluidParticle> ParticlesRead;
        [WriteOnly] public NativeArray<FluidParticle> ParticlesWrite;
        [ReadOnly] public FluidSimulationParameters Params;
        
        public void Execute(int index)
        {
            FluidParticle particle = ParticlesRead[index];
            float2 position = Positions[index].Position.xy;
            
            float2 pressureForce = float2.zero;
            float2 viscosityForce = float2.zero;
            
            // Add direct repulsion force for particles that are too close
            float2 directRepulsionForce = float2.zero;
            
            // Define the minimum allowed distance (effective particle radius)
            float minDistance = Params.ParticleRadius > 0 ? 
                               Params.ParticleRadius : 
                               0.2f * Params.SmoothingRadius; // Default if not set
            
            for (int i = 0; i < Positions.Length; i++)
            {
                if (i == index) continue;
                
                float2 otherPosition = Positions[i].Position.xy;
                float2 r = position - otherPosition;
                float d = math.length(r);
                
                // Skip if too far
                if (d >= Params.SmoothingRadius) continue;
                
                // Hard constraint repulsion for overlapping particles
                // This creates an extremely strong repulsion when particles are closer than minDistance
                if (d < minDistance)
                {
                    float overlap = minDistance - d;
                    float2 repulsionDir = r / math.max(d, 0.001f); // Normalized direction vector
                    
                    // Extremely strong repulsion with cubic falloff for hard constraint
                    // The strength can be controlled through the Inspector with OverlapPreventionStrength
                    float repulsionMagnitude = 50.0f * Params.OverlapPreventionStrength * 
                                              Params.PressureCoefficient * 
                                              math.pow(overlap / minDistance, 3);
                    
                    directRepulsionForce += repulsionDir * repulsionMagnitude;
                    
                    // Log serious overlaps
                    if (overlap > 0.5f * minDistance)
                    {
                        // This indicates a severe overlap that shouldn't happen if simulation is stable
                        // We can't log from jobs, but we could set a flag here if needed
                    }
                }
                // Regular pressure force calculation for normal separation
                else 
                {
                    float2 pressureGradient = SmoothingKernels.SpikyGradient(r, Params.SmoothingRadius);
                    if (math.lengthsq(pressureGradient) > 1e-10f)
                    {
                        // Improved pressure calculation using particle pressures
                        float combinedPressure = particle.Pressure + ParticlesRead[i].Pressure;
                        float pressureTerm = combinedPressure * 0.5f * ParticlesRead[i].Mass / ParticlesRead[i].Density;
                        pressureForce += pressureTerm * pressureGradient;
                    }
                    
                    // Viscosity force with better stability
                    float viscosityLaplacian = SmoothingKernels.ViscosityLaplacian(r, Params.SmoothingRadius);
                    viscosityForce += ParticlesRead[i].Mass * 
                                     (ParticlesRead[i].Velocity - particle.Velocity) / 
                                     ParticlesRead[i].Density * viscosityLaplacian;
                }
            }
            
            // Apply forces with better scaling
            float2 totalForce = float2.zero;
            
            // Apply pressure force (negative because it pushes particles apart)
            totalForce -= pressureForce;
            
            // Add direct repulsion (this is crucial for preventing stacking)
            totalForce += directRepulsionForce;
            
            // Add viscosity force
            totalForce += Params.ViscosityCoefficient * particle.Viscosity * viscosityForce;
            
            // Apply gravity
            totalForce += new float2(0, -Params.GravityForce) * particle.Density;
            
            // Clamp total force to prevent instability
            float maxForce = 4000.0f; // Increased for stronger repulsion
            float forceMagnitude = math.length(totalForce);
            if (forceMagnitude > maxForce)
            {
                totalForce = totalForce * (maxForce / forceMagnitude);
            }
            
            // Create a new particle with updated force
            FluidParticle updatedParticle = particle;
            updatedParticle.Force = totalForce;
            
            ParticlesWrite[index] = updatedParticle;
        }
    }

    [BurstCompile]
    private struct IntegrateParticlesJob : IJobParallelFor
    {
        [ReadOnly] public NativeArray<FluidParticle> Particles;
        public NativeArray<LocalTransform> Positions;
        public NativeArray<float2> NewVelocities;
        [ReadOnly] public FluidSimulationParameters Params;
        public float DeltaTime;
        
        public Unity.Mathematics.Random RandomGenerator;
        
        public void Execute(int index)
        {
            FluidParticle particle = Particles[index];
            LocalTransform transform = Positions[index];
            
            // Calculate acceleration from force with clamping to avoid extreme values
            float2 acceleration = particle.Force / math.max(particle.Density, 1e-5f);
            
            // Clamp acceleration to reasonable limits
            float maxAcceleration = 100.0f; // Adjust based on your simulation scale
            float accelerationMagnitude = math.length(acceleration);
            if (accelerationMagnitude > maxAcceleration)
            {
                acceleration = acceleration * (maxAcceleration / accelerationMagnitude);
            }
            
            // Update velocity with proper time step - using adaptive time step for stability
            // Reduce timestep if acceleration is very large
            float adaptiveDt = accelerationMagnitude > 50.0f ? DeltaTime * 0.5f : DeltaTime;
            float2 newVelocity = particle.Velocity + acceleration * adaptiveDt;
            
            // Apply velocity damping to increase stability
            newVelocity *= 0.99f; // Small damping factor
            
            // Clamp velocity to reasonable limits
            float maxVelocity = 40.0f; // Increased to allow for proper collision response
            float velocityMagnitude = math.length(newVelocity);
            if (velocityMagnitude > maxVelocity)
            {
                newVelocity = newVelocity * (maxVelocity / velocityMagnitude);
            }
            
            // Update position using semi-implicit Euler integration
            float2 newPosition = transform.Position.xy + newVelocity * adaptiveDt;
            
            // Handle boundary collisions with proper velocity response - enhanced collision strength
            bool collided = false;
            
            // Boundary collision using more aggressive response
            float collisionStrength = Params.BoundaryDamping * 1.5f;
            float penetrationResponse = 0.02f; // Push particles away from boundary by this amount
            
            if (newPosition.x < Params.BoundaryMin.x)
            {
                // Move particle away from boundary and apply strong reflection to velocity
                newPosition.x = Params.BoundaryMin.x + penetrationResponse;
                newVelocity.x = math.abs(newVelocity.x) * collisionStrength;
                collided = true;
            }
            else if (newPosition.x > Params.BoundaryMax.x)
            {
                newPosition.x = Params.BoundaryMax.x - penetrationResponse;
                newVelocity.x = -math.abs(newVelocity.x) * collisionStrength;
                collided = true;
            }
            
            if (newPosition.y < Params.BoundaryMin.y)
            {
                newPosition.y = Params.BoundaryMin.y + penetrationResponse;
                newVelocity.y = math.abs(newVelocity.y) * collisionStrength;
                collided = true;
            }
            else if (newPosition.y > Params.BoundaryMax.y)
            {
                newPosition.y = Params.BoundaryMax.y - penetrationResponse;
                newVelocity.y = -math.abs(newVelocity.y) * collisionStrength;
                collided = true;
            }
            
            // If collision occurred, add a very small random jitter to avoid particles aligning on boundaries
            if (collided)
            {
                // Create a new random state by hashing the particle index
                uint seed = RandomGenerator.state + (uint)(index * 747796405);
                Unity.Mathematics.Random random = new Unity.Mathematics.Random(seed);
                
                // Apply small jitter to avoid particles lining up on boundaries
                newPosition += random.NextFloat2(-0.001f, 0.001f);
            }
            
            // Update position with new coordinates
            transform.Position = new float3(newPosition.x, newPosition.y, 0);
            Positions[index] = transform;
            
            // Store new velocity for the update job
            NewVelocities[index] = newVelocity;
        }
    }

    [BurstCompile]
    private struct UpdateParticlesJob : IJobParallelFor
    {
        [ReadOnly] public NativeArray<FluidParticle> ParticlesRead;
        [WriteOnly] public NativeArray<FluidParticle> ParticlesWrite;
        [ReadOnly] public NativeArray<float2> NewVelocities;
        
        public void Execute(int index)
        {
            FluidParticle updatedParticle = ParticlesRead[index];
            updatedParticle.Velocity = NewVelocities[index];
            
            ParticlesWrite[index] = updatedParticle;
        }
    }
    
    protected override void OnUpdate()
    {
        if (fluidParticleQuery.IsEmpty)
            return;
        
        var fluidParams = SystemAPI.GetSingleton<FluidSimulationParameters>();
        
        // Use a fixed smaller timestep or adaptive timestep for stability
        float deltaTime = math.min(0.016f, UnityEngine.Time.deltaTime);
        float stableTimeStep = math.min(deltaTime, 0.005f); // Cap at 5ms for stability
        
        // Debug the boundary values and parameters to make sure they are reasonable
        Debug.Log($"Boundaries: Min({fluidParams.BoundaryMin.x}, {fluidParams.BoundaryMin.y}) Max({fluidParams.BoundaryMax.x}, {fluidParams.BoundaryMax.y})");
        Debug.Log($"Particle radius: {(fluidParams.ParticleRadius > 0 ? fluidParams.ParticleRadius : 0.2f * fluidParams.SmoothingRadius)}");
        Debug.Log($"Time step: {stableTimeStep}");
        Debug.Log($"Smoothing radius: {fluidParams.SmoothingRadius}");
        
        // Create native collections to store data during the simulation
        var allParticles = new NativeList<FluidParticle>(Allocator.TempJob);
        var allTransforms = new NativeList<LocalTransform>(Allocator.TempJob);
        var allEntities = new NativeList<Entity>(Allocator.TempJob);
        
        // Collect all particle data
        Entities
            .WithStoreEntityQueryInField(ref fluidParticleQuery)
            .ForEach((Entity entity, in FluidParticle particle, in LocalTransform transform) => {
                allParticles.Add(particle);
                allTransforms.Add(transform);
                allEntities.Add(entity);
            }).Run();
            
        int particleCount = allParticles.Length;
        
        if (particleCount == 0)
        {
            allParticles.Dispose();
            allTransforms.Dispose();
            allEntities.Dispose();
            return;
        }
            
        // Convert to NativeArray for job compatibility
        var particlesRead = new NativeArray<FluidParticle>(allParticles.AsArray(), Allocator.TempJob);
        var particlesAfterDensity = new NativeArray<FluidParticle>(particleCount, Allocator.TempJob);
        var particlesAfterForce = new NativeArray<FluidParticle>(particleCount, Allocator.TempJob);
        var particlesFinal = new NativeArray<FluidParticle>(particleCount, Allocator.TempJob);
        var transformArray = new NativeArray<LocalTransform>(allTransforms.AsArray(), Allocator.TempJob);
        var newVelocities = new NativeArray<float2>(particleCount, Allocator.TempJob);
        
        // Step 1: Calculate density and pressure
        var densityJob = new CalculateDensityJob
        {
            Positions = transformArray,
            ParticlesRead = particlesRead,
            ParticlesWrite = particlesAfterDensity,
            Params = fluidParams
        };
        
        // Step 2: Calculate forces based on density and pressure
        var forceJob = new CalculateForceJob
        {
            Positions = transformArray,
            ParticlesRead = particlesAfterDensity, // Read from previous job output
            ParticlesWrite = particlesAfterForce,  // Write to new buffer
            Params = fluidParams
        };
        
        // Step 3: Integrate particles to calculate new positions and velocities
        var integrateJob = new IntegrateParticlesJob
        {
            Particles = particlesAfterForce,     // Read from previous job output
            Positions = transformArray,
            NewVelocities = newVelocities,
            Params = fluidParams,
            DeltaTime = stableTimeStep,
            // Initialize random number generator with a stable seed from the main thread
            RandomGenerator = new Unity.Mathematics.Random((uint)(System.DateTimeOffset.Now.ToUnixTimeMilliseconds() & 0xFFFFFFFF))
        };
        
        // Step 4: Update particle velocities
        var updateParticlesJob = new UpdateParticlesJob
        {
            ParticlesRead = particlesAfterForce,  // Read from previous job output
            ParticlesWrite = particlesFinal,      // Final output
            NewVelocities = newVelocities
        };
        
        // Schedule jobs with dependencies
        JobHandle densityHandle = densityJob.Schedule(particleCount, 64);
        JobHandle forceHandle = forceJob.Schedule(particleCount, 64, densityHandle);
        JobHandle integrateHandle = integrateJob.Schedule(particleCount, 64, forceHandle);
        JobHandle updateHandle = updateParticlesJob.Schedule(particleCount, 64, integrateHandle);
        
        // Wait for jobs to complete
        updateHandle.Complete();
        
        // Create a command buffer to apply changes to the entities
        var ecb = new EntityCommandBuffer(Allocator.TempJob);
        
        // Apply changes to entities
        for (int i = 0; i < particleCount; i++)
        {
            ecb.SetComponent(allEntities[i], particlesFinal[i]);
            ecb.SetComponent(allEntities[i], transformArray[i]);
        }
        
        // Play back the command buffer
        ecb.Playback(EntityManager);
        
        // Clean up resources
        ecb.Dispose();
        particlesRead.Dispose();
        particlesAfterDensity.Dispose();
        particlesAfterForce.Dispose();
        particlesFinal.Dispose();
        transformArray.Dispose();
        newVelocities.Dispose();
        allParticles.Dispose();
        allTransforms.Dispose();
        allEntities.Dispose();
    }
}