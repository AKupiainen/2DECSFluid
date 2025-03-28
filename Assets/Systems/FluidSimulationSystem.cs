using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Transforms;

[UpdateInGroup(typeof(SimulationSystemGroup))]
public partial class FluidSimulationSystem : SystemBase
{
    private EntityQuery _fluidParticleQuery;
    
    protected override void OnCreate()
    {
        _fluidParticleQuery = GetEntityQuery(typeof(FluidParticle), typeof(LocalTransform));
    }
    
    [BurstCompile]
    private static class SmoothingKernels
    {
        public static float Poly6(float2 r, float h)
        {
            float r2 = math.lengthsq(r);
            float h2 = h * h;
            
            if (r2 >= h2)
            {
                return 0;
            }

            float coef = 4.0f / (math.PI * math.pow(h, 8));
            return coef * math.pow(h2 - r2, 3);
        }
        
        public static float2 SpikyGradient(float2 r, float h)
        {
            float d = math.length(r);
            
            if (d > h || d < 1e-6f)
            {
                return float2.zero;
            }

            float coef = -45.0f / (math.PI * math.pow(h, 5));
            
            float factor;
            if (d < 0.2f * h)
            {
                factor = coef * math.pow(h - d, 2) / math.max(d, 0.01f * h) * 3.0f;
            }
            else
            {
                factor = coef * math.pow(h - d, 2) / d;
            }
            
            return r * factor;
        }
        
        public static float ViscosityLaplacian(float2 r, float h)
        {
            float d = math.length(r);
            
            if (d >= h)
            {
                return 0;
            }

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
            int closeParticleCount = 0;
            
            for (int i = 0; i < Positions.Length; i++)
            {
                if (i == index)
                {
                    continue;
                }

                float2 otherPosition = Positions[i].Position.xy;
                float2 r = position - otherPosition;
                float d = math.length(r);
                
                if (d < 0.25f * Params.SmoothingRadius)
                {
                    closeParticleCount++;
                }
                
                density += ParticlesRead[i].Mass * SmoothingKernels.Poly6(r, Params.SmoothingRadius);
            }
            
            density += particle.Mass * SmoothingKernels.Poly6(float2.zero, Params.SmoothingRadius);
            
            if (closeParticleCount > 0)
            {
                density += closeParticleCount * particle.Mass * 2.0f;
            }
            
            density = math.max(density, Params.TargetDensity * 0.1f);
            
            float pressureDiff = density - Params.TargetDensity;
            float pressure;
            
            if (pressureDiff > 0)
            {
                pressure = Params.PressureCoefficient * pressureDiff * pressureDiff * 1.5f;
            }
            else
            {
                pressure = 0;
            }
            
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
            
            float2 directRepulsionForce = float2.zero;
            
            float minDistance = Params.ParticleRadius > 0 ? 
                               Params.ParticleRadius : 
                               0.2f * Params.SmoothingRadius;
            
            for (int i = 0; i < Positions.Length; i++)
            {
                if (i == index)
                {
                    continue;
                }

                float2 otherPosition = Positions[i].Position.xy;
                float2 r = position - otherPosition;
                float d = math.length(r);
                
                if (d >= Params.SmoothingRadius)
                {
                    continue;
                }
                
                if (d < minDistance)
                {
                    float overlap = minDistance - d;
                    float2 repulsionDir = r / math.max(d, 0.001f); 
                    
                    float repulsionMagnitude = 50.0f * Params.OverlapPreventionStrength * 
                                              Params.PressureCoefficient * 
                                              math.pow(overlap / minDistance, 3);
                    
                    directRepulsionForce += repulsionDir * repulsionMagnitude;
                    
                    if (overlap > 0.5f * minDistance)
                    {
                        // This indicates a severe overlap that shouldn't happen if simulation is stable
                        // We can't log from jobs, but we could set a flag here if needed
                    }
                }
                else 
                {
                    float2 pressureGradient = SmoothingKernels.SpikyGradient(r, Params.SmoothingRadius);
                    
                    if (math.lengthsq(pressureGradient) > 1e-10f)
                    {
                        float combinedPressure = particle.Pressure + ParticlesRead[i].Pressure;
                        float pressureTerm = combinedPressure * 0.5f * ParticlesRead[i].Mass / ParticlesRead[i].Density;
                        pressureForce += pressureTerm * pressureGradient;
                    }
                    
                    float viscosityLaplacian = SmoothingKernels.ViscosityLaplacian(r, Params.SmoothingRadius);
                    viscosityForce += ParticlesRead[i].Mass * 
                                     (ParticlesRead[i].Velocity - particle.Velocity) / 
                                     ParticlesRead[i].Density * viscosityLaplacian;
                }
            }
            
            float2 totalForce = float2.zero;
            totalForce -= pressureForce;
            
            totalForce += directRepulsionForce;
            totalForce += Params.ViscosityCoefficient * particle.Viscosity * viscosityForce;
            
            totalForce += new float2(0, -Params.GravityForce) * particle.Density;
            
            float maxForce = 4000.0f; 
            float forceMagnitude = math.length(totalForce);
            
            if (forceMagnitude > maxForce)
            {
                totalForce *= (maxForce / forceMagnitude);
            }
            
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
        
        public Random RandomGenerator;
        
        public void Execute(int index)
        {
            FluidParticle particle = Particles[index];
            LocalTransform transform = Positions[index];
            
            float2 acceleration = particle.Force / math.max(particle.Density, 1e-5f);
            
            float maxAcceleration = 100.0f; 
            float accelerationMagnitude = math.length(acceleration);
            
            if (accelerationMagnitude > maxAcceleration)
            {
                acceleration *= (maxAcceleration / accelerationMagnitude);
            }
            
            float adaptiveDt = accelerationMagnitude > 50.0f ? DeltaTime * 0.5f : DeltaTime;
            float2 newVelocity = particle.Velocity + acceleration * adaptiveDt;
            
            newVelocity *= 0.99f;
            
            float maxVelocity = 40.0f; 
            float velocityMagnitude = math.length(newVelocity);
            
            if (velocityMagnitude > maxVelocity)
            {
                newVelocity *= (maxVelocity / velocityMagnitude);
            }
            
            float2 newPosition = transform.Position.xy + newVelocity * adaptiveDt;
            
            bool collided = false;
            
            float collisionStrength = Params.BoundaryDamping * 1.5f;
            float penetrationResponse = 0.02f; 
            
            if (newPosition.x < Params.BoundaryMin.x)
            {
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
            
            if (collided)
            {
                uint seed = RandomGenerator.state + (uint)(index * 747796405);
                Random random = new(seed);
                
                newPosition += random.NextFloat2(-0.001f, 0.001f);
            }
            
            transform.Position = new float3(newPosition.x, newPosition.y, 0);
            Positions[index] = transform;
            
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
        if (_fluidParticleQuery.IsEmpty)
        {
            return;
        }

        FluidSimulationParameters fluidParams = SystemAPI.GetSingleton<FluidSimulationParameters>();
        
        float deltaTime = math.min(fluidParams.DeltaTime, UnityEngine.Time.deltaTime);
        float stableTimeStep = math.min(deltaTime, 0.005f);
        
        NativeList<FluidParticle> allParticles = new(Allocator.TempJob);
        NativeList<LocalTransform> allTransforms = new(Allocator.TempJob);
        NativeList<Entity> allEntities = new(Allocator.TempJob);
        
        Entities
            .WithStoreEntityQueryInField(ref _fluidParticleQuery)
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
        
        NativeArray<FluidParticle> particlesRead = new(allParticles.AsArray(), Allocator.TempJob);
        NativeArray<FluidParticle> particlesAfterDensity = new(particleCount, Allocator.TempJob);
        NativeArray<FluidParticle> particlesAfterForce = new(particleCount, Allocator.TempJob);
        NativeArray<FluidParticle> particlesFinal = new(particleCount, Allocator.TempJob);
        NativeArray<LocalTransform> transformArray = new(allTransforms.AsArray(), Allocator.TempJob);
        NativeArray<float2> newVelocities = new(particleCount, Allocator.TempJob);
        
        CalculateDensityJob densityJob = new()
        {
            Positions = transformArray,
            ParticlesRead = particlesRead,
            ParticlesWrite = particlesAfterDensity,
            Params = fluidParams
        };
        
        CalculateForceJob forceJob = new()
        {
            Positions = transformArray,
            ParticlesRead = particlesAfterDensity, 
            ParticlesWrite = particlesAfterForce,  
            Params = fluidParams
        };
        

        IntegrateParticlesJob integrateJob = new()
        {
            Particles = particlesAfterForce,    
            Positions = transformArray,
            NewVelocities = newVelocities,
            Params = fluidParams,
            DeltaTime = stableTimeStep,
            RandomGenerator = new Random((uint)(System.DateTimeOffset.Now.ToUnixTimeMilliseconds() & 0xFFFFFFFF))
        };
        
        UpdateParticlesJob updateParticlesJob = new()
        {
            ParticlesRead = particlesAfterForce,  
            ParticlesWrite = particlesFinal,     
            NewVelocities = newVelocities
        };
        
        JobHandle densityHandle = densityJob.Schedule(particleCount, 64);
        JobHandle forceHandle = forceJob.Schedule(particleCount, 64, densityHandle);
        JobHandle integrateHandle = integrateJob.Schedule(particleCount, 64, forceHandle);
        JobHandle updateHandle = updateParticlesJob.Schedule(particleCount, 64, integrateHandle);
        
        updateHandle.Complete();
        EntityCommandBuffer ecb = new(Allocator.TempJob);
        
        for (int i = 0; i < particleCount; i++)
        {
            ecb.SetComponent(allEntities[i], particlesFinal[i]);
            ecb.SetComponent(allEntities[i], transformArray[i]);
        }
        
        ecb.Playback(EntityManager);
        
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