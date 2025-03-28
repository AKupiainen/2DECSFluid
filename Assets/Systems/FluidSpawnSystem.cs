using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using Unity.Transforms;
using UnityEngine;

[UpdateInGroup(typeof(InitializationSystemGroup))]
public partial class FluidSpawnSystem : SystemBase
{
    private EntityQuery fluidParticleQuery;
    private float spawnTimer;
    private int particlesSpawned;
    private bool isSpawning;
    
    protected override void OnCreate()
    {
        base.OnCreate();
        fluidParticleQuery = GetEntityQuery(typeof(FluidParticle));
        RequireForUpdate<FluidSimulationParameters>(); // Fix for the previous error
        spawnTimer = 0f;
        particlesSpawned = 0;
        isSpawning = false;
    }

    protected override void OnUpdate()
    {
        var fluidParams = SystemAPI.GetSingleton<FluidSimulationParameters>();
        
        // Check if we need to start spawning
        if (!isSpawning && fluidParticleQuery.IsEmpty)
        {
            isSpawning = true;
        }
        
        // If not spawning or all particles spawned, exit
        if (!isSpawning || particlesSpawned >= fluidParams.ParticleCount)
        {
            return;
        }
        
        // Increment timer
        spawnTimer += SystemAPI.Time.DeltaTime;
        
        // Spawn rate: particles per second
        float spawnRate = 100f; // Adjust this value to control spawn speed
        
        // Number of particles to spawn this frame
        int particlesToSpawn = (int)(spawnTimer * spawnRate);
        if (particlesToSpawn <= 0)
            return;
        
        // Reset timer
        spawnTimer = 0f;
        
        // Calculate grid dimensions
        int particlesPerSide = (int)math.sqrt(fluidParams.ParticleCount);
        float spacing = fluidParams.SmoothingRadius * 0.5f;
        float startX = (fluidParams.BoundaryMin.x + fluidParams.BoundaryMax.x) * 0.5f - (particlesPerSide * spacing * 0.5f);
        float startY = (fluidParams.BoundaryMin.y + fluidParams.BoundaryMax.y) * 0.5f - (particlesPerSide * spacing * 0.5f);
        
        // Limit particles to spawn to remaining count
        particlesToSpawn = math.min(particlesToSpawn, fluidParams.ParticleCount - particlesSpawned);
        
        EntityCommandBuffer ecb = new EntityCommandBuffer(Allocator.Temp);
        
        for (int i = 0; i < particlesToSpawn; i++)
        {
            // Calculate the position in the grid
            int totalIndex = particlesSpawned + i;
            int x = totalIndex % particlesPerSide;
            int y = totalIndex / particlesPerSide;
            
            // Skip if we've gone beyond the grid bounds
            if (y >= particlesPerSide)
                continue;
            
            Entity particle = ecb.CreateEntity();
            
            ecb.AddComponent(particle, new FluidParticle
            {
                Mass = fluidParams.ParticleMass,
                Density = 0.0f,
                Pressure = 0.0f,
                Velocity = float2.zero,
                Force = float2.zero,
                Viscosity = fluidParams.InitialViscosity
            });
            
            float posX = startX + x * spacing;
            float posY = startY + y * spacing;
            
            ecb.AddComponent(particle, new LocalTransform
            {
                Position = new float3(posX, posY, 0),
                Rotation = quaternion.identity,
                Scale = 1.0f
            });
        }
        
        particlesSpawned += particlesToSpawn;
        
        ecb.Playback(EntityManager);
        ecb.Dispose();
    }
}