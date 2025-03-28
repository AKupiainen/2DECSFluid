using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using Unity.Transforms;
using Random = Unity.Mathematics.Random;

[UpdateInGroup(typeof(InitializationSystemGroup))]
public partial class FluidEmitterSystem : SystemBase
{
    private float _spawnTimer;
    private int _particlesSpawned;
    private bool _isEmitting;
    private Random _random;
    
    protected override void OnCreate()
    {
        base.OnCreate();
        
        RequireForUpdate<FluidSimulationParameters>();
        RequireForUpdate<FluidEmitterParameters>();
        
        _spawnTimer = 0f;
        _particlesSpawned = 0;
        _isEmitting = true;
        _random = Random.CreateFromIndex(1234);
    }

    protected override void OnUpdate()
    {
        FluidSimulationParameters fluidParams = SystemAPI.GetSingleton<FluidSimulationParameters>();
        FluidEmitterParameters emitterParams = SystemAPI.GetSingleton<FluidEmitterParameters>();
        
        if (_particlesSpawned >= emitterParams.MaxParticles && emitterParams.MaxParticles > 0)
        {
            _isEmitting = false;
        }
        
        if (!_isEmitting)
        {
            return;
        }
        
        _spawnTimer += SystemAPI.Time.DeltaTime;
        float timePerParticle = 1.0f / emitterParams.EmissionRate;
        
        if (_spawnTimer < timePerParticle)
        {
            return;
        }
        
        int particlesToEmit = (int)(_spawnTimer / timePerParticle);
        _spawnTimer %= timePerParticle;
        
        if (emitterParams.MaxParticles > 0)
        {
            particlesToEmit = math.min(particlesToEmit, emitterParams.MaxParticles - _particlesSpawned);
            
            if (particlesToEmit <= 0)
            {
                return;
            }
        }
        
        EntityCommandBuffer ecb = new(Allocator.Temp);
        
        for (int i = 0; i < particlesToEmit; i++)
        {
            Entity particle = ecb.CreateEntity();
            
            float2 position = GetEmissionPosition(emitterParams);
            float2 velocity = GetEmissionVelocity(emitterParams);
            
            ecb.AddComponent(particle, new FluidParticle
            {
                Mass = fluidParams.ParticleMass,
                Density = 0.0f,
                Pressure = 0.0f,
                Velocity = velocity,
                Force = float2.zero,
                Viscosity = fluidParams.InitialViscosity
            });
            
            ecb.AddComponent(particle, new LocalTransform
            {
                Position = new float3(position.x, position.y, 0),
                Rotation = quaternion.identity,
                Scale = 1.0f
            });
        }
        
        _particlesSpawned += particlesToEmit;
        
        ecb.Playback(EntityManager);
        ecb.Dispose();
    }
    
    private float2 GetEmissionPosition(FluidEmitterParameters emitterParams)
    {
        switch (emitterParams.EmissionShape)
        {
            case EmissionShape.Point:
                return emitterParams.EmitterPosition;
                
            case EmissionShape.Circle:
                float angle = _random.NextFloat(0, math.PI * 2);
                float radius = _random.NextFloat(0, emitterParams.EmissionRadius);
                return emitterParams.EmitterPosition + new float2(
                    math.cos(angle) * radius,
                    math.sin(angle) * radius
                );
                
            case EmissionShape.Rectangle:
                return emitterParams.EmitterPosition + new float2(
                    _random.NextFloat(-emitterParams.EmissionSize.x, emitterParams.EmissionSize.x) * 0.5f,
                    _random.NextFloat(-emitterParams.EmissionSize.y, emitterParams.EmissionSize.y) * 0.5f
                );
                
            default:
                return emitterParams.EmitterPosition;
        }
    }
    
    private float2 GetEmissionVelocity(FluidEmitterParameters emitterParams)
    {
        float angle;
        
        if (emitterParams.EmissionAngleSpread > 0)
        {
            angle = emitterParams.EmissionDirection + 
                   _random.NextFloat(-emitterParams.EmissionAngleSpread, emitterParams.EmissionAngleSpread) * 0.5f;
        }
        else
        {
            angle = emitterParams.EmissionDirection;
        }
        
        float speed;
        
        if (emitterParams.EmissionSpeedVariation > 0)
        {
            speed = emitterParams.EmissionSpeed * 
                    (1.0f + _random.NextFloat(-emitterParams.EmissionSpeedVariation, emitterParams.EmissionSpeedVariation));
        }
        else
        {
            speed = emitterParams.EmissionSpeed;
        }
        
        return new float2(
            math.cos(angle) * speed,
            math.sin(angle) * speed
        );
    }
}

public enum EmissionShape
{
    Point,
    Circle,
    Rectangle
}