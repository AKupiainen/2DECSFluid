using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using Unity.Transforms;
using UnityEngine;

public class FluidSimulationGizmos : MonoBehaviour
{
    [Header("Visualization Settings")]
    [SerializeField] private Color _boundaryColor = new(0.7f, 0.7f, 0.7f, 0.5f);
    [SerializeField] private Color _particleColor = new Color(0.2f, 0.6f, 1.0f, 0.8f);
    [SerializeField] private bool _colorByVelocity = true;
    [SerializeField] private bool _colorByDensity;
    [SerializeField] private float _particleSize = 0.5f;
    [SerializeField] private bool _showVelocityVectors = true;
    
    private EntityQuery _fluidParticleQuery;
    private EntityManager _entityManager;
    private FluidSimulationParameters _fluidParams;
    private bool _initialized;
    
    private void Start()
    {
        if (World.DefaultGameObjectInjectionWorld != null)
        {
            _entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
            _fluidParticleQuery = _entityManager.CreateEntityQuery(typeof(FluidParticle), typeof(LocalTransform));
            _initialized = true;
        }
    }
    
    private void OnDrawGizmos()
    {
        if (!_initialized && Application.isPlaying)
        {
            Start();
            
            if (!_initialized)
            {
                return;
            }
        }
        
        if (!TryGetFluidParameters())
        {
            return;
        }

        DrawBoundary();
        
        if (Application.isPlaying)
        {
            DrawParticles();
        }
    }
    
    private bool TryGetFluidParameters()
    {
        if (!_initialized || !_entityManager.IsQueryValid(_fluidParticleQuery))
        {
            return false;
        }

        EntityQuery paramQuery = _entityManager.CreateEntityQuery(typeof(FluidSimulationParameters));
        
        
        if (paramQuery.CalculateEntityCount() <= 0)
        {
            return false;
        }

        Entity entity = paramQuery.GetSingletonEntity();
        _fluidParams = _entityManager.GetComponentData<FluidSimulationParameters>(entity);
        return true;
    }
    
    private void DrawBoundary()
    {
        Gizmos.color = _boundaryColor;
        
        float2 min = _fluidParams.BoundaryMin;
        float2 max = _fluidParams.BoundaryMax;
        
        Gizmos.DrawLine(new Vector3(min.x, min.y, 0), new Vector3(max.x, min.y, 0));
        Gizmos.DrawLine(new Vector3(max.x, min.y, 0), new Vector3(max.x, max.y, 0));
        Gizmos.DrawLine(new Vector3(max.x, max.y, 0), new Vector3(min.x, max.y, 0));
        Gizmos.DrawLine(new Vector3(min.x, max.y, 0), new Vector3(min.x, min.y, 0));
    }
    
    private void DrawParticles()
    {
        if (!_initialized || !_entityManager.IsQueryValid(_fluidParticleQuery))
        {
            return;
        }
        
        NativeArray<Entity> entities = _fluidParticleQuery.ToEntityArray(Allocator.Temp);
        
        foreach (Entity entity in entities)
        {
            FluidParticle particle = _entityManager.GetComponentData<FluidParticle>(entity);
            LocalTransform transform = _entityManager.GetComponentData<LocalTransform>(entity);
            
            Color particleGizmoColor = _particleColor;
            
            if (_colorByVelocity)
            {
                float speed = math.length(particle.Velocity);
                float normalizedSpeed = math.clamp(speed / 10f, 0, 1);
                particleGizmoColor = Color.Lerp(Color.blue, Color.red, normalizedSpeed);
            }
            else if (_colorByDensity)
            {
                float normalizedDensity = math.clamp((particle.Density - 900) / 300f, 0, 1);
                particleGizmoColor = Color.Lerp(Color.green, Color.yellow, normalizedDensity);
            }
            
            Gizmos.color = particleGizmoColor;
            
            float radius = _particleSize * 0.5f * _fluidParams.SmoothingRadius;
            Gizmos.DrawSphere(transform.Position, radius);
            
            if (_showVelocityVectors && math.length(particle.Velocity) > 0.1f)
            {
                Vector3 velocityVector = new Vector3(particle.Velocity.x, particle.Velocity.y, 0) * 0.2f;
                Gizmos.DrawRay(transform.Position, velocityVector);
            }
        }
        
        entities.Dispose();
    }
}