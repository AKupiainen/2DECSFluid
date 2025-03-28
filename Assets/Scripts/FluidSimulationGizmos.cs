using Unity.Entities;
using Unity.Mathematics;
using Unity.Transforms;
using UnityEngine;

public class FluidSimulationGizmos : MonoBehaviour
{
    [Header("Visualization Settings")]
    public Color boundaryColor = new Color(0.7f, 0.7f, 0.7f, 0.5f);
    public Color particleColor = new Color(0.2f, 0.6f, 1.0f, 0.8f);
    public bool colorByVelocity = true;
    public bool colorByDensity = false;
    public float particleSize = 0.5f;
    public bool showVelocityVectors = true;
    
    // Cache for accessing ECS entities
    private EntityQuery fluidParticleQuery;
    private EntityManager entityManager;
    private FluidSimulationParameters fluidParams;
    private bool initialized = false;
    
    private void Start()
    {
        // Get references to ECS components
        if (World.DefaultGameObjectInjectionWorld != null)
        {
            entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
            fluidParticleQuery = entityManager.CreateEntityQuery(typeof(FluidParticle), typeof(LocalTransform));
            initialized = true;
        }
    }
    
    private void OnDrawGizmos()
    {
        // Make sure we're initialized
        if (!initialized && Application.isPlaying)
        {
            Start();
            if (!initialized) return;
        }
        
        // Try to get the fluid parameters
        if (!TryGetFluidParameters())
            return;
            
        // Draw the boundary
        DrawBoundary();
        
        // Draw particles if application is playing
        if (Application.isPlaying)
        {
            DrawParticles();
        }
    }
    
    private bool TryGetFluidParameters()
    {
        if (!initialized || !entityManager.IsQueryValid(fluidParticleQuery))
            return false;
            
        // Try to find entity with FluidSimulationParameters
        var paramQuery = entityManager.CreateEntityQuery(typeof(FluidSimulationParameters));
        if (paramQuery.CalculateEntityCount() <= 0)
            return false;
            
        var entity = paramQuery.GetSingletonEntity();
        fluidParams = entityManager.GetComponentData<FluidSimulationParameters>(entity);
        return true;
    }
    
    private void DrawBoundary()
    {
        Gizmos.color = boundaryColor;
        
        // Calculate the boundary corners
        float2 min = fluidParams.BoundaryMin;
        float2 max = fluidParams.BoundaryMax;
        
        // Draw bottom edge
        Gizmos.DrawLine(new Vector3(min.x, min.y, 0), new Vector3(max.x, min.y, 0));
        // Draw right edge
        Gizmos.DrawLine(new Vector3(max.x, min.y, 0), new Vector3(max.x, max.y, 0));
        // Draw top edge
        Gizmos.DrawLine(new Vector3(max.x, max.y, 0), new Vector3(min.x, max.y, 0));
        // Draw left edge
        Gizmos.DrawLine(new Vector3(min.x, max.y, 0), new Vector3(min.x, min.y, 0));
    }
    
    private void DrawParticles()
    {
        if (!initialized || !entityManager.IsQueryValid(fluidParticleQuery))
            return;
            
        // Get all entities with FluidParticle and LocalTransform components
        var entities = fluidParticleQuery.ToEntityArray(Unity.Collections.Allocator.Temp);
        
        foreach (var entity in entities)
        {
            var particle = entityManager.GetComponentData<FluidParticle>(entity);
            var transform = entityManager.GetComponentData<LocalTransform>(entity);
            
            // Determine color based on settings
            Color particleGizmoColor = particleColor;
            
            if (colorByVelocity)
            {
                float speed = math.length(particle.Velocity);
                float normalizedSpeed = math.clamp(speed / 10f, 0, 1);
                particleGizmoColor = Color.Lerp(Color.blue, Color.red, normalizedSpeed);
            }
            else if (colorByDensity)
            {
                float normalizedDensity = math.clamp((particle.Density - 900) / 300f, 0, 1);
                particleGizmoColor = Color.Lerp(Color.green, Color.yellow, normalizedDensity);
            }
            
            Gizmos.color = particleGizmoColor;
            
            // Draw a small sphere for each particle
            float radius = particleSize * 0.5f * fluidParams.SmoothingRadius;
            Gizmos.DrawSphere(transform.Position, radius);
            
            // Optionally draw velocity vectors
            if (showVelocityVectors && math.length(particle.Velocity) > 0.1f)
            {
                Vector3 velocityVector = new Vector3(particle.Velocity.x, particle.Velocity.y, 0) * 0.2f;
                Gizmos.DrawRay(transform.Position, velocityVector);
            }
        }
        
        // Dispose of the temporary array
        entities.Dispose();
    }
    
    // This ensures gizmos are drawn even when the GameObject isn't selected
    private void OnDrawGizmosSelected()
    {
        // Intentionally empty - we handle drawing in OnDrawGizmos
    }
}