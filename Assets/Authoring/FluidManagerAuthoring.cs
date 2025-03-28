using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

// This is the authoring component that appears in the Inspector
public class FluidSimulationParametersAuthoring : MonoBehaviour
{
    [Header("Simulation Core Settings")]
    [Tooltip("How far particles affect each other")]
    public float smoothingRadius = 1.0f;
    
    [Tooltip("Controls the minimum distance between particle centers (prevents overlap)")]
    public float particleRadius = 0.2f;
    
    [Tooltip("Target/rest density of the fluid")]
    public float targetDensity = 1.0f;
    
    [Header("Force Parameters")]
    [Tooltip("Strength of pressure forces")]
    public float pressureCoefficient = 200.0f;
    
    [Tooltip("Strength of viscosity forces")]
    public float viscosityCoefficient = 10.0f;
    
    [Tooltip("Strength of gravity")]
    public float gravityForce = 9.8f;
    
    [Header("Boundary Settings")]
    [Tooltip("Minimum boundary coordinates")]
    public Vector2 boundaryMin = new(-10, -10);
    
    [Tooltip("Maximum boundary coordinates")]
    public Vector2 boundaryMax = new(10, 10);
    
    [Range(0, 1)]
    [Tooltip("Damping when hitting boundaries (0 = no bounce, 1 = perfect bounce)")]
    public float boundaryDamping = 0.5f;
    
    [Header("Advanced Settings")]
    [Tooltip("Internal fixed timestep for simulation stability")]
    public float fixedDeltaTime = 0.005f;
    
    [Range(0.1f, 5.0f)]
    [Tooltip("Overlap prevention strength (higher = stronger repulsion)")]
    public float overlapPreventionStrength = 1.0f;

    [Range(0.1f, 5.0f)] [Tooltip("Mass of particle")]
    public float particleMass = 1.0f;
    
    [Range(0.1f, 5.0f)] [Tooltip("Mass of particle")]
    public float initialViscosity = 0.5f;
    
    [Tooltip("ParticleCount")]
    public int particleCount = 1000;
    
    // This class will handle the conversion to a runtime ECS component
    public class FluidSimulationParametersBaker : Baker<FluidSimulationParametersAuthoring>
    {
        public override void Bake(FluidSimulationParametersAuthoring authoring)
        {
            var entity = GetEntity(TransformUsageFlags.None);
            
            AddComponent(entity, new FluidSimulationParameters
            {
                DeltaTime = authoring.fixedDeltaTime,
                SmoothingRadius = authoring.smoothingRadius,
                ParticleRadius = authoring.particleRadius,
                PressureCoefficient = authoring.pressureCoefficient,
                ViscosityCoefficient = authoring.viscosityCoefficient,
                GravityForce = authoring.gravityForce,
                TargetDensity = authoring.targetDensity,
                BoundaryMin = new float2(authoring.boundaryMin.x, authoring.boundaryMin.y),
                BoundaryMax = new float2(authoring.boundaryMax.x, authoring.boundaryMax.y),
                BoundaryDamping = authoring.boundaryDamping,
                OverlapPreventionStrength = authoring.overlapPreventionStrength,
                ParticleMass = authoring.particleMass,
                InitialViscosity = authoring.initialViscosity,
                ParticleCount = authoring.particleCount
            });
        }
    }
}