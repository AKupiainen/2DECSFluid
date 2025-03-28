using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

public class FluidSimulationParametersAuthoring : MonoBehaviour
{
    [Header("Simulation Core Settings")]
    [Tooltip("How far particles affect each other")]
    [SerializeField] private float _smoothingRadius = 1.0f;
    
    [Tooltip("Controls the minimum distance between particle centers (prevents overlap)")]
    [SerializeField] private float _particleRadius = 0.2f;
    
    [Tooltip("Target/rest density of the fluid")]
    [SerializeField] private float _targetDensity = 1.0f;
    
    [Header("Force Parameters")]
    [Tooltip("Strength of pressure forces")]
    [SerializeField] private float _pressureCoefficient = 200.0f;
    
    [Tooltip("Strength of viscosity forces")]
    [SerializeField] private float _viscosityCoefficient = 10.0f;
    
    [Tooltip("Strength of gravity")]
    [SerializeField] private float _gravityForce = 9.8f;
    
    [Header("Boundary Settings")]
    [Tooltip("Minimum boundary coordinates")]
    [SerializeField] private Vector2 _boundaryMin = new(-10, -10);
    
    [Tooltip("Maximum boundary coordinates")]
    [SerializeField] private Vector2 _boundaryMax = new(10, 10);
    
    [Range(0, 1)]
    [Tooltip("Damping when hitting boundaries (0 = no bounce, 1 = perfect bounce)")]
    [SerializeField] private float _boundaryDamping = 0.5f;
    
    [Header("Advanced Settings")]
    [Tooltip("Internal fixed timestep for simulation stability")]
    [SerializeField] private float _fixedDeltaTime = 0.005f;
    
    [Range(0.1f, 5.0f)]
    [Tooltip("Overlap prevention strength (higher = stronger repulsion)")]
    [SerializeField] private float _overlapPreventionStrength = 1.0f;

    [Range(0.1f, 5.0f)] [Tooltip("Mass of particle")]
    [SerializeField] private float _particleMass = 1.0f;
    
    [Range(0.1f, 5.0f)] [Tooltip("Mass of particle")]
    [SerializeField] private float _initialViscosity = 0.5f;
    
    [Tooltip("ParticleCount")]
    [SerializeField] private int _particleCount = 1000;
    
    public class FluidSimulationParametersBaker : Baker<FluidSimulationParametersAuthoring>
    {
        public override void Bake(FluidSimulationParametersAuthoring authoring)
        {
            Entity entity = GetEntity(TransformUsageFlags.None);
            
            AddComponent(entity, new FluidSimulationParameters
            {
                DeltaTime = authoring._fixedDeltaTime,
                SmoothingRadius = authoring._smoothingRadius,
                ParticleRadius = authoring._particleRadius,
                PressureCoefficient = authoring._pressureCoefficient,
                ViscosityCoefficient = authoring._viscosityCoefficient,
                GravityForce = authoring._gravityForce,
                TargetDensity = authoring._targetDensity,
                BoundaryMin = new float2(authoring._boundaryMin.x, authoring._boundaryMin.y),
                BoundaryMax = new float2(authoring._boundaryMax.x, authoring._boundaryMax.y),
                BoundaryDamping = authoring._boundaryDamping,
                OverlapPreventionStrength = authoring._overlapPreventionStrength,
                ParticleMass = authoring._particleMass,
                InitialViscosity = authoring._initialViscosity,
                ParticleCount = authoring._particleCount
            });
        }
    }
}