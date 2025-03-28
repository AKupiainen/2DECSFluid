using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

public class FluidEmitterAuthoring : MonoBehaviour
{
    [Header("Emission Shape")]
    [SerializeField] private EmissionShape _emissionShape = EmissionShape.Circle;
    [SerializeField] private float _emissionRadius = 0.5f;
    [SerializeField] private Vector2 _emissionSize = new(1f, 1f);

    [Header("Emission Rate")]
    [SerializeField] private float _emissionRate = 100f;
    [SerializeField] private int _maxParticles = 1000; 

    [Header("Emission Direction")]
    [SerializeField] private float _emissionDirection;
    [SerializeField] private float _emissionAngleSpread = 30f; 

    [Header("Emission Speed")]
    [SerializeField] private float _emissionSpeed = 2f;
    [SerializeField][Range(0f, 1f)] private float _emissionSpeedVariation = 0.2f;

    private class FluidEmitterBaker : Baker<FluidEmitterAuthoring>
    {
        public override void Bake(FluidEmitterAuthoring authoring)
        {
            Entity entity = GetEntity(TransformUsageFlags.Dynamic);
            AddComponent(entity, new FluidEmitterParameters
            {
                EmitterPosition = new float2(authoring.transform.position.x, authoring.transform.position.y),
                EmissionShape = authoring._emissionShape,
                EmissionRadius = authoring._emissionRadius,
                EmissionSize = new float2(authoring._emissionSize.x, authoring._emissionSize.y),
                EmissionRate = authoring._emissionRate,
                MaxParticles = authoring._maxParticles,
                EmissionDirection = math.radians(authoring._emissionDirection),
                EmissionAngleSpread = math.radians(authoring._emissionAngleSpread),
                EmissionSpeed = authoring._emissionSpeed,
                EmissionSpeedVariation = authoring._emissionSpeedVariation
            });
        }
    }
}