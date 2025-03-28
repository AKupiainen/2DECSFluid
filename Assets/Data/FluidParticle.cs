using Unity.Entities;
using Unity.Mathematics;

public struct FluidParticle : IComponentData
{
    public float Mass;
    public float Density;
    public float Pressure;
    public float2 Velocity;
    public float2 Force;
    public float Viscosity;
}