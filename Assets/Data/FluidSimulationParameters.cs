using Unity.Entities;
using Unity.Mathematics;

public struct FluidSimulationParameters : IComponentData
{
    public float DeltaTime;
    public float SmoothingRadius;       
    public float ParticleRadius;        
    public float PressureCoefficient;   
    public float ViscosityCoefficient;  
    public float GravityForce;          
    public float TargetDensity;         
    public float2 BoundaryMin;          
    public float2 BoundaryMax;          
    public float BoundaryDamping;       
    public float OverlapPreventionStrength;
    public float ParticleMass;
    public float InitialViscosity;
    public int ParticleCount;
}