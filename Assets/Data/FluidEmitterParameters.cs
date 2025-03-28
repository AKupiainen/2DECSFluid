using Unity.Entities;
using Unity.Mathematics;

public struct FluidEmitterParameters : IComponentData
{
    public float2 EmitterPosition;
    public EmissionShape EmissionShape;
    public float EmissionRadius;      
    public float2 EmissionSize;      
    
    public float EmissionRate;        
    public int MaxParticles;        
    
    public float EmissionDirection;   
    public float EmissionAngleSpread; 
    public float EmissionSpeed;
    public float EmissionSpeedVariation;
}