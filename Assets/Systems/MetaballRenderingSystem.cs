using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using Unity.Transforms;
using UnityEngine;

[UpdateInGroup(typeof(PresentationSystemGroup))]
public partial class MetaballRenderingSystem : SystemBase
{
    private static readonly int DensityTex = Shader.PropertyToID("_DensityTex");
    private static readonly int FluidColor = Shader.PropertyToID("_FluidColor");
    private static readonly int DensityThreshold = Shader.PropertyToID("_DensityThreshold");
    private static readonly int Smoothness = Shader.PropertyToID("_Smoothness");
    private static readonly int FresnelPower = Shader.PropertyToID("_FresnelPower");
    private static readonly int Particles = Shader.PropertyToID("Particles");
    private static readonly int Result = Shader.PropertyToID("Result");
    private static readonly int ParticleCount = Shader.PropertyToID("ParticleCount");
    private static readonly int SmoothingRadius = Shader.PropertyToID("SmoothingRadius");
    private static readonly int ParticleRadius = Shader.PropertyToID("ParticleRadius");
    private static readonly int DensityMultiplier = Shader.PropertyToID("DensityMultiplier");
    private static readonly int BoundaryMin = Shader.PropertyToID("BoundaryMin");
    private static readonly int BoundaryMax = Shader.PropertyToID("BoundaryMax");
    private static readonly int BoundaryCenter = Shader.PropertyToID("BoundaryCenter");
    
    private EntityQuery _fluidParticleQuery;
    private Material _metaballMaterial;
    private RenderTexture _metaballTexture;
    private ComputeShader _metaballCompute;
    private int _metaballKernel;
    private GameObject _metaballQuad;
    private MeshRenderer _quadRenderer;
    
    protected override void OnCreate()
    {
        _fluidParticleQuery = GetEntityQuery(typeof(FluidParticle), typeof(LocalTransform));
        RequireForUpdate<FluidSimulationParameters>();
    }
    
    protected override void OnStartRunning()
    {
        _metaballCompute = Resources.Load<ComputeShader>("MetaballCompute");
        _metaballKernel = _metaballCompute.FindKernel("CSMain");
        _metaballMaterial = new Material(Shader.Find("Custom/MetaballShader"));
        
        CreateMetaballQuad();
        InitializeRenderTexture();
    }
    
    private void CreateMetaballQuad()
    {
        FluidSimulationParameters fluidParams = SystemAPI.GetSingleton<FluidSimulationParameters>();
        
        _metaballQuad = GameObject.CreatePrimitive(PrimitiveType.Quad);
        _metaballQuad.name = "MetaballRenderQuad";
        
        float2 boundarySize = fluidParams.BoundaryMax - fluidParams.BoundaryMin;
        float2 boundaryCenter = (fluidParams.BoundaryMax + fluidParams.BoundaryMin) * 0.5f;
        
        _metaballQuad.transform.position = new Vector3(boundaryCenter.x, boundaryCenter.y, 0);
        _metaballQuad.transform.localScale = new Vector3(boundarySize.x, boundarySize.y, 1);
        
        _quadRenderer = _metaballQuad.GetComponent<MeshRenderer>();
        _quadRenderer.material = _metaballMaterial;
        
        if (Camera.main == null)
        {
            GameObject cameraObj = new GameObject("MetaballCamera");
            Camera cam = cameraObj.AddComponent<Camera>();
            cam.orthographic = true;
            
            float maxDimension = Mathf.Max(boundarySize.x, boundarySize.y);
            cam.orthographicSize = maxDimension * 0.6f;
            
            cameraObj.transform.position = new Vector3(boundaryCenter.x, boundaryCenter.y, -10);
            cameraObj.transform.LookAt(_metaballQuad.transform);
        }
    }
    
    private void InitializeRenderTexture()
    {
        FluidSimulationParameters fluidParams = SystemAPI.GetSingleton<FluidSimulationParameters>();
        int resolution = fluidParams.MetaballResolution;
        
        if (_metaballTexture != null)
        {
            _metaballTexture.Release();
        }
        
        _metaballTexture = new RenderTexture(resolution, resolution, 0, RenderTextureFormat.RFloat)
        {
            enableRandomWrite = true,
            filterMode = FilterMode.Bilinear
        };
        
        _metaballTexture.Create();
        _metaballMaterial.SetTexture(DensityTex, _metaballTexture);
    }
    
    protected override void OnUpdate()
    {
        if (_metaballCompute == null || _metaballMaterial == null || _quadRenderer == null)
        {
            return;
        }

        if (_fluidParticleQuery.IsEmpty)
        {
            return;
        }

        FluidSimulationParameters fluidParams = SystemAPI.GetSingleton<FluidSimulationParameters>();
        
        if (_metaballTexture.width != fluidParams.MetaballResolution)
        {
            InitializeRenderTexture();
        }

        _metaballMaterial.SetColor(FluidColor, new Color(
            fluidParams.FluidColor.x, 
            fluidParams.FluidColor.y, 
            fluidParams.FluidColor.z, 
            fluidParams.FluidColor.w));
        _metaballMaterial.SetFloat(DensityThreshold, fluidParams.DensityThreshold);
        _metaballMaterial.SetFloat(Smoothness, fluidParams.Smoothness);
        _metaballMaterial.SetFloat(FresnelPower, fluidParams.FresnelPower);
        
        NativeList<float3> particlePositions = new(Allocator.TempJob);
        
        Entities
            .WithStoreEntityQueryInField(ref _fluidParticleQuery)
            .ForEach((in FluidParticle _, in LocalTransform transform) => {
                particlePositions.Add(transform.Position);
            }).Run();
        
        int particleCount = particlePositions.Length;
        
        if (particleCount == 0)
        {
            particlePositions.Dispose();
            return;
        }
        
        ComputeBuffer particleBuffer = new(particleCount, sizeof(float) * 3);
        particleBuffer.SetData(particlePositions.AsArray());
        
        _metaballCompute.SetBuffer(_metaballKernel, Particles, particleBuffer);
        _metaballCompute.SetTexture(_metaballKernel, Result, _metaballTexture);
        _metaballCompute.SetInt(ParticleCount, particleCount);
        _metaballCompute.SetFloat(SmoothingRadius, fluidParams.SmoothingRadius);
        _metaballCompute.SetFloat(ParticleRadius, fluidParams.ParticleRadius);
        _metaballCompute.SetFloat(DensityMultiplier, 50.0f);
        
        float2 boundaryCenter = (fluidParams.BoundaryMax + fluidParams.BoundaryMin) * 0.5f;
        _metaballCompute.SetFloats(BoundaryMin, fluidParams.BoundaryMin.x, fluidParams.BoundaryMin.y);
        _metaballCompute.SetFloats(BoundaryMax, fluidParams.BoundaryMax.x, fluidParams.BoundaryMax.y);
        _metaballCompute.SetFloats(BoundaryCenter, boundaryCenter.x, boundaryCenter.y);
        
        int resolution = fluidParams.MetaballResolution;
        int threadGroupsX = Mathf.CeilToInt(resolution / 8.0f);
        int threadGroupsY = Mathf.CeilToInt(resolution / 8.0f);
        
        _metaballCompute.Dispatch(_metaballKernel, threadGroupsX, threadGroupsY, 1);
        
        particleBuffer.Release();
        particlePositions.Dispose();
    }
    
    protected override void OnDestroy()
    {
        if (_metaballTexture != null)
        {
            _metaballTexture.Release();
            _metaballTexture = null;
        }
        
        if (_metaballMaterial != null)
        {
            Object.Destroy(_metaballMaterial);
        }
        
        if (_metaballQuad != null)
        {
            Object.Destroy(_metaballQuad);
        }
    }
}