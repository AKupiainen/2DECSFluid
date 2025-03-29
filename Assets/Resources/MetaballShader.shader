Shader "Custom/MetaballShader"
{
    Properties
    {
        _DensityTex ("Density Texture", 2D) = "white" {}
        _WaterColor ("Water Color", Color) = (0.2, 0.5, 0.9, 0.4)
        _HighlightColor ("Highlight Color", Color) = (0.7, 0.9, 1.0, 0.6)
        _DensityThreshold ("Density Threshold", Range(0.001, 1.0)) = 0.15
        _WaterSpeed ("Water Speed", Range(0.1, 5.0)) = 1.2
        _WaveScale ("Wave Scale", Range(1.0, 20.0)) = 6.0
        _CartoonEdge ("Cartoon Edge Threshold", Range(0.01, 0.5)) = 0.1
        _CausticIntensity ("Caustic Intensity", Range(0.0, 2.0)) = 0.5
        _Transparency ("Transparency", Range(0.0, 1.0)) = 0.5
        _EdgeThickness ("Edge Thickness", Range(0.0, 0.1)) = 0.02
        _EdgeSoftness ("Edge Softness", Range(0.0, 0.5)) = 0.05
    }
    
    SubShader
    {
        Tags { "Queue"="Transparent" "RenderType"="Transparent" }
        LOD 100
        
        Blend SrcAlpha OneMinusSrcAlpha
        ZWrite Off
        
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            
            #include "UnityCG.cginc"
            
            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
                float3 normal : NORMAL;
            };
            
            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
                float3 viewDir : TEXCOORD1;
                float3 worldNormal : TEXCOORD2;
                float4 worldPos : TEXCOORD3;
            };
            
            sampler2D _DensityTex;
            float4 _WaterColor;
            float4 _HighlightColor;
            float _DensityThreshold;
            float _WaterSpeed;
            float _WaveScale;
            float _CartoonEdge;
            float _CausticIntensity;
            float _Transparency;
            float _EdgeThickness;
            float _EdgeSoftness;
            float _DebugView;
            
            float noise(float2 p) {
                return frac(sin(dot(p, float2(12.9898, 78.233))) * 43758.5453);
            }
            
            float cartoonWave(float2 uv, float time) {
              
                float wave1 = sin(uv.x * _WaveScale + time * 0.5) * 0.5 + 0.5;
                float wave2 = sin(uv.y * _WaveScale * 0.7 - time * 0.4) * 0.5 + 0.5;
                float wave3 = sin((uv.x + uv.y) * _WaveScale * 0.5 + time * 0.3) * 0.5 + 0.5;
                
                float wave = (wave1 + wave2 + wave3) / 3.0;
                
                return smoothstep(0.4 - _CartoonEdge, 0.4 + _CartoonEdge, wave);
            }
            
            float cartoonCaustic(float2 uv, float time) {
                
                float c1 = sin(uv.x * 10.0 + time) * sin(uv.y * 10.0 - time * 0.7);
                float c2 = sin((uv.x + uv.y) * 8.0 - time * 0.5) * cos((uv.x - uv.y) * 8.0 + time * 0.3);
                
                float caustic = (c1 + c2) * 0.5 + 0.5;
                return smoothstep(0.5 - _CartoonEdge, 0.5 + _CartoonEdge, caustic) * _CausticIntensity;
            }
            
            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;
                o.viewDir = normalize(WorldSpaceViewDir(v.vertex));
                o.worldNormal = UnityObjectToWorldNormal(v.normal);
                o.worldPos = mul(unity_ObjectToWorld, v.vertex);
                return o;
            }
            
            fixed4 frag (v2f i) : SV_Target
            {
                float time = _Time.y * _WaterSpeed;
                float density = tex2D(_DensityTex, i.uv).r;
                
                float edgeFactor = smoothstep(_DensityThreshold - _EdgeSoftness, _DensityThreshold + _EdgeSoftness, density);
                if (edgeFactor <= 0.01)
                    discard; 
                
                float wavePattern = cartoonWave(i.uv, time);
                float causticPattern = cartoonCaustic(i.uv, time);
                
                float viewFresnel = 1.0 - saturate(dot(i.worldNormal, i.viewDir));
                float edge = smoothstep(0.6, 0.8, viewFresnel);
                
                float densityEdge = abs(density - _DensityThreshold) < _EdgeThickness ? 1.0 : 0.0;
                
                float3 waterColor = lerp(_WaterColor.rgb, _HighlightColor.rgb, wavePattern);
                
                waterColor += _HighlightColor.rgb * causticPattern;
                waterColor = lerp(waterColor, _HighlightColor.rgb, edge * 0.7);
    
                waterColor = lerp(waterColor, _HighlightColor.rgb, densityEdge * 0.5);
                
                float alpha = lerp(_WaterColor.a * _Transparency, _HighlightColor.a, 
                                 edge * 0.5 + densityEdge * 0.5) * edgeFactor;
                
                return float4(waterColor, alpha);
            }
            ENDCG
        }
    }
}