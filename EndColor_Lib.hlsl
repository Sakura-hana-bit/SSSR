
return 1;
#ifndef EndColor_LIB
#define EndColor_LIB
}

inline float4 SamplScreenRT(int RTindex, float4 uv)
{
    float2 bUV = ViewportUVToBufferUV(uv.xy);

    return SceneTextureLookup(bUV, RTindex, false);
}
float pow2(float x)
{
    return x * x;
}

float2 pow2(float2 x)
{
    return x * x;
}

float3 pow2(float3 x)
{
    return x * x;
}

float4 pow2(float4 x)
{
    return x * x;
}

float pow3(float x)
{
    return x * x * x;
}

float2 pow3(float2 x)
{
    return x * x * x;
}

float3 pow3(float3 x)
{
    return x * x * x;
}

float4 pow3(float4 x)
{
    return x * x * x;
}

float pow4(float x)
{
    float xx = x * x;
    return xx * xx;
}

float2 pow4(float2 x)
{
    float2 xx = x * x;
    return xx * xx;
}

float3 pow4(float3 x)
{
    float3 xx = x * x;
    return xx * xx;
}

float4 pow4(float4 x)
{
    float4 xx = x * x;
    return xx * xx;
}
float D_GGX_FORUE(float NoH, float Roughness)
{
    Roughness = pow4(Roughness);
    float D = (NoH * Roughness - NoH) * NoH + 1;
    return Roughness / (3.1415926 * pow2(D));
}
float Vis_SmithGGXCorrelated(float NoL, float NoV, float Roughness)
{
    float a = pow2(Roughness);
    float LambdaV = NoV * sqrt((-NoL * a + NoL) * NoL + a);
    float LambdaL = NoL * sqrt((-NoV * a + NoV) * NoV + a);
    return (0.5 / (LambdaL + LambdaV)) / PI;
}
float SSR_BRDF(float3 V, float3 L, float3 N, float Roughness)
{
    float3 H = normalize(L + V);

    float NoH = max(dot(N, H), 0);
    float NoL = max(dot(N, L), 0);
    float NoV = max(dot(N, V), 0);

    float D = D_GGX_FORUE(NoH, Roughness);
    float G = Vis_SmithGGXCorrelated(NoL, NoV, Roughness);

    return max(0, D * G);
}

float4 test(float2 Neighbor_UV, texture2D _SSR_RayCastRT, SamplerState _SSR_RayCastRTSampler)
{
    return _SSR_RayCastRT.Sample(_SSR_RayCastRTSampler, Neighbor_UV);
}

static const int2 offset[9] = { int2(-2.0, -2.0), int2(0.0, -2.0), int2(2.0, -2.0), int2(-2.0, 0.0), int2(0.0, 0.0), int2(2.0, 0.0), int2(-2.0, 2.0), int2(0.0, 2.0), int2(2.0, 2.0) };
float4 Spatiofilter_SingleSPP
    (
    float2 UV,
                              float2 BlueNoise,
                              float3 ViewPos, //in camera space is mean View Dir
                              float3 ViewNormal,
                              float Roughness,

                              int _SSR_NumResolver,
                              float2 _SSR_ScreenSize,

                              texture2D _SSR_RayCastRT,
                              SamplerState _SSR_RayCastRTSampler,
                              texture2D _SSR_PosRT,
                              SamplerState _SSR_PosRTSampler
)
{
    float2x2 OffsetRotationMatrix = float2x2(BlueNoise.x, BlueNoise.y, -BlueNoise.y, -BlueNoise.x);

    float NumWeight = 0;
    float Weight = 0;
    float2 Offset_UV, Neighbor_UV;
    float4 SampleColor = 0, ReflecttionColor = 0;
    float4 HitUV_PDF = 0;
    float4 d = 0;
    for (int i = 0; i < _SSR_NumResolver; i++)
    {
        Offset_UV = mul(OffsetRotationMatrix, offset[i] * (1 / _SSR_ScreenSize.xy));
        Neighbor_UV = UV + Offset_UV;

         HitUV_PDF = _SSR_RayCastRT.Sample(_SSR_RayCastRTSampler, Neighbor_UV);
        float3 Hit_ViewPos = _SSR_PosRT.Sample(_SSR_PosRTSampler, Neighbor_UV, 1);

		///SpatioSampler
        
        Weight = SSR_BRDF(normalize(-ViewPos), normalize(Hit_ViewPos - ViewPos), ViewNormal, Roughness) / max(1e-5, HitUV_PDF.a);
       // SampleColor.rgb = tex2Dlod(_SSR_SceneColor_RT, float4(HitUV_PDF.rg, 0, 0)).rgb;
        SampleColor.rgb = SamplScreenRT(0, float4(HitUV_PDF.rg, 0, 0)).rgb;
        d = float4(SampleColor.rgb, 0);
        SampleColor.rgb /= 1 + Luminance(SampleColor.rgb);
        SampleColor.a = HitUV_PDF.b;

        ReflecttionColor += SampleColor * Weight;
        NumWeight += Weight;
    }

    ReflecttionColor /= NumWeight;
    ReflecttionColor.rgb /= 1 - Luminance(ReflecttionColor.rgb);
    ReflecttionColor = max(1e-5, ReflecttionColor);
	//ReflecttionColor.a = tex2D(_SSR_RayMask_RT, UV).r;
    return float4(ReflecttionColor);
}

void EndColor_MarkFunc()
{
return;
#endif