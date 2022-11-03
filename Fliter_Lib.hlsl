return 1;
#ifndef FLITER_LB
#define FLITER_LB
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
inline half Luma4(half3 Color)
{
    return (Color.g * 2) + (Color.r + Color.b);
}
inline half HdrWeight4(half3 Color, half Exposure)
{
    return rcp(Luma4(Color) * Exposure + 4);
}

#define AA_Filter
inline void ResolverAABB(Texture2D currColor, SamplerState _SSR_TemporalScalesampler, half Sharpness, half ExposureScale, half AABBScale, half2 uv, half2 TexelSize, inout half Variance, inout half4 MinColor, inout half4 MaxColor, inout half4 FilterColor)
{
    const int2 SampleOffset[9] = { int2(-1.0, -1.0), int2(0.0, -1.0), int2(1.0, -1.0), int2(-1.0, 0.0), int2(0.0, 0.0), int2(1.0, 0.0), int2(-1.0, 1.0), int2(0.0, 1.0), int2(1.0, 1.0) };
    half4 SampleColors[9];

    for (uint i = 0; i < 9; i++)
    {
#if AA_BicubicFilter
            half4 BicubicSize = half4(TexelSize, 1.0 / TexelSize);
            SampleColors[i] = Texture2DSampleBicubic(currColor, uv + ( SampleOffset[i] / TexelSize), BicubicSize.xy, BicubicSize.zw);
#else
       // SampleColors[i] = tex2D(currColor, uv + (SampleOffset[i] / TexelSize));
        SampleColors[i] = currColor.Sample(_SSR_TemporalScalesampler, uv + (SampleOffset[i] / TexelSize));
#endif
    }

#ifdef AA_Filter
        half SampleWeights[9];
        for(uint j = 0; j < 9; j++) {
            SampleWeights[j] = HdrWeight4(SampleColors[j].rgb, ExposureScale);
        }

        half TotalWeight = 0;
        for(uint k = 0; k < 9; k++) {
            TotalWeight += SampleWeights[k];
        }  
        SampleColors[4] = (SampleColors[0] * SampleWeights[0] + SampleColors[1] * SampleWeights[1] + SampleColors[2] * SampleWeights[2] +  SampleColors[3] * SampleWeights[3] + SampleColors[4] * SampleWeights[4] + SampleColors[5] * SampleWeights[5] +  SampleColors[6] * SampleWeights[6] + SampleColors[7] * SampleWeights[7] + SampleColors[8] * SampleWeights[8]) / TotalWeight;
#endif

    half4 m1 = 0.0;
    half4 m2 = 0.0;
    for (uint x = 0; x < 9; x++)
    {
        m1 += SampleColors[x];
        m2 += SampleColors[x] * SampleColors[x];
    }

    half4 mean = m1 / 9.0;
    half4 stddev = sqrt((m2 / 9.0) - pow2(mean));
        
    MinColor = mean - AABBScale * stddev;
    MaxColor = mean + AABBScale * stddev;

    FilterColor = SampleColors[4];
    MinColor = min(MinColor, FilterColor);
    MaxColor = max(MaxColor, FilterColor);

    half4 TotalVariance = 0;
    for (uint z = 0; z < 9; z++)
    {
        TotalVariance += pow2(Luminance(SampleColors[z].xyz) - Luminance(mean.xyz));
    }
    Variance = saturate((TotalVariance / 9) * 256);
    Variance *= FilterColor.a;
}

float4 Temporalfilter_SingleSPP(float2 Velocity, half2 UV,half3 WorldNormal,
                               float2 _SSR_ScreenSize, float _SSR_TemporalScale, float _SSR_TemporalWeight,
                               Texture2D _SSR_Spatial_RT, SamplerState _SSR_Spatial_RTsampler,
                               Texture2D _SSR_TemporalPrev_RT, SamplerState _SSR_TemporalPrev_RTsampler,
                                
)
{

	/////Get AABB ClipBox
    half SSR_Variance = 0;
    half4 SSR_CurrColor = 0;
    half4 SSR_MinColor, SSR_MaxColor;
    ResolverAABB(_SSR_Spatial_RT, _SSR_Spatial_RTsampler, 0, 10, _SSR_TemporalScale, UV, _SSR_ScreenSize.xy, SSR_Variance, SSR_MinColor, SSR_MaxColor, SSR_CurrColor);

	/////Clamp TemporalColor
    half4 SSR_PrevColor = _SSR_TemporalPrev_RT.Sample(_SSR_TemporalPrev_RTsampler, UV-Velocity);
	//half4 SSR_PrevColor = Bilateralfilter(_SSR_TemporalPrev_RT, UV - Velocity, _SSR_ScreenSize.xy);
    SSR_PrevColor = clamp(SSR_PrevColor, SSR_MinColor, SSR_MaxColor);

	/////Combine TemporalColor
    half Temporal_BlendWeight = saturate(_SSR_TemporalWeight * (1 - length(Velocity) * 8));
    half4 ReflectionColor = lerp(SSR_CurrColor, SSR_PrevColor, Temporal_BlendWeight);

    return ReflectionColor;
}

void Fliter_MarkFUNC()
{
    return;
    
#endif