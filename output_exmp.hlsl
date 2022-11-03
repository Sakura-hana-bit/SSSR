//return 1;
//#ifndef TRACE_LIB
//#define TRACE_LIB
//}
//inline half GetScreenFadeBord(half2 pos, half value)
//{
//    half borderDist = min(1 - max(pos.x, pos.y), min(pos.x, pos.y));
//    return saturate(borderDist > value ? 1 : borderDist / value);
//}


//inline float4 SamplScreenRT(int RTindex, float4 uv)
//{
//    float2 bUV = ViewportUVToBufferUV(uv.xy);

//    return SceneTextureLookup(bUV, RTindex, false);
//}

//inline half distanceSquared(half2 A, half2 B)
//{
//    A -= B;
//    return dot(A, A);
//}

//inline half distanceSquared(half3 A, half3 B)
//{
//    A -= B;
//    return dot(A, A);
//}
//void swap(inout half v0, inout half v1)
//{
//    half temp = v0;
//    v0 = v1;
//    v1 = temp;
//}


//float distanceSquared(float2 p1, float2 p0)
//{
//    return distanceSquared(p1.x, p0.x) + distanceSquared(p1.y, p0.y);

//}

//bool intersectsDepthBuffer(half rayZMin, half rayZMax, half sceneZ, half layerThickness)
//{
//    return (rayZMax >= sceneZ - layerThickness) && (rayZMin <= sceneZ);
//}


//void rayIterations(int RTindex, in bool traceBehind_Old, in bool traceBehind, inout half2 P, inout half stepDirection, inout half end, inout int stepCount, inout int maxSteps, inout bool intersecting,
//                   inout half sceneZ, inout half2 dP, inout half3 Q, inout half3 dQ, inout half k, inout half dk,
//                   inout half rayZMin, inout half rayZMax, inout half prevZMaxEstimate, inout bool permute, inout half2 hitPixel,
//                   half2 invSize, inout half layerThickness)
//{
//    bool stop = intersecting;
    
//    for (; (P.x * stepDirection) <= end && stepCount < maxSteps && !stop; P += dP, Q.z += dQ.z, k += dk, stepCount += 1)
//    {
//        rayZMin = prevZMaxEstimate;
//        rayZMax = (dQ.z * 0.5 + Q.z) / (dk * 0.5 + k);
//        prevZMaxEstimate = rayZMax;

//        if (rayZMin > rayZMax)
//        {
//            swap(rayZMin, rayZMax);
//        }

//        hitPixel = permute ? P.yx : P;
//        // sceneZ = tex2Dlod(forntDepth, half4(hitPixel * invSize, 0, 0)).r;
//        //sceneZ = -LinearEyeDepth(sceneZ);
//        sceneZ = SamplScreenRT(RTindex, half4(hitPixel * invSize, 0, 0));
//        bool isBehind = (rayZMin <= sceneZ);

//        if (traceBehind_Old == 1)
//        {
//            intersecting = isBehind && (rayZMax >= sceneZ - layerThickness);
//        }
//        else
//        {
//            intersecting = (rayZMax >= sceneZ - layerThickness);
//        }

//        stop = traceBehind ? intersecting : isBehind;
//    }
//    P -= dP, Q.z -= dQ.z, k -= dk;
//}



//#define UE_TRACE
//bool Linear2D_Trace(int RTindex,
//                             half3 csOrigin,
//                             half3 csDirection,
//                             half2 csZBufferSize,
//                             half jitter,
//                             int maxSteps,
//                             half layerThickness,
//                             half traceDistance,
//                             in out half2 hitPixel,
//                             int stepSize,
//                             bool traceBehind,
//                             in out half3 csHitPoint,
//                             in out half stepCount)
//{

//    half2 invSize = half2(1 / csZBufferSize.x, 1 / csZBufferSize.y);
//    hitPixel = half2(-1, -1);

//    half nearPlaneZ = -0.01;
//    half rayLength = ((csOrigin.z + csDirection.z * traceDistance) > nearPlaneZ) ? ((nearPlaneZ - csOrigin.z) / csDirection.z) : traceDistance;
//    half3 csEndPoint = csDirection * rayLength + csOrigin;
    
//#ifdef UE_TRACE
//    half4 H0 = mul(half4(csOrigin, 1),View.WorldToClip);
//    half4 H1 = mul(half4(csEndPoint, 1),View.WorldToClip);
//#else
//    half4 H0 = mul(projectMatrix, half4(csOrigin, 1));
//    half4 H1 = mul(projectMatrix, half4(csEndPoint, 1));
//#endif
 
//    half k0 = 1 / H0.w;
//    half k1 = 1 / H1.w;
//    half2 P0 = H0.xy * k0;
//    half2 P1 = H1.xy * k1;
//    half3 Q0 = csOrigin * k0;
//    half3 Q1 = csEndPoint * k1;

//#if 1
//    half yMax = csZBufferSize.y - 0.5;
//    half yMin = 0.5;
//    half xMax = csZBufferSize.x - 0.5;
//    half xMin = 0.5;
//    half alpha = 0;

//    if (P1.y > yMax || P1.y < yMin)
//    {
//        half yClip = (P1.y > yMax) ? yMax : yMin;
//        half yAlpha = (P1.y - yClip) / (P1.y - P0.y);
//        alpha = yAlpha;
//    }
//    if (P1.x > xMax || P1.x < xMin)
//    {
//        half xClip = (P1.x > xMax) ? xMax : xMin;
//        half xAlpha = (P1.x - xClip) / (P1.x - P0.x);
//        alpha = max(alpha, xAlpha);
//    }

//    P1 = lerp(P1, P0, alpha);
//    k1 = lerp(k1, k0, alpha);
//    Q1 = lerp(Q1, Q0, alpha);
//#endif

//    P1 = (distanceSquared(P0, P1) < 0.0001) ? P0 + half2(0.01, 0.01) : P1;
//    half2 delta = P1 - P0;
//    bool permute = false;

//    if (abs(delta.x) < abs(delta.y))
//    {
//        permute = true;
//        delta = delta.yx;
//        P1 = P1.yx;
//        P0 = P0.yx;
//    }

//    half stepDirection = sign(delta.x);
//    half invdx = stepDirection / delta.x;
//    half2 dP = half2(stepDirection, invdx * delta.y);
//    half3 dQ = (Q1 - Q0) * invdx;
//    half dk = (k1 - k0) * invdx;
    
//    dP *= stepSize;
//    dQ *= stepSize;
//    dk *= stepSize;
//    P0 += dP * jitter;
//    Q0 += dQ * jitter;
//    k0 += dk * jitter;

//    half3 Q = Q0;
//    half k = k0;
//    half prevZMaxEstimate = csOrigin.z;
//    stepCount = 0;
//    half rayZMax = prevZMaxEstimate, rayZMin = prevZMaxEstimate;
//    half sceneZ = 100000;
//    half end = P1.x * stepDirection;
//    bool intersecting = intersectsDepthBuffer(rayZMin, rayZMax, sceneZ, layerThickness);
//    half2 P = P0;
//    int originalStepCount = 0;

//    bool traceBehind_Old = true;
//    rayIterations(RTindex, traceBehind_Old, traceBehind, P, stepDirection, end, originalStepCount, maxSteps, intersecting, sceneZ, dP, Q, dQ, k, dk, rayZMin, rayZMax, prevZMaxEstimate, permute, hitPixel, invSize, layerThickness);

//    stepCount = originalStepCount;
//    Q.xy += dQ.xy * stepCount;
//    csHitPoint = Q * (1 / k);
//    return intersecting;
//}

//void Linear_2DTrace_SingleSPP( out half4 RayHit_PDF,
//float3 Ray_Origin_VS, float3 ViewNormal, float3 Ray_Dir_VS,float4 H,
//float2 _SSR_ScreenSize, float Jitter, float _SSR_NumSteps_Linear, float _SSR_CullBack,
//float _SSR_Thickness, float _SSR_TraceDistance,
//inout float4 Ray_HitUV, int _SSR_TraceBehind,
//inout float4 Ray_HitPoint, inout float Ray_NumMarch,
//int _SSR_BackwardsRay, inout float Mask, int _SSR_RayStepSize,
//half _SSR_ScreenFade
//)
//{
    
//    //-----Consten Property-------------------------------------------------------------------------
//    half Ray_HitMask = 0.0, Ray_NumMarch = 0.0;
//    half2 Ray_HitUV = 0.0;
//    half3 Ray_HitPoint = 0.0;
    
    
//    //-----BackwardRay-----------------------------------------------------------------------------
   
//    [BRANCH]
//    if (_SSR_BackwardsRay == 0 && Ray_Dir_VS.z < 0)
//    {
//        RayHit_PDF = 0;
//        Mask = 0;
//        return;
//    }

//    //-----Ray Trace-----------------------------------------------------------------------------
//    half Ray_Bump = max(0.01 * Ray_Origin_VS.z, 0.001);
    
//    bool Hit = Linear2D_Trace(1, Ray_Origin_VS + ViewNormal * Ray_Bump, Ray_Dir_VS, _SSR_ScreenSize, Jitter, _SSR_NumSteps_Linear, _SSR_Thickness, _SSR_TraceDistance, Ray_HitUV, _SSR_RayStepSize, _SSR_TraceBehind == 1, Ray_HitPoint, Ray_NumMarch);
//    Ray_HitUV /= _SSR_ScreenSize;

   
//    [BRANCH]
//    if (Hit)
//    {
//        Ray_HitMask = Square(1 - max(2 * half(Ray_NumMarch) / half(_SSR_NumSteps_Linear) - 1, 0));
//        Ray_HitMask *= saturate(((_SSR_TraceDistance - dot(Ray_HitPoint - Ray_Origin_VS, Ray_Dir_VS))));

//        if (_SSR_CullBack < 1)
//        {
//            // half3 Ray_HitNormal_WS = tex2Dlod(_CameraGBufferTexture2, half4(Ray_HitUV, 0, 0)).rgb * 2 - 1;
//            half3 Ray_HitNormal_WS = SamplScreenRT(8, half4(Ray_HitUV, 0, 0)).rgb * 2 - 1;
            
//           // half3 Ray_Dir_WS = mul(_SSR_CameraToWorldMatrix, half4(Ray_Dir_VS, 0)).xyz;
//            half3 Ray_Dir_WS = mul(half4(Ray_Dir_VS, 0), View.ViewToWorld);
            
//            if (dot(Ray_HitNormal_WS, Ray_Dir_WS) > 0)
//                Ray_HitMask = 0;
//        }
//    }
    
//    Mask = Square(Ray_HitMask * GetScreenFadeBord(Ray_HitUV, _SSR_ScreenFade));
//    RayHit_PDF = half4(Ray_HitUV, Mask, H.a);
//   // Mask = Square(Ray_HitMask * GetScreenFadeBord(Ray_HitUV, _SSR_ScreenFade));
//}



//void TRACELIB_MarkFunc()
//{
//return ;

//#endif 
int func()
{
//temp Param
    float Mask;
    half4 RayHit_PDF;
    
//RayMarhing
Linear_2DTrace_SingleSPP(RayHit_PDF,
Ray_Origin_VS, ViewNormalVS, Ray_Dir_VS, H,
 _SSR_ScreenSize, Jitter, _SSR_NumSteps_Linear,
 _SSR_CullBack,
 _SSR_Thickness, _SSR_TraceDistance,
 _SSR_TraceBehind, _SSR_BackwardsRay,
 Mask, _SSR_RayStepSize,
 _SSR_ScreenFade
);
    
//param out put
    Mask_out = Mask;

    RayHit_PDF_out = RayHit_PDF;

    
    return 1;
}

int endFunc()
{
    float4 EndColor;
    
    EndColor = Spatiofilter_SingleSPP
    (
 UV,
 BlueNoise,
 ViewPos.xyz,
 ViewNormal.xyz,
 Roughness.r,
 (int) _SSR_NumResolver,
 _SSR_ScreenSize.xy,
 _SSR_RayCastRT,
  _SSR_RayCastRTSampler,
 _SSR_PosRT,
_SSR_PosRTSampler
);
    
    return EndColor;

}



