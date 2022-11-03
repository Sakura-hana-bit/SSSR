return 1;
#ifndef TRACE_LIB
#define TRACE_LIB
}
inline half GetScreenFadeBord(half2 pos, half value)
{
    half borderDist = min(1 - max(pos.x, pos.y), min(pos.x, pos.y));
    return saturate(borderDist > value ? 1 : borderDist / value);
}


inline float4 SamplScreenRT(int RTindex, float4 uv)
{
    float2 bUV = ViewportUVToBufferUV(uv.xy);

    return SceneTextureLookup(bUV, RTindex, false);
}

inline half distanceSquared(half2 A, half2 B)
{
    A -= B;
    return dot(A, A);
}

inline half distanceSquared(half3 A, half3 B)
{
    A -= B;
    return dot(A, A);
}
void swap(inout half v0, inout half v1)
{
    half temp = v0;
    v0 = v1;
    v1 = temp;
}




bool intersectsDepthBuffer(half rayZMin, half rayZMax, half sceneZ, half layerThickness)
{
    return (rayZMax >= sceneZ - layerThickness); // && (rayZMin <= sceneZ);
}
float2 ScreenPosToUv10(float2 Pos)
{
    return float2(Pos.x * 0.5 + 0.5, Pos.y * -0.5 + 0.5);

}
bool rayIterationsForUE(int RTindex, in bool traceBehind_Old, in bool traceBehind, inout half2 P, inout half stepDirection, inout half end, inout int stepCount, inout int maxSteps, inout bool intersecting,
                   inout half sceneZ, inout half2 dP, inout half3 Q, inout half3 dQ, inout half k, inout half dk,
                   inout half rayZMin, inout half rayZMax, inout half prevZMaxEstimate, inout bool permute, inout half2 hitPixel,
                   half2 invSize, inout half layerThickness)
{
    bool stop = intersecting;
    stepCount++;
   // maxSteps
    for (; (P.x * stepDirection) <= end && stepCount < maxSteps && !stop; P += dP, Q.z += dQ.z, k += dk, stepCount += 1)
    {
        rayZMin = prevZMaxEstimate;
        rayZMax = (dQ.z * 0.5 + Q.z) / (dk * 0.5 + k);
        prevZMaxEstimate = rayZMax;

        if (rayZMin > rayZMax)
        {
            swap(rayZMin, rayZMax);
        }

        hitPixel = permute ? P.yx : P;
        
        // sceneZ = tex2Dlod(forntDepth, half4(hitPixel * invSize, 0, 0)).r;
        //sceneZ = -LinearEyeDepth(sceneZ);
        
        sceneZ = SamplScreenRT(RTindex, half4(ScreenPosToUv10(hitPixel * invSize), 0, 0)).r;
        bool isBehind = (rayZMin <= sceneZ);

        if (traceBehind_Old == 1)
        {
            intersecting = isBehind && (rayZMax >= sceneZ - layerThickness);
        }
        else
        {
            intersecting = (rayZMax >= sceneZ - layerThickness);
        }

        stop = traceBehind ? intersecting : isBehind;
    }
    P -= dP, Q.z -= dQ.z, k -= dk;
    return stop;
}

bool rayIterations(int RTindex, in bool traceBehind_Old, in bool traceBehind, inout half2 P, inout half stepDirection, inout half end, inout int stepCount, inout int maxSteps, inout bool intersecting,
                   inout half sceneZ, inout half2 dP, inout half3 Q, inout half3 dQ, inout half k, inout half dk,
                   inout half rayZMin, inout half rayZMax, inout half prevZMaxEstimate, inout bool permute, inout half2 hitPixel,
                   half2 invSize, inout half layerThickness)
{
    bool stop = intersecting;
    stepCount++;
   // maxSteps
    for (; (P.x * stepDirection) <= end && stepCount < maxSteps && !stop; P += dP, Q.z += dQ.z, k += dk, stepCount += 1)
    {
        rayZMin = prevZMaxEstimate;
        rayZMax = (dQ.z * 0.5 + Q.z) / (dk * 0.5 + k);
        prevZMaxEstimate = rayZMax;

        if (rayZMin > rayZMax)
        {
            swap(rayZMin, rayZMax);
        }

        hitPixel = permute ? P.yx : P;

        
        sceneZ = SamplScreenRT(RTindex, half4(ScreenPosToUv10(hitPixel * invSize), 0, 0)).r;
        bool isBehind = (rayZMin <= sceneZ);

        if (traceBehind_Old == 1)
        {
            intersecting = isBehind && (rayZMax >= sceneZ - layerThickness);
        }
        else
        {
            intersecting = (rayZMax >= sceneZ - layerThickness);
        }

        stop = traceBehind ? intersecting : isBehind;
    }
    P -= dP, Q.z -= dQ.z, k -= dk;
    return stop;

}



//design for UE,And use lerp to repalce camera space postion iter-----------------
/*
为UE设计，二次步进迭代，使用lerp替代原来的步进方案。
*/
bool rayIterations_LerpForUE(
inout float2 hitUV,
float useX,
float maxLen,

float2 startFrag,
float2 endFrag,
float2 staruv,
float4 startView,
float4 endView,

float steps,
float2 increment, //frist Steps len
float ScaleDelta,
float2 deltaDir,
float thickness,
float2 texSize,
  float originalStepCount
)
{
    float4 positionTo = startView;
    float2 frag = startFrag;
    
    //hit flag ----------------------------------------------------------
    float search0 = 0;
    float search1 = 0;

    int hit0 = 0;
    int hit1 = 0;
    
     //hit param ----------------------------------------------------------
    float viewDistance = startView.z;
    float depth = thickness;
    float deltaX = deltaDir.x;
    float deltaY = deltaDir.y;
    float2 uv = staruv;
    
    float i = 0;
    for (i = 0; i < int(ScaleDelta); ++i)
    {

        originalStepCount++;
        frag += increment;
        uv.xy = frag / texSize;
        
        
       [branch]
        if (uv.x > 1 || uv.y > 1)
        {
            hitUV = 0;

            return 0;
        }
    //    positionTo = texture(positionTexture, uv.xy);
        positionTo = SamplScreenRT(1, half4(uv.xy, 0, 0)).r;
        search1 =
      lerp
        ((frag.y - startFrag.y) / deltaY
        , (frag.x - startFrag.x) / deltaX
        , useX
        );

        search1 = clamp(search1, 0.0, 1.0);

        //get now pixel View space Depth---------------------------------------------
        viewDistance = (startView.z * endView.z) / lerp(endView.z, startView.z, search1);
        
        [branch]
        if (viewDistance < 0.00001)
        {
            hitUV = 0;

            return 0;
        }
        
        depth = viewDistance - positionTo.z;

        
        
        
        //start epual depth -------------------------------------------------------
        if (depth > 0 && depth < thickness)
        {
            hit0 = 1;
            break;
        }
        else
        {
            search0 = search1;
        }
    }

    search1 = search0 + ((search1 - search0) / 2.0);

    steps *= hit0;

    for (i = 0; i < steps; ++i)
    {
        frag = lerp(startFrag.xy, endFrag.xy, search1);
        uv.xy = frag / texSize;
        positionTo = SamplScreenRT(1, half4(uv.xy, 0, 0)).r;

        viewDistance = (startView.z * endView.z) / lerp(endView.z, startView.z, search1);
        depth = viewDistance - positionTo.r;

        if (depth > 0 && depth < thickness)
        {
            hit1 = 1;
            search1 = search0 + ((search1 - search0) / 2);
        }
        else
        {
            float temp = search1;
            search1 = search1 + ((search1 - search0) / 2);
            search0 = temp;
        }
    }

    hitUV = uv;
    return hit1;

}




#define UE_TRACE
bool Linear2D_Trace(int RTindex,
                             half3 csOrigin,
                             half3 csDirection,
                             half2 csZBufferSize,
                             half jitter,
                             int maxSteps,
                             half layerThickness,
                             half traceDistance,
                             in out half2 hitPixel,
                             float stepSize,
                             bool traceBehind,
                             in out half3 csHitPoint,
                             in out half stepCount,
                             float firstSampleStepScale //should be 0-1,
                       
)
{

    half2 invSize = half2(1 / csZBufferSize.x, 1 / csZBufferSize.y);
    hitPixel = half2(-1, -1);

    half nearPlaneZ = -0.01;
    half rayLength = ((csOrigin.z + csDirection.z * traceDistance) < nearPlaneZ) ? ((nearPlaneZ - csOrigin.z) / csDirection.z) : traceDistance;
    half3 csEndPoint = csDirection * rayLength + csOrigin;
    
#ifdef UE_TRACE
    half4 H0 = mul(half4(csOrigin, 1), View.ViewToClip);
    half4 H1 = mul(half4(csEndPoint, 1), View.ViewToClip);
#else
    half4 H0 = mul(projectMatrix, half4(csOrigin, 1));
    half4 H1 = mul(projectMatrix, half4(csEndPoint, 1));
#endif
 
    half k0 = 1 / H0.w;
    half k1 = 1 / H1.w;
    half2 P0 = H0.xy * k0;
    half2 P1 = H1.xy * k1;
    half3 Q0 = csOrigin; // * k0;
    half3 Q1 = csEndPoint; // * k1;


    
    P0 = ScreenPosToUv10(P0);
    P1 = ScreenPosToUv10(P1);
    P0 = clamp(P0, 0, 1);
    P1 = clamp(P1, 0, 1);
    float2 uv = P0;
    
    k0 /= invSize;
    k1 /= invSize;
    P0 /= invSize;
    P1 /= invSize;

    float2 wq = P0;
    float2 wq1 = P1;
    
#if 0
    half yMax = csZBufferSize.y - 0.5;
    half yMin = 0.5;
    half xMax = csZBufferSize.x - 0.5;
    half xMin = 0.5;
    half alpha = 0;

    if (P1.y > yMax || P1.y < yMin)
    {
        half yClip = (P1.y > yMax) ? yMax : yMin;
        half yAlpha = (P1.y - yClip) / (P1.y - P0.y);
        alpha = yAlpha;
    }
    if (P1.x > xMax || P1.x < xMin)
    {
        half xClip = (P1.x > xMax) ? xMax : xMin;
        half xAlpha = (P1.x - xClip) / (P1.x - P0.x);
        alpha = max(alpha, xAlpha);
    }

    P1 = lerp(P1, P0, alpha);
    k1 = lerp(k1, k0, alpha);
    Q1 = lerp(Q1, Q0, alpha);
#endif

 //make the point bigger than 0--------------------------------------
    P1 = (distanceSquared(P0, P1) < 0.0001) ? P0 + half2(0.01, 0.01) : P1;
 
 //Check Dirctor steps size--------------------------------------------
    float2 delta = P1 - P0;
    float2 deltaXY = delta;
    bool permute = abs(delta.x) >= abs(delta.y) ? 1.0 : 0.0;
    float deltaScaled = lerp(abs(delta.y), abs(delta.x), permute) * clamp(firstSampleStepScale, 0.0, 1.0);
    delta = delta / max(deltaScaled, 0.001);
    

    float2 dP = delta; //increment of Pixel 
//Use lerp to get the crrurnt position depth-------------------------------------

   
    
    half2 P = P0;
    //int originalStepCount = 0;

    bool traceBehind_Old = true;
  //  rayIterations(RTindex, traceBehind_Old, traceBehind, P, stepDirection, end, originalStepCount, maxSteps, intersecting, sceneZ, dP, Q, dQ, k, dk, rayZMin, rayZMax, prevZMaxEstimate, permute, hitPixel, invSize, layerThickness);
    
    float2 hitUV = 0;
    float originalStepCount = 0;
    bool hit = rayIterations_LerpForUE(hitUV, (float) permute, rayLength, P0, P1, uv,
    float4(csOrigin, 1), float4(csEndPoint, 1), maxSteps, dP, deltaScaled, deltaXY, layerThickness, csZBufferSize, originalStepCount);
    
    stepCount = originalStepCount;
  //  Q.xy += dQ.xy * stepCount;
  //  csHitPoint = Q * (1 / k);

   // hitPixel = P0 + dP;
    hitPixel = hitUV; // / csZBufferSize; //hitUV;
    return hit;
}

void Linear_2DTrace_SingleSPP(out half4 RayHit_PDF,
float3 Ray_Origin_VS, float3 ViewNormal, float3 Ray_Dir_VS, float4 H,
float2 _SSR_ScreenSize, float Jitter, float _SSR_NumSteps_Linear, float _SSR_CullBack,
float _SSR_Thickness, float _SSR_TraceDistance, int _SSR_TraceBehind,
float _SSR_BackwardsRay, inout float Mask, float _SSR_RayStepSize,
half _SSR_ScreenFade, float firstSampleStepScale,float roughness
)
{
    if (roughness > 0.875)
    {
        RayHit_PDF = 0;
        Mask = 0;
        return;
    }

    
    	//-----Consten Property-------------------------------------------------------------------------
    half Ray_HitMask = 0.0, Ray_NumMarch = 0.0;
    half2 Ray_HitUV = 0.0;
    half3 Ray_HitPoint = 0.0;

    //-----BackwardRay-----------------------------------------------------------------------------
   
    [branch]
    if (_SSR_BackwardsRay == 0 && Ray_Dir_VS.z < 0)
    {
        RayHit_PDF = 0;
        Mask = 0;
        return;
    }

    //-----Ray Trace-----------------------------------------------------------------------------
    half Ray_Bump = max(0.01 * Ray_Origin_VS.z, 0.001);
    
    bool Hit = Linear2D_Trace(1, Ray_Origin_VS + ViewNormal * Ray_Bump, Ray_Dir_VS, _SSR_ScreenSize, Jitter, _SSR_NumSteps_Linear, _SSR_Thickness, _SSR_TraceDistance, Ray_HitUV, _SSR_RayStepSize, _SSR_TraceBehind == 1, Ray_HitPoint, Ray_NumMarch, firstSampleStepScale);
    //Ray_HitUV.xy /= _SSR_ScreenSize;

   
    [branch]
    if (Hit)
    {
        Ray_HitMask = Square(1 - max(2 * half(Ray_NumMarch) / half(_SSR_NumSteps_Linear) - 1, 0));
        Ray_HitMask *= saturate(((_SSR_TraceDistance - dot(Ray_HitPoint - Ray_Origin_VS, Ray_Dir_VS))));

        if (_SSR_CullBack < 1)
        {
            // half3 Ray_HitNormal_WS = tex2Dlod(_CameraGBufferTexture2, half4(Ray_HitUV, 0, 0)).rgb * 2 - 1;
            half3 Ray_HitNormal_WS = SamplScreenRT(8, half4(Ray_HitUV, 0, 0)).rgb * 2 - 1;
            
           // half3 Ray_Dir_WS = mul(_SSR_CameraToWorldMatrix, half4(Ray_Dir_VS, 0)).xyz;
            half3 Ray_Dir_WS = mul(half4(Ray_Dir_VS, 0), (ResolvedView.CameraViewToTranslatedWorld)).rgb;
            
            if (dot(Ray_HitNormal_WS, Ray_Dir_WS) > 0)
                Ray_HitMask = 0;
        }
    }
    
    

    Mask = Square(Ray_HitMask * GetScreenFadeBord(Ray_HitUV, _SSR_ScreenFade));
 //   RayHit_PDF = 1 / _SSR_RayStepSize;//Hit;
    RayHit_PDF = half4(Ray_HitUV, Mask, H.w);
  
   // RayHit_PDF = float4(Ray_Origin_VS + ViewNormal * _SSR_RayStepSize, 0);
   // Mask = Square(Ray_HitMask * GetScreenFadeBord(Ray_HitUV, _SSR_ScreenFade));
}



void TRACELIB_MarkFunc()
{
    return;

#endif 