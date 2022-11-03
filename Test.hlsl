
if(roughness>0.9f)return 0;
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
 _SSR_ScreenFade,firstSampleStepScale,
roughness
);
    
//param out put
    Mask_out =
Mask;

    RayHit_PDF_out = RayHit_PDF;

    
    return 1;