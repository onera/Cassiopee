#version 400 compatibility
//
// Fragment shader for producing clouds (mostly sunny)
//
/*varying float LightIntensity; 
varying vec3 MCposition;
varying vec4 initColor;
varying vec4 vertex;*/
in V2F_OUT
{
    vec4 position;
    vec4 mv_position;
    vec4 mvp_position;
    vec4 view_normal;
    vec4 nrm_view_normal;
    vec4 color;
    vec4 vdata1, vdata2, vdata3, vdata4;
} v2f_out;

uniform sampler3D Noise;
uniform vec3 Offset;
uniform vec3 CloudColor;
uniform int shadow;
uniform sampler2D ShadowMap;

void main (void)
{
    float LightIntensity = v2f_out.vdata1.w;
    vec3  MCposition     = v2f_out.vdata1.xyz;
    vec4 initColor       = v2f_out.color;
    vec4 vertex          = v2f_out.position;
    
    vec3 SkyColor = vec3(initColor.r, initColor.g, initColor.b); 
    vec4 noisevec = texture3D(Noise, 1.2 * (vec3(0.5) + MCposition + Offset));

    float intensity = (noisevec[0] + noisevec[1] +
                       noisevec[2] + noisevec[3]) * 1.7;

    intensity = 1.95 * abs(2.0 * intensity - 1.0);
    intensity = clamp(intensity, 0.0, 1.0);
    vec3 color = mix(CloudColor, SkyColor, intensity) * LightIntensity;
    color = clamp(color, 0.0, 1.0);

    float shadowValue = 1.;
    if (shadow > 0)
    {
    // Coords -> texCoords
    vec4 ShadowCoord = gl_TextureMatrix[0] * vertex;
    vec4 shadowCoordinateW = ShadowCoord / ShadowCoord.w;

    // Used to lower moire pattern and self-shadowing
    //shadowCoordinateW.z -= 0.00001;
    shadowCoordinateW.z -= (abs(LightIntensity*2./3.)+0.1)*0.00001;

    // Z buffer du point dans la texture rendu du pt de vue de la lumiere
    float distanceFromLight = texture2D(ShadowMap, shadowCoordinateW.st).r;
    float s = shadowCoordinateW.s;
    float t = shadowCoordinateW.t;      
    if (ShadowCoord.w > 0.0 && s > 0.001 && s < 0.999 && t > 0.001 && t < 0.999)
      shadowValue = distanceFromLight < shadowCoordinateW.z ? 0.5 : 1.0;
    }


    color = shadowValue * color;
    gl_FragColor = vec4(color, initColor.a);
    gl_FragColor.a = initColor.a;
}
