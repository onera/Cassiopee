#version 400 compatibility
//
// Envmap shader
//
const vec3 Xunitvec = vec3(1.0, 0.0, 0.0);
const vec3 Yunitvec = vec3(0.0, 1.0, 0.0);

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
//varying vec4 color;
//varying vec4 vertex;
//varying vec3 Normal;
//varying vec3 EyeDir;

uniform float MixRatio;
uniform float intensity;
uniform sampler2D EnvMap;
uniform int shadow;
uniform sampler2D ShadowMap;

void main (void)
{
    vec3 EyeDir = v2f_out.mv_position.xyz;
    vec3 BaseColor = vec3(v2f_out.color.r, v2f_out.color.g, v2f_out.color.b);
    vec3 LightPos  = gl_LightSource[0].position.xyz;
    vec3 L = normalize(LightPos - EyeDir);
    float LightIntensity = dot(v2f_out.nrm_view_normal.xyz, L);
    vec3 N = v2f_out.nrm_view_normal.xyz;// Normal;
    if (LightIntensity < 0.) 
    { N = -N; LightIntensity = - LightIntensity; }

    // Compute reflection vector
    vec3 reflectDir = reflect(EyeDir, N);

    // Compute altitude and azimuth angles
    vec2 index;

    index.y = dot(normalize(reflectDir), Yunitvec);
    reflectDir.y = 0.0;
    index.x = dot(normalize(reflectDir), Xunitvec) * 0.5;

    // Translate index values into proper range
    if (reflectDir.z >= 0.0)
        index = (index + 1.0) * 0.5;
    else
    {
        index.t = (index.t + 1.0) * 0.5;
        index.s = (-index.s) * 0.5 + 1.0;
    }
    
    // if reflectDir.z >= 0.0, s will go from 0.25 to 0.75
    // if reflectDir.z <  0.0, s will go from 0.75 to 1.25, and
    // that's OK, because we've set the texture to wrap.
  
    // Do a lookup into the environment map.
    vec3 envColor = vec3(texture2D(EnvMap, index));

    // Add lighting to base color and mix
    vec3 base = LightIntensity * intensity * BaseColor;
    envColor = mix(envColor, base, MixRatio);
    
    vec4 col = vec4(envColor, v2f_out.color.a);

     float shadowValue = 1.;
     if (shadow > 0)
     {
     // Coords -> texCoords
     vec4 ShadowCoord = gl_TextureMatrix[0] * v2f_out.position;//vertex;
     vec4 shadowCoordinateW = ShadowCoord / ShadowCoord.w;

     // Used to lower moire pattern and self-shadowing
     //shadowCoordinateW.z -= 0.00001;
     float dotNL = dot(N, L);
     shadowCoordinateW.z -= (abs(dotNL)+0.1)*0.00001;

     // Z buffer du point dans la texture rendu du pt de vue de la lumiere
     float distanceFromLight = texture2D(ShadowMap, shadowCoordinateW.st).r;
     float s = shadowCoordinateW.s;
     float t = shadowCoordinateW.t;      
     if (ShadowCoord.w > 0.0 && s > 0.001 && s < 0.999 && t > 0.001 && t < 0.999)
       shadowValue = distanceFromLight < shadowCoordinateW.z ? 0.5 : 1.0;
     }

    gl_FragColor = shadowValue * col;
    gl_FragColor.a = v2f_out.color.a;
}
