#version 400 compatibility
// Gooch (hand drawing) shader

/*varying float NdotL;
varying vec3 ReflectVec;
varying vec3 ViewVec;
varying vec4 color;
varying vec3 ecPos;
varying vec3 tnorm;
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

uniform float specularFactor;
uniform float exponent;
uniform int shadow;
uniform sampler2D ShadowMap;

void main (void)
{
    float NdotL     = v2f_out.vdata2.w;
    vec3 ReflectVec = v2f_out.vdata1.xyz;
    vec3 ViewVec    = v2f_out.vdata2.xyz;
    vec4 color      = v2f_out.color;
    vec3 ecPos      = v2f_out.mv_position.xyz;
    vec3 tnorm      = v2f_out.nrm_view_normal.xyz;
    vec4 vertex     = v2f_out.position;

    vec4 colorf;
    vec3 SurfaceColor = vec3(color.r, color.g, color.b);
    vec3 WarmColor = vec3(color.r-0.1, color.g-0.1, color.b-0.1);
    vec3 CoolColor = vec3(color.r+0.2, color.g+0.2, color.b+0.2);
    float DiffuseWarm = 0.45;
    float DiffuseCool = 0.45;

    vec3 kcool   = min(CoolColor + DiffuseCool*SurfaceColor, 1.0);
    vec3 kwarm   = min(WarmColor + DiffuseWarm*SurfaceColor, 1.0); 
    vec3 kfinal  = mix(kcool, kwarm, NdotL);

    vec3 nreflect = normalize(ReflectVec);
    vec3 nview    = normalize(ViewVec);

    float spec    = max(dot(nreflect, nview), 0.0);
    spec          = pow(spec, 32.0*(2.-specularFactor));
       
    colorf = vec4(min(kfinal+spec, 1.0), 1.0);

    // silhouette
    float angle = dot(tnorm, normalize(ecPos));
    angle = abs(angle);
    angle = pow(angle, -1.1);
    if (angle > 3.-exponent)  {colorf = vec4(0.,0.,0.,1.);}

    float shadowValue = 1.;
    if (shadow > 0)
    {
    // Coords -> texCoords
    vec4 ShadowCoord = gl_TextureMatrix[0] * vertex;
    vec4 shadowCoordinateW = ShadowCoord / ShadowCoord.w;

    // Used to lower moire pattern and self-shadowing
    //shadowCoordinateW.z -= 0.00001;
    shadowCoordinateW.z -= (abs(NdotL)+0.1)*0.00001;

    // Z buffer du point dans la texture rendu du pt de vue de la lumiere
    float distanceFromLight = texture2D(ShadowMap, shadowCoordinateW.st).r;
    float s = shadowCoordinateW.s;
    float t = shadowCoordinateW.t;      
    if (ShadowCoord.w > 0.0 && s > 0.001 && s < 0.999 && t > 0.001 && t < 0.999)
      shadowValue = distanceFromLight < shadowCoordinateW.z ? 0.5 : 1.0;
    }

    gl_FragColor = shadowValue * colorf;
    gl_FragColor.a = color.a;
}
