#version 400 compatibility
// Flat shader (no light)

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

uniform int shadow;
uniform sampler2D ShadowMap;

void main (void)
{ 
  vec3 Nv     = v2f_out.view_normal.xyz;
  vec3 P      = v2f_out.mv_position.xyz;
  vec4 vertex = v2f_out.position;

  float shadowValue = 1.;
  if (shadow > 0)
  {
  vec3 N = normalize(Nv);
  vec3 L = normalize(gl_LightSource[0].position.xyz-P);
    
  // Coords -> texCoords
  vec4 ShadowCoord = gl_TextureMatrix[0] * vertex;
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
  gl_FragColor = shadowValue * v2f_out.color;
  gl_FragColor.a = v2f_out.color.a;
}
