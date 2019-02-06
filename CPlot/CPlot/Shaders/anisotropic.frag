#version 400 compatibility
// Metal (anisotropic)
/*varying vec3 MCposition;
varying vec4 color;
varying vec3 normal;
varying vec3 tangent;
varying vec3 light;
varying vec3 eye;
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

uniform vec3 intensity;
uniform float bump;
uniform sampler3D Noise;
uniform int shadow;
uniform sampler2D ShadowMap;

void main(void) 
{
  vec3 MCposition = v2f_out.vdata1.xyz;
  vec4 color      = v2f_out.color;
  vec3 normal     = v2f_out.nrm_view_normal.xyz;
  vec3 tangent    = v2f_out.vdata2.xyz;
  vec3 light      = v2f_out.vdata3.xyz;
  vec3 eye        = -v2f_out.mv_position.xyz;
  vec4 vertex     = v2f_out.position;

  vec3 ambient = 0.00001 * vec3(gl_LightSource[0].ambient);
  vec3 diffuse = 0.1 * vec3(gl_LightSource[0].diffuse); // 0.2
  vec3 specular = 0.7 * vec3(gl_LightSource[0].specular); // 0.6
  float exponent = 0.4 * gl_FrontMaterial.shininess; // 0.5

  vec3 n = normalize(normal);
  vec3 t = normalize(tangent);
  
  if (bump > 0.) // bump map, modifie tangent
  {  
     // granite noise
     /*
     vec4  noisevec = texture3D(Noise, MCposition);
     float val = min(1.0, noisevec[3] * 18.0);
     noisevec = texture3D(Noise, MCposition+vec3(0.001,0.,0.));
     float nx = min(1.0, noisevec[3] * 18.0) - val;
     noisevec = texture3D(Noise, MCposition+vec3(0.,0.001,0.));
     float ny = min(1.0, noisevec[3] * 18.0) - val;
     noisevec = texture3D(Noise, MCposition+vec3(0.,0.,0.001));
     float nz = min(1.0, noisevec[3] * 18.0) - val;
     */

     // cloud noise
     vec4  noisevec = texture3D(Noise, MCposition);
     float val = (noisevec[0] + noisevec[1] + noisevec[2]+ noisevec[3])*1.7;
     val = abs(2.0 * val - 1.0);
     noisevec = texture3D(Noise, MCposition+vec3(0.001,0.,0.));
     float val2 = (noisevec[0] + noisevec[1] + noisevec[2] + noisevec[3])*1.7;
     val2 = abs(2.0 * val2 - 1.0);
     float nx = val2 - val;
     noisevec = texture3D(Noise, MCposition+vec3(0.,0.001,0.));
     val2 = (noisevec[0] + noisevec[1] + noisevec[2] + noisevec[3])*1.7;
     val2 = abs(2.0 * val2 - 1.0);
     float ny = val2 - val;
     noisevec = texture3D(Noise, MCposition+vec3(0.,0.,0.001));
     val2 = (noisevec[0] + noisevec[1] + noisevec[2] + noisevec[3])*1.7;
     val2 = abs(2.0 * val2 - 1.0);
     float nz = val2 - val;

     // Theoriquement juste 
     //vec3 bumpN = n + bump * 0.5 * gl_NormalMatrix * vec3(nx, ny ,nz);
     //bumpN = normalize(bumpN);
     //t = t-dot(bumpN,t)*bumpN;
     //t = normalize(t);

     // joli
     vec3 bumpN = bump * 1.5 * gl_NormalMatrix * vec3(nx, ny ,nz);
     t = t-dot(bumpN,t)*bumpN;
     t = normalize(t);
  }

  vec3 l = normalize(light);
  vec3 e = normalize(eye);
  float dl = dot(t, l);
  float de = dot(t, e);
  float dr = sqrt(1.0 - dl*dl);

  //vec3 colorl = vec3(color.r, color.g, color.b)*0.25;
  //vec3 col = ambient + intensity * diffuse * max(0.0, dr) + specular * intensity * pow(max(dr*sqrt(1.0 - de*de)-dl*de,0.), exponent);
  //colorl += col;
  //colorl = clamp(colorl, 0., 1.);

  vec3 colorl = vec3(color.r, color.g, color.b);
  vec3 col = colorl*0.25 + ambient + diffuse * colorl * max(0.0, dr) + specular * intensity * pow(max(dr*sqrt(1.0 - de*de)-dl*de,0.), exponent);
  colorl = clamp(col, 0., 1.);

  vec4 color2 = vec4(colorl, color.a);

  float shadowValue = 1.;
  if (shadow > 0)
  {
  vec3 L = normalize(gl_LightSource[0].position.xyz+eye);
    
  // Coords -> texCoords
  vec4 ShadowCoord = gl_TextureMatrix[0] * vertex;
  vec4 shadowCoordinateW = ShadowCoord / ShadowCoord.w;

  // Used to lower moire pattern and self-shadowing
  //shadowCoordinateW.z -= 0.00001;
  float dotNL = dot(n, L);
  shadowCoordinateW.z -= (abs(dotNL)+0.1)*0.00001;

  // Z buffer du point dans la texture rendu du pt de vue de la lumiere
  float distanceFromLight = texture2D(ShadowMap, shadowCoordinateW.st).r;
  float s = shadowCoordinateW.s;
  float t = shadowCoordinateW.t;      
  if (ShadowCoord.w > 0.0 && s > 0.001 && s < 0.999 && t > 0.001 && t < 0.999)
      shadowValue = distanceFromLight < shadowCoordinateW.z ? 0.5 : 1.0;
  }

   gl_FragColor = shadowValue * color2;
   gl_FragColor.a = color.a;
}
