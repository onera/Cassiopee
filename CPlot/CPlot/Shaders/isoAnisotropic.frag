// Shader iso + anisotropic
//
// iso + metal
varying vec4 color;
varying vec3 normal;
varying vec3 tangent;
varying vec3 light;
varying vec3 eye;
varying vec4 vertex;

uniform vec3 intensity;
uniform sampler1D colormap;

// Number of iso
uniform float niso;
uniform float alpha; // colormap range
uniform float beta;
uniform float blend;
uniform float edgeStyle;
uniform int shadow;
uniform sampler2D ShadowMap;

void main()
{
  float f, fs;
  int vali;
  f = color.r; f = alpha*f + beta;
  fs = f;
  vali = int(f*niso);
  f = float(vali)/niso;
  vec3 val; 
  val = vec3(texture1D(colormap, f));

  float df, borne, borne1, borne2;
  df = fwidth(fs);
  borne = edgeStyle*df;
  if (fs-f < borne) 
  { 
     borne1 = 0.8*borne; borne2 = 0.2*borne;
     if (fs-f > borne1) { df = (fs-f-borne1)/borne2; val = val * df; }
     else if (fs-f < borne2) { df = (-fs+f+borne2)/borne2; val = val * df; }
     else { val = vec3(0.); }
  }
  vec4 color2 = vec4(val.r, val.g, val.b, blend); 

  // anisotropic light
  vec3 ambient = 0.001 * vec3(gl_LightSource[0].ambient);
  vec3 diffuse = 0.003 * vec3(gl_LightSource[0].diffuse);
  vec3 specular = 0.7 * vec3(gl_LightSource[0].specular);
  float exponent = 0.5 * gl_FrontMaterial.shininess;
 
  vec3 n = normalize(normal);
  vec3 t = normalize(tangent);
  vec3 l = normalize(light);
  vec3 e = normalize(eye);
  float dl = dot(t, l);
  float de = dot(t, e);
  float dr = sqrt(1.0 - dl*dl);
  vec3 colorl = vec3(color2.r, color2.g, color2.b)*0.4;
  vec3 col = ambient + diffuse * intensity * max(0.0, dr) + specular * intensity * pow(max(dr*sqrt(1.0 - de*de)-dl*de,0.) , exponent);

  colorl += col;
  colorl = clamp(colorl, 0., 1.);  
    
  color2 = vec4(colorl, color2.a);
    
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
  gl_FragColor.a = blend;
}
