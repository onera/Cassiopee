// Flat iso shader (no light)
varying vec3 Nv;
varying vec3 P;
varying vec4 color;
varying vec4 vertex;
// Number of iso
uniform float niso;
// edge style
uniform float edgeStyle;
uniform sampler1D colormap;
uniform float alpha; // colormap range
uniform float beta;
uniform float blend;
uniform int shadow;
uniform sampler2D ShadowMap;

void main (void)
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

  gl_FragColor = shadowValue * color2;
  gl_FragColor.a = blend;
}
