// Shader iso colored lines

varying vec4 color;
varying vec3 Nv;
varying vec3 P;
varying vec4 vertex;

// Number of iso
uniform float niso;
// light is on?
uniform int lightOn;
// edge style
uniform float edgeStyle;
uniform sampler1D colormap;
uniform float alpha; // colormap range
uniform float beta;
uniform float blend;
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

  vec3 bgColor;
  if (lightOn == 0) bgColor = vec3(1.,1.,1.);
  else bgColor = vec3(0.6, 0.6, 0.6);

  float df, borne, borne1, borne2;
  //df = color.g;
  df = fwidth(fs);
  borne = (edgeStyle+1.)*df;
  vec4 color2;
  vec3 N = normalize(Nv);
  vec3 L = normalize(gl_LightSource[0].position.xyz-P);

  if (fs-f > borne) 
  {
    color2 = vec4(bgColor.r, bgColor.g, bgColor.b, blend);
    // light on!
    if (lightOn == 1)
    {
      vec3 E, R;
      vec4 Iamb, Idiff, Ispec;      
      E = normalize(-P);
      if (dot(N, L) < 0.) N = -N;

      R = normalize(-reflect(L,N));
      Iamb = gl_LightSource[0].ambient; 
      Idiff = gl_LightSource[0].diffuse*max(dot(N,L), 0.0);
      Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0),0.2 * gl_FrontMaterial.shininess);

      vec4 col = Iamb + color2 * Idiff + Ispec;
      color2 = clamp(col, 0., 1.);
      color2.a = blend;
    }
  } 
  else 
  {
     borne1 = 0.8*borne; borne2 = 0.2*borne;
     if (fs-f > borne1) 
     { df = (fs-f-borne1)/borne2; val = mix(val, bgColor, df); }
     else if (fs-f < borne2) 
     { df = (-fs+f+borne2)/borne2; val = mix(val, bgColor, df); }
     color2 = vec4(val.r, val.g, val.b, 1.);
  }

  float shadowValue = 1.;
  if (shadow > 0)
  {
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

  df = color2.a;
  gl_FragColor = shadowValue * color2;
  gl_FragColor.a = df;
}

