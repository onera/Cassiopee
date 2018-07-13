// Shader iso + gooch
// = color isolines + blended background
varying vec4 color;

// Number of iso
uniform float niso;
// edge style
uniform float edgeStyle;
uniform sampler1D colormap;
uniform float alpha; // colormap range
uniform float beta;
uniform float blend;

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

  vec3 bgColor = vec3(1.,1.,1.);

  float df, borne, borne1, borne2;
  df = fwidth(fs);
  borne = (edgeStyle+1.)*df;
  vec4 color2; 
  if (fs-f < borne)
  { 
     borne1 = 0.8*borne; borne2 = 0.2*borne;
     if (fs-f > borne1) { df = (fs-f-borne1)/borne2; val = mix(val, bgColor, df); }
     else if (fs-f < borne2) { df = (-fs+f+borne2)/borne2; val = mix(val, bgColor, df); }
    color2 = vec4(val.r, val.g, val.b,  blend);
  } 
  else color2 = vec4(bgColor.r, bgColor.g, bgColor.b, 0.);
 
  gl_FragColor = color2;
}

