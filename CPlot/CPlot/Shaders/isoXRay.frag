/* iso + XRay */ 
varying vec3 N;
varying vec3 I;
varying vec4 color;
// Number of iso
uniform float niso;
// edge style
uniform float edgeStyle;
uniform sampler1D colormap;
uniform float alpha; // colormap range
uniform float beta;
uniform float blend;
uniform float EdgeFalloff;

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

  float opacity = dot(normalize(N), normalize(I));
  opacity = abs(opacity);
  opacity = 1.0 - pow(opacity, EdgeFalloff);
  opacity = clamp(opacity, 0., 1.);
    
  gl_FragColor = color2;
  gl_FragColor.a = opacity*blend;
}
