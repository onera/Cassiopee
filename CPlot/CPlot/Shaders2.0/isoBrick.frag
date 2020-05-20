// Shader iso + brick
// = color bands + bump mapped edges
varying vec4 color;
varying vec3 Nv;
varying vec3 P;

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
  borne = (edgeStyle+5.)*df;
  if (fs-f < borne)
  { 
     borne1 = 0.8*borne; borne2 = 0.2*borne;
     if (fs-f > borne1) { df = (fs-f-borne1)/borne2; val = val * df; }
     else if (fs-f < borne2) { df = (-fs+f+borne2)/borne2; val = val * df; }
     else { val = vec3(1.); }
  } 
  
  vec4 color2 = vec4(val.r, val.g, val.b, blend);
  vec3 N = normalize(Nv);
  vec3 L = normalize(gl_LightSource[0].position.xyz-P);

  // light on!
  if (lightOn == 1)
  {
    vec3 E, R;
    vec4 Iamb, Idiff, Ispec;    
    E = normalize(-P);
    L = normalize(gl_LightSource[0].position.xyz-P);
    if (dot(N, L) < 0.) N = -N;
    R = normalize(-reflect(L,N));
    Iamb = gl_LightSource[0].ambient; 
    Idiff = gl_LightSource[0].diffuse*max(dot(N,L), 0.0);
    Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0),0.2 * gl_FrontMaterial.shininess);

    if (fs-f < borne) { Ispec=vec4(0.); Iamb = vec4(0.); }

    vec4 col = Iamb + color2 * Idiff + Ispec;
    color2 = clamp(col, 0., 1.);
    color2.a = blend;
  }

  gl_FragColor = color2;
}

