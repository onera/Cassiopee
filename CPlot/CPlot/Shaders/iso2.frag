// Shader continuous-colormap

varying vec4 color;
varying vec3 Nv;
varying vec3 P;

// light is on?
uniform int lightOn;
//colormap
uniform sampler1D colormap;
uniform float alpha;
uniform float beta;
uniform float blend;

void main()
{
  float r, g, b, f;
  f = color.r;
  vec3 val;
  f = alpha*f+beta;
  f = clamp(f, 0.0f, 1.0f);
  val = vec3(texture1D(colormap, f));
  
  vec4 color2 = vec4(val.r, val.g, val.b, blend);

  // light on!
  if (lightOn == 1)
  {
    vec3 N, E, L, R;
    vec4 Iamb, Idiff, Ispec;

    N = normalize(Nv);
    E = normalize(-P);
    L = normalize(gl_LightSource[0].position.xyz-P);
    if (dot(N, L) < 0.) N = -N;

    R = normalize(-reflect(L,N));
    Iamb = gl_LightSource[0].ambient; 
    Idiff = gl_LightSource[0].diffuse*max(dot(N,L), 0.0);
    Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0),0.2 * gl_FrontMaterial.shininess);

    vec4 col = Iamb + color2 * Idiff + Ispec;
    color2 = clamp(col, 0., 1.);
    color2.a = blend;
  }
  
  gl_FragColor = color2;
}
