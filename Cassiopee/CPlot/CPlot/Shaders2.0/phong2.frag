/*
    Phong (two sides)
*/
varying vec3 Nv;
varying vec3 P;
varying vec4 color;

void main (void)
{ 
  vec3 N = normalize(Nv);
  vec3 E = normalize(-P);
  vec3 L = normalize(gl_LightSource[0].position.xyz-P);
  float dotNL = dot(N, L);
  if (dotNL < 0.) { N = -N; dotNL = -dotNL; }

  vec3 R = normalize(-reflect(L,N));
  vec4 Iamb = gl_LightSource[0].ambient; 
  vec4 Idiff = gl_LightSource[0].diffuse*max(dotNL, 0.0);
  vec4 Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0),0.25 * gl_FrontMaterial.shininess);
  vec4 col = Iamb + color * Idiff + Ispec;

  gl_FragColor = clamp(col, 0., 1.);
  gl_FragColor.a = color.a;
}
