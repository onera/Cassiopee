#version 400 compatibility
/*
    Phong (two sides)
*/

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

void main (void)
{ 
  vec3 Nv = v2f_out.view_normal.xyz;
  vec3 P  = v2f_out.mv_position.xyz;

  vec3 N = normalize(Nv);
  vec3 E = normalize(-P);
  vec3 L = normalize(gl_LightSource[0].position.xyz-P);
  float dotNL = dot(N, L);
  if (dotNL < 0.) { N = -N; dotNL = -dotNL; }

  vec3 R = normalize(-reflect(L,N));
  vec4 Iamb = gl_LightSource[0].ambient; 
  vec4 Idiff = gl_LightSource[0].diffuse*max(dotNL, 0.0);
  vec4 Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0),0.25 * gl_FrontMaterial.shininess);
  vec4 col = Iamb + v2f_out.color * Idiff + Ispec;

  gl_FragColor = clamp(col, 0., 1.);
  gl_FragColor.a = v2f_out.color.a;
}
