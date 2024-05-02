// isos + brick

varying vec3 Nv;
varying vec3 P;
varying vec4 color;

void main()
{ 
  P = vec3(gl_ModelViewMatrix * gl_Vertex);
  Nv = gl_NormalMatrix * gl_Normal;
  color = gl_Color;
  gl_Position = ftransform();
}
