// Display direct r/g/b + light

varying vec3 Nv;
varying vec3 P;
varying vec4 color;
varying vec4 vertex;

void main()
{ 
  P = vec3(gl_ModelViewMatrix * gl_Vertex);
  Nv = gl_NormalMatrix * gl_Normal;
  vertex = gl_Vertex;
  color = gl_Color;
  gl_Position = ftransform();
}

