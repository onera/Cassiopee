// isos + granite

varying vec3 Nv;
varying vec3 P;
varying vec4 color;
varying vec3 MCposition;
varying vec4 vertex;
uniform float Scale;

void main()
{ 
  P = vec3(gl_ModelViewMatrix * gl_Vertex);
  Nv = gl_NormalMatrix * gl_Normal;
  vertex = gl_Vertex;
  MCposition = vec3(gl_Vertex) * Scale;
  color = gl_Color;
  gl_Position = ftransform();
}

