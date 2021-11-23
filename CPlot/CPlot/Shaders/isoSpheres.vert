// iso + Sphere billboarding

varying vec3 P;
varying vec4 color;

void main()
{
  P = vec3(gl_ModelViewMatrix * gl_Vertex);
  gl_TexCoord[0] = gl_MultiTexCoord0;
  color = gl_Color;
  gl_Position = ftransform();
}
