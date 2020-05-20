// Sphere billboarding
// Remplace un quad par une sphere

varying vec4 initColor;
varying vec3 P;

void main()
{
  P = vec3(gl_ModelViewMatrix * gl_Vertex);
  gl_TexCoord[0] = gl_MultiTexCoord0;
  initColor = gl_Color;
  gl_Position = ftransform();
}
