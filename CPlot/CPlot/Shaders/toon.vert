// Toon shader
varying vec3 normal, lightDir;
varying vec4 color;

void main()
{
  vec4 ecPos;
  ecPos = vec4(gl_ModelViewMatrix * gl_Vertex);
  lightDir = normalize(vec3(gl_LightSource[0].position) - ecPos.xyz);
  normal = normalize(gl_NormalMatrix * gl_Normal);
  
  color = gl_Color;
  gl_Position = ftransform();
}