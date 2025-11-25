//
// Glass shader
//
varying vec3  Normal;
varying vec3  EyeDir;
varying vec4  EyePos;
varying vec4 color;

void main(void) 
{
  Normal = normalize(gl_NormalMatrix * gl_Normal);
  vec4 pos = gl_ModelViewMatrix * gl_Vertex;
  EyeDir = pos.xyz;
  EyePos = gl_ModelViewProjectionMatrix * gl_Vertex;
  color = gl_Color;
  gl_Position = ftransform();   
}