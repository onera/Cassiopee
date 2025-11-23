// Gooch (hand drawing) shader
varying float NdotL;
varying vec3 ReflectVec;
varying vec3 ViewVec;
varying vec3 ecPos;
varying vec3 tnorm;
varying vec4 color;
varying vec4 vertex;

void main(void)
{
  ecPos = vec3(gl_ModelViewMatrix * gl_Vertex);
  tnorm = normalize(gl_NormalMatrix * gl_Normal);
  vec3 lightVec = normalize(gl_LightSource[0].position.xyz - ecPos);
  ReflectVec = normalize(reflect(-lightVec, tnorm));
  ViewVec = normalize(-ecPos);
  NdotL = (dot(lightVec, tnorm) + 1.0) * 0.5;
  vertex = gl_Vertex;
  color = gl_Color;
  gl_Position = ftransform();
}