// isos + anisotropic

varying vec4 color;
varying vec3 normal;
varying vec3 tangent;
varying vec3 light;
varying vec3 eye;
varying vec4 vertex;

void main()
{ 
  vec3 t = vec3(1.,0.,0.);
  normal = normalize(gl_NormalMatrix * gl_Normal);

  //vec3 t = vec3(normal.z,normal.z,-normal.x-normal.y);
  //t = normalize(t);

  tangent = normalize(t - dot(t,normal)*normal);
    
  vec3 P = vec3(gl_ModelViewMatrix * gl_Vertex);
  vec3 lightPosition = gl_LightSource[0].position.xyz;
  light = lightPosition - P;
  eye = -P;
  vertex = gl_Vertex;
  color = gl_Color;
  gl_Position = ftransform();
}

