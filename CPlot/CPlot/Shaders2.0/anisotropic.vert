//
// Anisotropic light (brushed metal)
//
varying vec3 MCposition;
varying vec4 color;
varying vec3 normal;
varying vec3 tangent;
varying vec3 light;
varying vec3 eye;
varying vec4 vertex;
uniform float Scale;

void main(void) 
{
    MCposition = vec3(gl_Vertex) * Scale;
    vec3 t = vec3(1.,0.1,.1);
    t = normalize(t);
    normal = normalize(gl_NormalMatrix * gl_Normal);

    //vec3 t = vec3(normal.z,normal.z,-normal.x-normal.y);
    //t = normalize(t);

    tangent = t - dot(t,normal)*normal; // strand dir
    //if (dot(tangent, tangent) < 1.e-12)
    //{
    //  t = vec3(0,1.,0.);
    //  tangent = t - dot(t,normal)*normal;
    //}   
    tangent = normalize(tangent);
    //bitangent = normalize(cross(tangent, normal));
    
    vec3 P = vec3(gl_ModelViewMatrix * gl_Vertex);
    vec3 lightPosition = gl_LightSource[0].position.xyz;
    light = lightPosition - P;
    eye = -P;
    vertex = gl_Vertex;
    color = gl_Color;
    gl_Position = ftransform();
}
