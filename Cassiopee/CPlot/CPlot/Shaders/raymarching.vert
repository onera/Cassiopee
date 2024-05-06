//
// Vertex shader for ray marching through a 3D texture
//
varying vec4 color;
varying vec3 cam;
uniform vec3 dx;
uniform vec3 xo;
uniform vec3 poscam;

void main() 
{
    gl_TexCoord[0] = gl_MultiTexCoord0;
    color = gl_Color;
    //cam = vec3(gl_ModelViewMatrix * gl_Vertex);
    //vec4 eye = vec4(2.,0.5,0.5,0);
    vec4 caml = vec4(poscam, 0); // modelview
    //vec4 cam2 = transpose(gl_ModelViewMatrix)*caml;     // world
  
    cam = vec3(caml);
    //cam = dx*cam + xo;
    cam.x = dx.x * (cam.x - xo.x);
    cam.y = dx.y * (cam.y - xo.y);
    cam.z = dx.z * (cam.z - xo.z); // texture coord.

    gl_Position = ftransform();
}