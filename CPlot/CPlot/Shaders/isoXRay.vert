/* isos + XRay */
	 
varying vec3 N;
varying vec3 I;
varying vec4 color;

void main()
{
    // Position in clip space
    I = vec3(gl_ModelViewMatrix * gl_Vertex);
 
    // Normal transform (transposed model-view inverse)
    N = gl_NormalMatrix * gl_Normal;
	
    color = gl_Color;
    gl_Position = ftransform();
 
}	 