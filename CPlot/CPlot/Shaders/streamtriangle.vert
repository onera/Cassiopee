#version 150 compatibility

out Vertex {
    vec4 P0;// Coordonnées dans l'espace réel
    vec4 vP;// Coordonnées dans l'espace caméra
    vec4 e1;
    vec4 e3;
    vec4 position;
    vec4 normal;
    vec4 color;
    vec4 translation;
} vertex;
uniform float scale;
uniform int fix_length;

void main()
{
    gl_Position = ftransform();
    vertex.position = gl_Vertex;
    vertex.P0    = gl_Position;
    vertex.vP    = gl_ModelViewMatrix * gl_Vertex;
    vertex.color = gl_Color;
    vertex.translation = gl_ModelViewProjectionMatrix*(5.E-2*scale*vec4(gl_Normal,0.));
    if ( vertex.translation.z > 0. ) 
	    vertex.translation = vec4(vertex.translation.xy, -vertex.translation.z, vertex.translation.w);
    vec4 transpos = vec4(gl_Color.xyz,0.) - vec4(0.5,0.5,0.5,0.0);
    transpos = (1-fix_length)*transpos + fix_length*normalize(transpos);
    float s = 0.5*scale;
    bool bx = (transpos.x > 1.E-6)||(transpos.x < -1.E-6);
    bool by = (transpos.y > 1.E-6)||(transpos.y < -1.E-6);
    bool bz = (transpos.z > 1.E-6)||(transpos.z < -1.E-6);
    vertex.e1 = s*vec4((abs(transpos.x)+abs(transpos.y)+abs(transpos.z)),0.,0.,0.);
    /*vertex.e1 = by ? gl_ModelViewProjectionMatrix*(s*vec4(transpos.y,-transpos.x,0.,0.)) :
	             gl_ModelViewProjectionMatrix*(s*vec4(transpos.z,0.,-transpos.x,0.));*/
    vertex.normal = vec4(gl_NormalMatrix * gl_Normal,0.);
    vertex.e3 = gl_ModelViewProjectionMatrix*(scale*vec4(transpos.xyz,0.));
}
