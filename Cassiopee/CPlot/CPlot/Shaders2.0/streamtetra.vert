#version 150 compatibility

out Vertex {
    vec4 P0; // Coordonnée du point dans l'espace image
    vec4 vP; // Coordonnée du point dans l'espace caméra
    vec4 e1; // Direction du premier triangle
    vec4 e2; // Direction du second triangle orthogonal au premier
    vec4 e3; // Direction suivant le vecteur "vitesse"
    vec4 position;// Position du point dans l'espace réel
    vec4 normal; // Normal à la surface
    vec4 color;  // Une des trois couleurs à donner au vecteur vitesse
    vec4 translation; // Petite translation pour que le vecteur ne soit pas
                      // confondu avec le triangle portant le champ
} vertex;
uniform float scale;// Taille du vecteur
uniform int fix_length;// Normalisation (= 1) du champs ou non ( =0 )

void main()
{
    gl_Position = ftransform();
    vertex.position = gl_Vertex;
    vertex.P0 = gl_Position;
    vertex.vP    = gl_ModelViewMatrix * gl_Vertex;
    vertex.color    = gl_Color;
    vertex.translation = gl_ModelViewProjectionMatrix*(5.E-3*scale*vec4(gl_Normal,0.));
    if ( vertex.translation.z > 0. ) 
        vertex.translation = vec4(vertex.translation.xy, -vertex.translation.z, vertex.translation.w);
    vec4 transpos = vec4(gl_Color.xyz,0.) - vec4(0.5,0.5,0.5,0.0);
    transpos = (1-fix_length)*transpos + fix_length*normalize(transpos);
    float s = 0.2*scale;
    bool bx = (transpos.x > 1.E-6)||(transpos.x < -1.E-6);
    bool by = (transpos.y > 1.E-6)||(transpos.y < -1.E-6);
    bool bz = (transpos.z > 1.E-6)||(transpos.z < -1.E-6);

    vertex.e1 = by ? gl_ModelViewProjectionMatrix*(s*vec4(transpos.y,-transpos.x,0.,0.)) :
	             gl_ModelViewProjectionMatrix*(s*vec4(transpos.z,0.,-transpos.x,0.));
    vertex.e2 = bx ? gl_ModelViewProjectionMatrix*(s*vec4(transpos.z,0.,-transpos.x,0.)) :
	             gl_ModelViewProjectionMatrix*(s*vec4(0.,transpos.z,-transpos.y,0.));
    vertex.normal = vec4(gl_NormalMatrix * gl_Normal,0.);
    vertex.e3 = gl_ModelViewProjectionMatrix*(scale*vec4(transpos.xyz,0.));
}
