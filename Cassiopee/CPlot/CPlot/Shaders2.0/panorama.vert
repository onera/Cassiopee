// panorama - 360 texture from 6 cube textures of cube map
varying vec4 color;
void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	color = gl_Color;
	gl_Position = ftransform();
}