// isos + gooch
#version 150 compatibility

varying vec4 color;

void main()
{ 
  color = gl_Color;
  gl_Position = ftransform();
}
