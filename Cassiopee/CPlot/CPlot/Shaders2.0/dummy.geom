#version 150
layout(triangles) in;
layout (triangle_strip, max_vertices=6) out;

//layout (std140) uniformMatrices {
//  mat4 projModelViewMatrix;
//  mat3 normalMatrix;
//};

in VertexData {
  //vec2 texCoord;
  vec3 normal;
  } VertexIn[];

out VertexData {
  //vec2 texCoord;
  vec3 normal;
} VertexOut;

void main()
{
  for (int i = 0; i < gl_in.length(); i++)
  {
    // copy attributes
    gl_Position = gl_in[i].gl_Position;
    VertexOut.normal = VertexIn[i].normal;
    //VertexOut.texCoord = VertexIn[i].texCoord;
    EmitVertex();
  }
  EndPrimitive();
  for (int i = 0; i < gl_in.length(); i++)
  {
    // copy attributes
    gl_Position = gl_in[i].gl_Position+vec4(1,0,0,0);
    VertexOut.normal = VertexIn[i].normal;
    //VertexOut.texCoord = VertexIn[i].texCoord;
    EmitVertex();
  }
  EndPrimitive();
}
/*
#version 410
in vec3 position;
in vec3 normal;

outVertexData {
  vec3 normal;
} VertexOut;

void main()
{
  VertexOut.normal = normal;
  gl_Position = vec4(position, 1.);
}
*/
