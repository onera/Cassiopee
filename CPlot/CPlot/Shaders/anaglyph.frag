// Fragment shader to perform Analyphic 3D conversion of two textures 
// for the left and right eyes - Monochrome
 
uniform sampler2D leftEyeTexture;
uniform sampler2D rightEyeTexture;
varying vec4 color;
vec2 vTexCoord;

void main(void)
{
      // monochrome anaglyph
      vTexCoord = gl_TexCoord[0].xy;
      vec4 leftFrag = texture2D(leftEyeTexture, vTexCoord);
      vec4 rightFrag = texture2D(rightEyeTexture, vTexCoord);
      gl_FragColor.r = leftFrag.r*0.299 + leftFrag.g*0.587 + leftFrag.b*0.114;
      gl_FragColor.g = rightFrag.r*0.299 + rightFrag.g*0.587 + rightFrag.b*0.114;
      gl_FragColor.b = rightFrag.r*0.299 + rightFrag.g*0.587 + rightFrag.b*0.114;
      gl_FragColor.a = color.a;
}