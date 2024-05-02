// Fragment shader to perform Analyphic 3D conversion of two textures 
// for the left and right eyes - Color
 
uniform sampler2D leftEyeTexture;
uniform sampler2D rightEyeTexture;
varying vec4 color;

void main(void)
{
      // Color anaglyph
      vec2 vTexCoord = gl_TexCoord[0].xy;
      vec4 leftFrag = texture2D(leftEyeTexture, vTexCoord);
      vec4 rightFrag = texture2D(rightEyeTexture, vTexCoord);

      // Standard
      //gl_FragColor.r = leftFrag.g*0.7 + leftFrag.b*0.3;
      //gl_FragColor.g = rightFrag.g*1.;
      //gl_FragColor.b = rightFrag.b*1.;
      //gl_FragColor.a = color.a;

      // Optimise par Dubois
      gl_FragColor.r = 0.4154*leftFrag.r + 0.4710*leftFrag.g + 0.1669*leftFrag.b -0.0109*rightFrag.r -0.0364* rightFrag.g -0.00360*rightFrag.b;
      gl_FragColor.g = -0.0458*leftFrag.r - 0.0484*leftFrag.g - 0.0257*leftFrag.b +0.3756*rightFrag.r +0.733* rightFrag.g +0.0111*rightFrag.b;
      gl_FragColor.b = -0.0547*leftFrag.r - 0.0615*leftFrag.g + 0.0128*leftFrag.b -0.0651*rightFrag.r -0.1287* rightFrag.g +1.2971*rightFrag.b;
      gl_FragColor.a = color.a;
}