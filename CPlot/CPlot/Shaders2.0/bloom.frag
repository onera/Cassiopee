// Depth of field shader
uniform sampler2D FrameBuffer;
uniform sampler2D depthMap;

uniform int windowWidth, windowHeight;
uniform float focalDepth; // position de la focale
uniform float radius; // taille du rayon de blur

vec2 poisson0, poisson1, poisson2, poisson3, poisson4;
vec2 poisson5, poisson6, poisson7;
vec2 maxCoC = vec2(radius, 2.*radius);

vec2 pixelSizeHigh;

vec4 dof(vec2 coords)
{
     // pas de blur au point focal
     // puis lineaire avant et apres sur dist en depth
     float dist = 0.02/(2.*radius);
     vec4 finalColor;
     float discRadius, discRadiusLow, centerDepth;

     centerDepth = texture2D(depthMap, coords).r;
     centerDepth = (centerDepth - focalDepth)/dist;
     centerDepth = clamp(centerDepth,-1.,1.);
     centerDepth = centerDepth*0.5 + 0.5;

     discRadius = abs(centerDepth*maxCoC.y - maxCoC.x);
     discRadiusLow = discRadius * 0.5;

     finalColor = vec4(0.0,0.0,0.0,0.0);
     float finalBlur = 0.;
     float blur, tapBlur;
    
     // -- 
     vec2 coordHigh = coords + (pixelSizeHigh*poisson0*discRadius);
     vec4 tapHigh = texture2D(FrameBuffer, coordHigh);
     float blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;
      
     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson1*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson2*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson3*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson4*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson5*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson6*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson7*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     return finalColor/finalBlur;
}

// Hue : in 0,360
// Saturation: 0-1
// L: in 0-1
vec4 HSL(vec4 color)
{
     float r = color.r;
     float g = color.g;
     float b = color.b;
     float Cmax = max(r,max(g,b));
     float Cmin = min(r,min(g,b));
     float delta = Cmax-Cmin;
     float H = 0.;
     if (delta == 0.) H = 0.;
     else if (Cmax == r) H = mod((g-b)/delta, 6.) * 60.;
     else if (Cmax == g) H = ((b-r)/delta+2.) * 60.;
     else H = ((r-g)/delta+4.) * 60.;
     float L = (Cmax+Cmin)*0.5;
     float S = 0.;
     if (delta == 0.) S = 0.;
     else S = delta/(1.-abs(2*L-1.));
     return vec4(H/360.,S,L,color.a);
}

// Color is HSL
vec4 RGB(vec4 color)
{
     float H = color.r * 360.;
     float S = color.g;
     float L = color.b;
     float C = (1. - abs(2*L-1.))*S;
     float X = C*(1.-abs(mod(H/60.,2.) -1.));
     float m = L - 0.5*C;
     float r1,g1,b1;
     if (H < 60.) { r1 = C; g1 = X; b1 = 0.; }
     else if (H < 120.) { r1 = X; g1 = C; b1 = 0.; }
     else if (H < 180.) { r1 = 0.; g1 = C; b1 = X; }
     else if (H < 240.) { r1 = 0.; g1 = X; b1 = C; }
     else if (H < 300.) { r1 = X; g1 = 0.; b1 = C; }
     else { r1 = C; g1 = 0.; b1 = X; }
     
     return vec4(r1+m,g1+m,b1+m,color.a);
}

// bloom
vec4 bloom(vec2 coords)
{
     float weight[5];
     weight[0] = 0.227027; weight[1] = 0.1945946; weight[2] = 0.1216216; weight[3] = 0.054054; weight[4] = 0.016216;

     vec4 color = texture(FrameBuffer, coords);

     float brightness = dot(color.rgb, vec3(0.2126, 0.7152, 0.0722));
     vec3 result;
     if (brightness > 1.)
     {
     vec2 texOffset = 1.0 / textureSize(FrameBuffer, 0); // gets size of single texel
     result = texture(FrameBuffer, coords).rgb * weight[0];

     for (int i = 1; i < 5; ++i)
     {
          result += texture(FrameBuffer, coords + vec2(texOffset.x * i, 0.0)).rgb * weight[i];
          result += texture(FrameBuffer, coords - vec2(texOffset.x * i, 0.0)).rgb * weight[i];
     }
     for(int i = 1; i < 5; ++i)
     {
          result += texture(FrameBuffer, coords + vec2(0.0, texOffset.y * i)).rgb * weight[i];
          result += texture(FrameBuffer, coords - vec2(0.0, texOffset.y * i)).rgb * weight[i];
     }
     
     // blending
     result = color.rgb + result;

     // tone mapping
     //const float exposure = 2.;
     //result = vec3(1.0) - exp(-result * exposure);

     // gamma correction
     const float gamma = 1.;
     result = pow(result, vec3(1.0 / gamma));
     }
     else
     result = texture(FrameBuffer, coords).rgb;
     return vec4(result, color.a);
}

void main(){

	pixelSizeHigh[0] = 1.0/float(windowWidth);
	pixelSizeHigh[1] = 1.0/float(windowHeight);

	// poisson-distributed positions
     poisson0 = vec2( 0.0, 0.0);
     poisson1 = vec2( 0.527837,-0.08586);
     poisson2 = vec2(-0.040088, 0.536087);
     poisson3 = vec2(-0.670445,-0.179949);
     poisson4 = vec2(-0.419418,-0.616039);
     poisson5 = vec2( 0.440453,-0.639399);
     poisson6 = vec2(-0.757088, 0.349334);
     poisson7 = vec2( 0.574619, 0.685879);

	vec4 col = dof(gl_TexCoord[0].xy);

     // by pass
     vec2 coord = gl_TexCoord[0].xy;
     col = texture2D(FrameBuffer, coord)+0.000001*col;
     vec4 hsl = HSL(col);
     // augmente la saturation
     hsl.g = hsl.g*0.8;
     // augmente la luminosite
     hsl.b = hsl.b*1.4;
     vec4 col2 = RGB(hsl);

     // Bloom
     vec4 col3 = bloom(coord)+0.000001*col2;
     gl_FragColor = col3;
}
