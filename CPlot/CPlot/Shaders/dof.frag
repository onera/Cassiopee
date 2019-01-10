#version 130
// Depth of field shader (post shader)
uniform sampler2D FrameBuffer;
uniform sampler2D depthMap;

uniform float focalDepth; // position de la focale
uniform float radius; // taille du rayon de blur
uniform float ext; // extension du blur
uniform float gamma; // gamma correction
uniform float sobelThreshold; // threshold for sobel

vec2 poisson0, poisson1, poisson2, poisson3, poisson4;
vec2 poisson5, poisson6, poisson7;
vec2 maxCoC = vec2(radius, 2.*radius);

vec2 pixelSizeHigh;

vec4 sobelColor = vec4(0,0,0,1);

vec4 dof(vec2 coords)
{
     // pas de blur au point focal
     // puis lineaire avant et apres sur dist en depth
     float dist = 0.02/(2.*radius);
     vec4 finalColor;
     float discRadius, discRadiusLow, centerDepth;

     centerDepth = texture2D(depthMap, coords).r;
     centerDepth = (centerDepth - focalDepth)/dist;
     centerDepth = pow(centerDepth, ext+0.0001);
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

vec4 sobel(vec2 coords)
{
     vec4 he = vec4(0.);
     he -= texture2D(FrameBuffer, vec2(coords.x - pixelSizeHigh.x, coords.y - pixelSizeHigh.y));
     he -= texture2D(FrameBuffer, vec2(coords.x - pixelSizeHigh.x, coords.y                  ))*2.;
     he -= texture2D(FrameBuffer, vec2(coords.x - pixelSizeHigh.x, coords.y + pixelSizeHigh.y));
     he += texture2D(FrameBuffer, vec2(coords.x + pixelSizeHigh.x, coords.y - pixelSizeHigh.y));
     he += texture2D(FrameBuffer, vec2(coords.x + pixelSizeHigh.x, coords.y                  ))*2.;
     he += texture2D(FrameBuffer, vec2(coords.x + pixelSizeHigh.x, coords.y + pixelSizeHigh.y));

     vec4 ve = vec4(0.);
     ve -= texture2D(FrameBuffer, vec2(coords.x - pixelSizeHigh.x, coords.y - pixelSizeHigh.y));
     ve -= texture2D(FrameBuffer, vec2(coords.x                  , coords.y - pixelSizeHigh.y))*2.;
     ve -= texture2D(FrameBuffer, vec2(coords.x + pixelSizeHigh.x, coords.y - pixelSizeHigh.y));
     ve += texture2D(FrameBuffer, vec2(coords.x - pixelSizeHigh.x, coords.y + pixelSizeHigh.y));
     ve += texture2D(FrameBuffer, vec2(coords.x                  , coords.y + pixelSizeHigh.y))*2.;
     ve += texture2D(FrameBuffer, vec2(coords.x + pixelSizeHigh.x, coords.y + pixelSizeHigh.y));

     vec3 e = sqrt(he.rgb * he.rgb + ve.rgb*ve.rgb);
     return vec4(e, 1.);
}

void main()
{
     pixelSizeHigh = 1.0/textureSize(FrameBuffer, 0);
 
	// poisson-distributed positions
     poisson0 = vec2( 0.0, 0.0);
     poisson1 = vec2( 0.527837,-0.08586);
     poisson2 = vec2(-0.040088, 0.536087);
     poisson3 = vec2(-0.670445,-0.179949);
     poisson4 = vec2(-0.419418,-0.616039);
     poisson5 = vec2( 0.440453,-0.639399);
     poisson6 = vec2(-0.757088, 0.349334);
     poisson7 = vec2( 0.574619, 0.685879);

     // color
     vec4 color1, color2, color;
     if (radius > 0) color1 = dof(gl_TexCoord[0].xy);
     else color1 = texture2D(FrameBuffer, gl_TexCoord[0].xy);

     if (sobelThreshold > 0) color2 = sobel(gl_TexCoord[0].xy);
     else color2 = vec4(0.);

     // mix colors
     if (sobelThreshold > 0)
     {
        float ct = max(color2.r, color2.g);
        ct = max(ct, color2.b);
        if (ct > sobelThreshold) color = sobelColor;
        else 
        color = mix(color1, 1.-color2, 0.3); // pastel colors
        //color = mix(color1,sobelColor,(1./sobelThreshold)*ct+1.);
     }
     else color = color1;

     // gamma correction
     vec4 res = pow(color, vec4(1.0 / gamma));
     res.a = color.a;
     
	gl_FragColor = res;
}
