/* XRay shader */ 
varying vec3 N;
varying vec3 I;
varying vec4 color;

uniform float EdgeFalloff;

void main()
{
    float opacity = dot(normalize(N), normalize(I));
    opacity = abs(opacity);
    opacity = 1.0 - pow(opacity, EdgeFalloff);
    opacity = clamp(opacity, 0., 1.);
    
    gl_FragColor = color;
    gl_FragColor.a = opacity;
}
