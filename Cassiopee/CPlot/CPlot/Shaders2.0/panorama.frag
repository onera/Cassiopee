// panorama - 360 texture from 6 cube textures of a cube map

uniform sampler2D cube_left;
uniform sampler2D cube_right;
uniform sampler2D cube_bottom;
uniform sampler2D cube_top;
uniform sampler2D cube_back;
uniform sampler2D cube_front;

#define M_PI 3.1415926535897932384626433832795

void main(void) 
{
    vec2 texcoord;
    texcoord = gl_TexCoord[0].xy;
    vec4 color;

    float theta = -M_PI + texcoord.x * 2 * M_PI; // between -pi and pi
    float phi = -M_PI/2. + texcoord.y * M_PI; // between -pi/2 and pi/2

    float x = cos(phi) * sin(theta);
    float y = sin(phi);
    float z = cos(phi) * cos(theta);

    float scale;
    vec2 px;
    
    if (abs(x) >= abs(y) && abs(x) >= abs(z)) 
    {
        if (x < 0.0) 
        {
            scale = -1.0 / x;
            px.x = ( z*scale + 1.0) / 2.0;
            px.y = ( y*scale + 1.0) / 2.0;
            color = texture2D(cube_left, px);
        }
        else 
        {
            scale = 1.0 / x;
            px.x = (-z*scale + 1.0) / 2.0;
            px.y = ( y*scale + 1.0) / 2.0;
            color = texture2D(cube_right, px);
        }
    }
    else if (abs(y) >= abs(z))
    {
        if (y < 0.0) 
        {
            scale = -1.0 / y;
            px.x = ( x*scale + 1.0) / 2.0;
            px.y = ( z*scale + 1.0) / 2.0;  
            color = texture2D(cube_top, px);
        }
        else 
        {
            scale = 1.0 / y;
            px.x = ( x*scale + 1.0) / 2.0;
            px.y = (-z*scale + 1.0) / 2.0;
            color = texture2D(cube_bottom, px);
        }
    }
    else 
    {
        if (z < 0.0) 
        {
            scale = -1.0 / z;
            px.x = (-x*scale + 1.0) / 2.0;
            px.y = ( y*scale + 1.0) / 2.0;
            color = texture2D(cube_back, px);
        }
        else
        {
            scale = 1.0 / z;
            px.x = ( x*scale + 1.0) / 2.0;      
            px.y = ( y*scale + 1.0) / 2.0;
            color = texture2D(cube_front, px);
        }
    }
    gl_FragColor = color;
}
