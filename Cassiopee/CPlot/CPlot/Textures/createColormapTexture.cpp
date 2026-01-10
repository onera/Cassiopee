/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "../Data.h"

//=============================================================================
// Cree une texture 1D pour la colormap
//=============================================================================
E_Int Data::createColormapTexture()
{
  glGenTextures(1, &_texColormap);

  GLenum target = GL_TEXTURE_1D;
  GLenum filter = GL_LINEAR;
  //GLenum address = GL_CLAMP_TO_BORDER;
  GLenum address = GL_CLAMP_TO_EDGE;
  glBindTexture(target, _texColormap);
  glTexParameteri(target, GL_TEXTURE_MAG_FILTER, filter);
  glTexParameteri(target, GL_TEXTURE_MIN_FILTER, filter);
  glTexParameteri(target, GL_TEXTURE_WRAP_S, address);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  return 1;
}

// rgb->hsv: rgb: entre 0 et 1, hsv h en degres,s,v entre 0 et 1
void Data::rgb2hsv(float r, float g, float b, float& h, float& s, float& v)
{
    float min, max, delta;
    min = r < g ? r : g;
    min = min < b ? min : b;
    max = r > g ? r : g;
    max = max > b ? max : b;
    v = max;
    delta = max - min;
    if (delta < 0.00001)
    {
      s = 0.; h = 0.; 
      return;
    }
    if (max > 0.0) 
    { 
      s = (delta / max);  
    } 
    else 
    {
      // if max is 0, then r = g = b = 0              
      // s = 0, h is undefined
      s = 0.0;
      h = 0.; 
      return;
    }
    if (r >= max)                           // > is bogus, just keeps compilor happy
      h = ( g - b ) / delta;        // between yellow & magenta
    else
    if (g >= max)
      h = 2.0 + ( b - r ) / delta;  // between cyan & yellow
    else
      h = 4.0 + ( r - g ) / delta;  // between magenta & cyan

    h *= 60.0;                              // degrees
    if (h < 0.0) h += 360.0;
}

void Data::hsv2rgb(float h, float s, float v, float& r, float& g, float& b)
{
    float hh, p, q, t, ff;
    long   i;
    
    if (s <= 0.0) 
    {       // < is bogus, just shuts up warnings
      r = v; g = v; b = v;
      return;
    }
    hh = h;
    if (hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = v * (1.0 - s);
    q = v * (1.0 - (s * ff));
    t = v * (1.0 - (s * (1.0 - ff)));

    switch(i) {
    case 0:
        r = v; g = t; b = p;
        break;
    case 1:
        r = q; g = v; b = p;
        break;
    case 2:
        r = p; g = v; b = t;
        break;

    case 3:
        r = p; g = q; b = v;
        break;
    case 4:
        r = t; g = p; b = v;
        break;
    case 5:
    default:
        r = v; g = p; b = q;
        break;
    }
    return;     
}

//=============================================================================
/*
  Calcul la color map de type donne
  La stocke dans la texture colormap utilise dans le shader iso.
*/
//=============================================================================
#define clamp(a, b, c) MIN(MAX(a,b),c)
void Data::fillColormapTexture(E_Int type)
{
  float r1 = ptrState->colormapR1;
  float g1 = ptrState->colormapG1;
  float b1 = ptrState->colormapB1;
  float r2 = ptrState->colormapR2;
  float g2 = ptrState->colormapG2;
  float b2 = ptrState->colormapB2;
  float r3 = ptrState->colormapR3;
  float g3 = ptrState->colormapG3;
  float b3 = ptrState->colormapB3;
  float check = 0.;
  if (type == 1 || type == 2) check = r1+g1+b1+r2+b2+g2+r3+g3+b3;
  else if (type == 5 || type == 6)
  {
    int s = ptrState->colormapSize-1;
    float* pr = ptrState->colormapR;
    float* pg = ptrState->colormapG;
    float* pb = ptrState->colormapB;  
    check = pr[0]+pg[0]+pb[0]+pr[s]+pg[s]+pb[s]+s/2.+pr[s/2]+pg[s/2]+pb[s/2];
  }
  if (type == _texColormapType && check == _texColormapMinMax) return; 

  E_Int w = 200; // discretisation texture
  float* image = new float[w * 3];
  float f;
  float dx = 1./(w-1.);
  float r, g, b;

  switch (type)
  {
    case 0: // blue to red colormap
    {
      for (E_Int i = 0; i < w; i++)
      {
        f = i*dx;
        if (f < 0.25) 
        {
          // bleu, vert augmente
          r = 0.; g = 4.*f; b = 1.;
        }
        else if (f < 0.5)
        {
          // bleu diminue
          r = 0.; g = 1.; b = 4.*(0.5-f);
        }
        else if (f < 0.75)
        {
          // rouge augmente
          r = 4*(f-0.5); g = 1.; b = 0.;
        }
        else
        {
          // vert diminue, on finit rouge
          r = 1.; g = 4.*(1.-f); b = 0.;
        }
        //r = 4.*(f-0.5); b = 4.*(0.5-f); g = 2.-ABS(2.-4.*f);
        //r = clamp(r, 0., 1.); b = clamp(b, 0., 1.); g = clamp(g, 0., 1.);
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }
  
    case 1: // Bi-color interpolation R G B de c1, c2 
    {
      for (E_Int i = 0; i < w; i++)
      {
        f = i*dx;
        r = (1.-f)*r1+f*r2;
        g = (1.-f)*g1+f*g2;
        b = (1.-f)*b1+f*b2;
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }
  
    case 2: // Bi-color interpolation H S V de c1 a c2
    {
      float h1,s1,v1,h2,v2,s2,h,s,v,ro,go,bo;
      float delta, delta1, delta2;
      rgb2hsv(r1,g1,b1,h1,s1,v1);
      rgb2hsv(r2,g2,b2,h2,s2,v2);
      delta = fabs(h2-h1);
      delta1 = fabs(h2-h1-360.);
      delta2 = fabs(h2-h1+360.);
      if (delta1 < delta) h2 += -360.;
      else if (delta2 < delta) h2 += 360.;
      
      for (E_Int i = 0; i < w; i++)
      {
        f = i*dx;
        h = (1.-f)*h1+f*h2;
        s = (1.-f)*s1+f*s2;
        v = (1.-f)*v1+f*v2;
        if (h < 0) h += 360.;
        else if (h > 360.) h += -360.;
        hsv2rgb(h,s,v,ro,go,bo);
        image[3*i] = (float)ro; image[3*i+1] = (float)go; image[3*i+2] = (float)bo;
      }
      break;
    }

    case 7:  // diverging colormap
    {
      for (E_Int i = 0; i < w; i++)
      {
        f = i*dx;
        if (f < 0.5)
        {
          if (f < 0.25) 
          {
            if (f < 0.125)
            {
              if (f < 0.03125) 
              {r = 0.23137254902+1.12941176471*f;
                g = 0.298039215686+1.7568627451*f;
                b = 0.752941176471+239.937254902*f; }
              else if (f < 0.0625)
              {r = 0.23137254902+1.12941176471*f;
                g = 0.298039215686+1.7568627451*f;
                b = 15.6588235294-237.050980392*f;  }
              else if (f < 0.09375) 
              {r = 0.223529411765+1.25490196078*f;
                g = 0.305882352941+1.63137254902*f;
                b = 0.764705882353+1.25490196078*f; }
              else 
              {r = 0.211764705882+1.38039215686*f;
                g = 0.305882352941+1.63137254902*f;
                b = 0.776470588235+1.12941176471*f;  }
            }
            else
            {   
              if (f < 0.15625) 
              { r = 0.227450980392+1.25490196078*f;
                g = 0.321568627451+1.50588235294*f;
                b = 0.807843137255+0.878431372549*f; }
              else if (f < 0.1875) 
              { r = 0.207843137255+1.38039215686*f;
                g = 0.321568627451+1.50588235294*f;
                b = 0.827450980392+0.752941176471*f; }
              else if (f < 0.21875) 
              { r = 0.207843137255+1.38039215686*f;
                g = 0.345098039216+1.38039215686*f;
                b = 0.874509803922+0.501960784314*f;  }
              else 
              { r = 0.207843137255+1.38039215686*f;
                g = 0.345098039216+1.38039215686*f;
                b = 0.901960784314+0.376470588235*f;  } 
            }
          }
          else
          { 
            if (f < 0.375)
            {
              if (f < 0.28125) 
              { r = 0.207843137255+1.38039215686*f;
                g = 0.407843137255+1.12941176471*f;
                b = 0.964705882353+0.125490196078*f;  }
              else if (f < 0.3125) 
              { r = 0.207843137255+1.38039215686*f;
                g = 0.407843137255+1.12941176471*f;
                b = 1.0;  }
              else if (f < 0.34375) 
              { r = 0.207843137255+1.38039215686*f;
                g = 0.486274509804+0.878431372549*f;
                b = 1.07843137255-0.250980392157*f;  }
              else 
              { r = 0.250980392157+1.25490196078*f;
                g = 0.486274509804+0.878431372549*f;
                b = 1.16470588235-0.501960784314*f;  }
            }
            else
            {
              if (f < 0.40625) 
              { r = 0.250980392157+1.25490196078*f;
                g = 0.580392156863+0.627450980392*f;
                b = 1.21176470588-0.627450980392*f;  }
              else if (f < 0.4375) 
              { r = 0.250980392157+1.25490196078*f;
                g = 0.63137254902+0.501960784314*f;
                b = 1.26274509804-0.752941176471*f;  }
              else if (f < 0.46875) 
              { r = 0.305882352941+1.12941176471*f;
                g = 0.741176470588+0.250980392157*f;
                b = 1.37254901961-1.00392156863*f;  }
              else 
              { r = 0.364705882353+1.00392156863*f;
                g = 0.858823529412;
                b = 1.43137254902-1.12941176471*f;  }
            }     
          }
        }
        else
        {
          if (f < 0.75) 
          {
            if (f < 0.625)
            {
              if (f < 0.53125) 
              { r = 0.364705882353+1.00392156863*f;
                g = 0.733333333333+0.250980392157*f;
                b = 1.61960784314-1.50588235294*f; }
              else if (f < 0.5625) 
              { r = 0.43137254902+0.878431372549*f;
                g = 1.2-0.627450980392*f;
                b = 1.61960784314-1.50588235294*f;  }
              else if (f < 0.59375) 
              { r = 0.572549019608+0.627450980392*f;
                g = 1.2-0.627450980392*f;
                b = 1.61960784314-1.50588235294*f; }
              else 
              { r = 0.647058823529+0.501960784314*f;
                g = 1.34901960784-0.878431372549*f;
                b = 1.61960784314-1.50588235294*f;  }
            }
            else
            {
              if (f < 0.65625) 
              { r = 0.803921568627+0.250980392157*f;
                g = 1.42745098039-1.00392156863*f;
                b = 1.69803921569-1.63137254902*f;  }
              else if (f < 0.6875) 
              { r = 0.96862745098;
                g = 1.50980392157-1.12941176471*f;
                b = 1.61568627451-1.50588235294*f; }
              else if (f < 0.71875) 
              { r = 0.96862745098;
                g = 1.59607843137-1.25490196078*f;
                b = 1.70196078431-1.63137254902*f; }
              else 
              { r = 1.23921568627-0.376470588235*f;
                g = 1.6862745098-1.38039215686*f;
                b = 1.61176470588-1.50588235294*f;  }
            }
          }
          else
          {
            if (f < 0.875)
            {
              if (f < 0.78125) 
              { r = 1.23921568627-0.376470588235*f;
                g = 1.78039215686-1.50588235294*f;
                b = 1.61176470588-1.50588235294*f; }
              else if (f < 0.8125) 
              { r = 1.43529411765-0.627450980392*f;
                g = 1.87843137255-1.63137254902*f;
                b = 1.61176470588-1.50588235294*f; }
              else if (f < 0.84375) 
              { r = 1.63921568627-0.878431372549*f;
                g = 1.98039215686-1.7568627451*f;
                b = 1.50980392157-1.38039215686*f; }
              else 
              { r = 1.63921568627-0.878431372549*f;
                g = 2.0862745098-1.88235294118*f;
                b = 1.50980392157-1.38039215686*f;  }
            }
            else
            {
              if (f < 0.900625) 
              { r = 1.85882352941-1.12941176471*f;
                g = 2.19607843137-2.00784313725*f;
                b = 1.50980392157-1.38039215686*f; }
              else if (f < 0.9375) 
              { r = 1.97254901961-1.25490196078*f;
                g = 4.2431372549-4.26666666667*f;
                b = 1.39607843137-1.25490196078*f;  }
              else if (f < 0.96875) 
              { r = 2.09019607843-1.38039215686*f;
                g = 2.83137254902-2.76078431373*f;
                b = 1.27843137255-1.12941176471*f;  }
              else 
              { r = 2.21176470588-1.50588235294*f;
                g = 4.53333333333-4.51764705882*f;
                b = 1.27843137255-1.12941176471*f;  }
            }
          }
        }
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }
    case 3: // Tri-color interpolation R G B de c1, c3, c2 
    {
      for (E_Int i = 0; i < w/2; i++)
      {
        f = i*dx*2.;
        r = (1.-f)*r1+f*r3;
        g = (1.-f)*g1+f*g3;
        b = (1.-f)*b1+f*b3;
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      for (int i = w/2; i < w; i++)
      {
        f = i*dx*2.-1.;
        r = (1.-f)*r3+f*r2;
        g = (1.-f)*g3+f*g2;
        b = (1.-f)*b3+f*b2;
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }
    case 4: // Tri-color interpolation H S V de c1,c3 a c2
    {
      float h1,s1,v1,h2,v2,s2,h3,v3,s3,h,s,v,ro,go,bo;
      float delta, delta1, delta2,h3s;
      rgb2hsv(r1,g1,b1,h1,s1,v1);
      rgb2hsv(r2,g2,b2,h2,s2,v2);
      rgb2hsv(r3,g3,b3,h3,s3,v3);

      h3s = h3;
      delta = fabs(h3-h1);
      delta1 = fabs(h3-h1-360.);
      delta2 = fabs(h3-h1+360.);
      if (delta1 < delta) h3 += -360.;
      else if (delta2 < delta) h3 += 360.;
      //printf("%f %f %f->%f %f %f\n",r1,g1,b1,h1,s1,v1);
      //printf("%f %f %f->%f %f %f\n",r3,g3,b3,h3,s3,v3);

      for (E_Int i = 0; i < w/2; i++)
      {
        f = i*dx*2.;
        h = (1.-f)*h1+f*h3;
        s = (1.-f)*s1+f*s3;
        v = (1.-f)*v1+f*v3;
        if (h < 0) h += 360.;
        else if (h > 360.) h += -360.;
        hsv2rgb(h,s,v,ro,go,bo);
        image[3*i] = (float)ro; image[3*i+1] = (float)go; image[3*i+2] = (float)bo;
      }
      h3 = h3s;
      delta = fabs(h2-h3);
      delta1 = fabs(h2-h3-360.);
      delta2 = fabs(h2-h3+360.);
      if (delta1 < delta) h2 += -360.;
      else if (delta2 < delta) h2 += 360.;
      //printf("%f %f %f->%f %f %f\n",r3,g3,b3,h3,s3,v3);
      //printf("%f %f %f->%f %f %f\n",r2,g2,b2,h2,s2,v2);

      for (int i = w/2; i < w; i++)
      {
        f = i*dx*2.-1.;
        h = (1.-f)*h3+f*h2;
        s = (1.-f)*s3+f*s2;
        v = (1.-f)*v3+f*v2;
        if (h < 0) h += 360.;
        else if (h > 360.) h += -360.;
        hsv2rgb(h,s,v,ro,go,bo);
        image[3*i] = (float)ro; image[3*i+1] = (float)go; image[3*i+2] = (float)bo;
      }
      break;
    }

    case 5: // Multi-color interpolation R G B
    {
      E_Int size = ptrState->colormapSize;
      double dsize = 1./(size-1.);
      float* pr = ptrState->colormapR;
      float* pg = ptrState->colormapG;
      float* pb = ptrState->colormapB;
      E_Int i0, i1;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        r = (1.-f)*pr[i0]+f*pr[i1];
        g = (1.-f)*pg[i0]+f*pg[i1];
        b = (1.-f)*pb[i0]+f*pb[i1];        
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }
    case 6: // Multi-color interpolation H S V
    {
      E_Int size = ptrState->colormapSize;
      double dsize = 1./(size-1.);
      float* pr = ptrState->colormapR;
      float* pg = ptrState->colormapG;
      float* pb = ptrState->colormapB;
      E_Int i0, i1;
      float h0,s0,v0,h1,v1,s1;
      float h,s,v,ro,go,bo;
      float delta, delta1, delta2;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        rgb2hsv(pr[i0],pg[i0],pb[i0], h0,s0,v0);
        rgb2hsv(pr[i1],pg[i1],pb[i1], h1,s1,v1);
        delta = fabs(h1-h0);
        delta1 = fabs(h1-h0-360.);
        delta2 = fabs(h1-h0+360.);
        if (delta1 < delta) h1 += -360.;
        else if (delta2 < delta) h1 += 360.;
        h = (1.-f)*h0+f*h1;
        s = (1.-f)*s0+f*s1;
        v = (1.-f)*v0+f*v1;
        if (h < 0) h += 360.;
        else if (h > 360.) h += -360.;
        hsv2rgb(h,s,v,ro,go,bo);
        image[3*i] = ro; image[3*i+1] = go; image[3*i+2] = bo;
      }
      break;
    }

    case 8: // Explicitely called viridis
    {
      if (_colormapBViridis == NULL) initViridis();
      E_Int size = _colormapSizeViridis;
      double dsize = 1./(size-1.);
      float* pr = _colormapRViridis;
      float* pg = _colormapGViridis;
      float* pb = _colormapBViridis;
      E_Int i0, i1;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        r = (1.-f)*pr[i0]+f*pr[i1];
        g = (1.-f)*pg[i0]+f*pg[i1];
        b = (1.-f)*pb[i0]+f*pb[i1];        
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }

    case 9: // Explicitely called inferno
    {
      if (_colormapBInferno == NULL) initInferno();
      E_Int size = _colormapSizeInferno;
      double dsize = 1./(size-1.);
      float* pr = _colormapRInferno;
      float* pg = _colormapGInferno;
      float* pb = _colormapBInferno;
      E_Int i0, i1;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        r = (1.-f)*pr[i0]+f*pr[i1];
        g = (1.-f)*pg[i0]+f*pg[i1];
        b = (1.-f)*pb[i0]+f*pb[i1];        
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }

    case 10: // Explicitely called magma
    {
      if (_colormapBMagma == NULL) initMagma();
      E_Int size = _colormapSizeMagma;
      double dsize = 1./(size-1.);
      float* pr = _colormapRMagma;
      float* pg = _colormapGMagma;
      float* pb = _colormapBMagma;
      E_Int i0, i1;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        r = (1.-f)*pr[i0]+f*pr[i1];
        g = (1.-f)*pg[i0]+f*pg[i1];
        b = (1.-f)*pb[i0]+f*pb[i1];        
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }

    case 11: // Explicitely called plasma
    {
      if (_colormapBPlasma == NULL) initPlasma();
      E_Int size = _colormapSizePlasma;
      double dsize = 1./(size-1.);
      float* pr = _colormapRPlasma;
      float* pg = _colormapGPlasma;
      float* pb = _colormapBPlasma;
      E_Int i0, i1;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        r = (1.-f)*pr[i0]+f*pr[i1];
        g = (1.-f)*pg[i0]+f*pg[i1];
        b = (1.-f)*pb[i0]+f*pb[i1];        
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }

    case 12: // Explicitely called jet
    {
      if (_colormapBJet == NULL) initJet();
      E_Int size = _colormapSizeJet;
      double dsize = 1./(size-1.);
      float* pr = _colormapRJet;
      float* pg = _colormapGJet;
      float* pb = _colormapBJet;
      E_Int i0, i1;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        r = (1.-f)*pr[i0]+f*pr[i1];
        g = (1.-f)*pg[i0]+f*pg[i1];
        b = (1.-f)*pb[i0]+f*pb[i1];        
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }

    case 13: // Explicitely called greys
    {
      if (_colormapBGreys == NULL) initGreys();
      E_Int size = _colormapSizeGreys;
      double dsize = 1./(size-1.);
      float* pr = _colormapRGreys;
      float* pg = _colormapGGreys;
      float* pb = _colormapBGreys;
      E_Int i0, i1;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        r = (1.-f)*pr[i0]+f*pr[i1];
        g = (1.-f)*pg[i0]+f*pg[i1];
        b = (1.-f)*pb[i0]+f*pb[i1];        
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }

    case 14: // Explicitely called nice blue
    {
      r1 = 0.; g1 = 0.; b1 = 0.;
      r2 = 1.; g2 = 1.; b2 = 1.;
      r3 = 0.; g3 = 0.38; b3 = 0.647;
      for (E_Int i = 0; i < w/2; i++)
      {
        f = i*dx*2.;
        r = (1.-f)*r1+f*r3;
        g = (1.-f)*g1+f*g3;
        b = (1.-f)*b1+f*b3;
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      for (E_Int i = w/2; i < w; i++)
      {
        f = i*dx*2.-1.;
        r = (1.-f)*r3+f*r2;
        g = (1.-f)*g3+f*g2;
        b = (1.-f)*b3+f*b2;
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }

    case 15: // Explicitely called greens
    {
      if (_colormapBGreens == NULL) initGreens();
      E_Int size = _colormapSizeGreens;
      double dsize = 1./(size-1.);
      float* pr = _colormapRGreens;
      float* pg = _colormapGGreens;
      float* pb = _colormapBGreens;
      E_Int i0, i1;

      for (E_Int i = 0; i < w; i++)
      {
        i0 = (size-1.)*i*dx;
        i1 = MIN(i0+1, size-1);
        f = (i*dx-i0*dsize)/dsize;
        r = (1.-f)*pr[i0]+f*pr[i1];
        g = (1.-f)*pg[i0]+f*pg[i1];
        b = (1.-f)*pb[i0]+f*pb[i1];        
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }

    default: // invalid colormap 
    {
      for (E_Int i = 0; i < w; i++)
      {
        f = i*dx;
        r = 0.; g = 0.; b = 0.;
        image[3*i] = r; image[3*i+1] = g; image[3*i+2] = b;
      }
      break;
    }
  }

  // Commut
  glBindTexture(GL_TEXTURE_1D, _texColormap);
  glTexImage1D(GL_TEXTURE_1D, 0, 3, w, 0, 
               GL_RGB, GL_FLOAT, image);

  _texColormapType = type;
  if (type == 1 || type == 2)
    _texColormapMinMax = r1+g1+b1+r2+b2+g2+r3+g3+b3;
  else if (type == 5 || type == 6)
  {
    E_Int s = ptrState->colormapSize-1;
    float* pr = ptrState->colormapR;
    float* pg = ptrState->colormapG;
    float* pb = ptrState->colormapB;  
    _texColormapMinMax = pr[0]+pg[0]+pb[0]+pr[s]+pg[s]+pb[s]+s/2.+pr[s/2]+pg[s/2]+pb[s/2];
  }
  delete [] image;
}
