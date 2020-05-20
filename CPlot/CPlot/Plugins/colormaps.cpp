/*    
    Copyright 2013-2020 Onera.

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
#include <cmath>

//=============================================================================
// Colormaps plugins
//=============================================================================
// Attention: si on utilise les SHADERS, les colormaps utilisees sont 
// dans createColormapTextures

// Parameters:
// IN: f -> scalar field value
// OUT: r, g, b -> red, green, blue color values.

//============================================================================
/*
  Blue to red colormap (par defaut pour les isos)
*/
//============================================================================
void colBlueToRed(Data* d, double f, float* r, float* g, float* b)
{
  if (f < 0.25) 
  {
    // bleu, vert augmente
    *r = 0.;
    *g = 4.*f;
    *b = 1.;
  }
  else if (f < 0.5)
  {
    // bleu diminue
    *r = 0.;
    *g = 1.;
    *b = 4.*(0.5-f);
  }
  else if (f < 0.75)
  {
    // rouge augmente
    *r = 4*(f-0.5);
    *g = 1.;
    *b = 0.;
  }
  else
  {
    // vert diminue, on finit rouge
    *r = 1.;
    *g = 4.*(1.-f);
    *b = 0.;
  }
  //printf("%f %f %f %f\n", f, *r, *g, *b);
}

//============================================================================
/*
  Green to red colormap (par defaut pour les maillages)
  Dans cette colormap, on evite le (0.1,0.1,1.) utilise pour la selection.
*/
//============================================================================
void colGreenToRed(Data* d, double f, float* r, float* g, float* b)
{
  if (f < 0.25)
  {
    // vert, bleu augmente
    *r = 0.;
    *g = 0.9;
    *b = 4.*f;
  }
  else if (f < 0.5)
  {
    // vert diminue (evite le bleu de selection)
    *r = 0.4;
    *g = 4.*(0.5-f);
    *b = 0.9;
  }
  else if (f < 0.75)
  {
    // rouge augmente
    *r = 4.*(f-0.5);
    *g = 0.9;
    *b = 4.*(0.75-f);
  }
  else
  {
    *r = 0.9;
    *g = 4.*(f-0.75);
    *b = 0.;
  }

  //printf("%f %f %f %f\n", f, *r, *g, *b);
}

//============================================================================
/*
  Bi-color colormap with RGB interpolation
*/
//============================================================================
void col2RGB(Data* d, double f, float* r, float* g, float* b)
{
  double r1 = d->ptrState->colormapR1;
  double g1 = d->ptrState->colormapG1;
  double b1 = d->ptrState->colormapB1;
  double r2 = d->ptrState->colormapR2;
  double g2 = d->ptrState->colormapG2;
  double b2 = d->ptrState->colormapB2;
  *r = (1.-f)*r1+f*r2;
  *g = (1.-f)*g1+f*g2;
  *b = (1.-f)*b1+f*b2;
}

//============================================================================
/*
  Bi-color colormap with HSV interpolation
*/
//============================================================================
void col2HSV(Data* d, double f, float* r, float* g, float* b)
{
  double r1 = d->ptrState->colormapR1;
  double g1 = d->ptrState->colormapG1;
  double b1 = d->ptrState->colormapB1;
  double r2 = d->ptrState->colormapR2;
  double g2 = d->ptrState->colormapG2;
  double b2 = d->ptrState->colormapB2;
  double h1,s1,v1,h2,v2,s2,h,s,v,ro,go,bo;
  double delta, delta1, delta2;

  d->rgb2hsv(r1,g1,b1,h1,s1,v1);
  d->rgb2hsv(r2,g2,b2,h2,s2,v2);
  delta = fabs(h2-h1);
  delta1 = fabs(h2-h1-360.);
  delta2 = fabs(h2-h1+360.);
  if (delta1 < delta) h2 += -360.;
  else if (delta2 < delta) h2 += 360.;

  h = (1.-f)*h1+f*h2;
  s = (1.-f)*s1+f*s2;
  v = (1.-f)*v1+f*v2;
  if (h < 0) h += 360.;
  else if (h > 360.) h += -360.;
  d->hsv2rgb(h,s,v,ro,go,bo);
  *r = (float)ro;
  *g = (float)go;
  *b = (float)bo;
}

//=============================================================================
/*
  Diverging colormap
*/
//=============================================================================
void diverging(Data* d, double f, float* r, float* g, float* b)
{
  if (f < 0.5)
  {
    if (f < 0.25) 
    {
      if (f < 0.125)
      {
        if (f < 0.03125) 
        {*r = 0.23137254902+1.12941176471*f;
          *g = 0.298039215686+1.7568627451*f;
          *b = 0.752941176471+239.937254902*f; }
        else if (f < 0.0625) 
        {*r = 0.23137254902+1.12941176471*f;
          *g = 0.298039215686+1.7568627451*f;
          *b = 15.6588235294-237.050980392*f;  }
        else if (f < 0.09375) 
        {*r = 0.223529411765+1.25490196078*f;
          *g = 0.305882352941+1.63137254902*f;
          *b = 0.764705882353+1.25490196078*f; }
        else 
        {*r = 0.211764705882+1.38039215686*f;
          *g = 0.305882352941+1.63137254902*f;
          *b = 0.776470588235+1.12941176471*f;  }
      }
      else
      {   
        if (f < 0.15625) 
        { *r = 0.227450980392+1.25490196078*f;
          *g = 0.321568627451+1.50588235294*f;
          *b = 0.807843137255+0.878431372549*f; }
        else if (f < 0.1875) 
        { *r = 0.207843137255+1.38039215686*f;
          *g = 0.321568627451+1.50588235294*f;
          *b = 0.827450980392+0.752941176471*f; }
        else if (f < 0.21875) 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.345098039216+1.38039215686*f;
          *b = 0.874509803922+0.501960784314*f;  }
        else 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.345098039216+1.38039215686*f;
          *b = 0.901960784314+0.376470588235*f;  } 
      }
    }
    else
    { 
      if (f < 0.375)
      {
        if (f < 0.28125) 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.407843137255+1.12941176471*f;
          *b = 0.964705882353+0.125490196078*f;  }
        else if (f < 0.3125) 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.407843137255+1.12941176471*f;
          *b = 1.0+0.0*f;  }
        else if (f < 0.34375) 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.486274509804+0.878431372549*f;
          *b = 1.07843137255-0.250980392157*f;  }
        else 
        {*r = 0.250980392157+1.25490196078*f;
          *g = 0.486274509804+0.878431372549*f;
          *b = 1.16470588235-0.501960784314*f;  }
      }
      else
      {
        if (f < 0.40625) 
        {*r = 0.250980392157+1.25490196078*f;
          *g = 0.580392156863+0.627450980392*f;
          *b = 1.21176470588-0.627450980392*f;  }
        else if (f < 0.4375) 
        {*r = 0.250980392157+1.25490196078*f;
          *g = 0.63137254902+0.501960784314*f;
          *b = 1.26274509804-0.752941176471*f;  }
        else if (f < 0.46875) 
        {*r = 0.305882352941+1.12941176471*f;
          *g = 0.741176470588+0.250980392157*f;
          *b = 1.37254901961-1.00392156863*f;  }
        else 
        {*r = 0.364705882353+1.00392156863*f;
          *g = 0.858823529412+0.0*f;
          *b = 1.43137254902-1.12941176471*f;  }
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
        { *r = 0.364705882353+1.00392156863*f;
          *g = 0.733333333333+0.250980392157*f;
          *b = 1.61960784314-1.50588235294*f; }
        else if (f < 0.5625) 
        {*r = 0.43137254902+0.878431372549*f;
          *g = 1.2-0.627450980392*f;
          *b = 1.61960784314-1.50588235294*f;  }
        else if (f < 0.59375) 
        { *r = 0.572549019608+0.627450980392*f;
          *g = 1.2-0.627450980392*f;
          *b = 1.61960784314-1.50588235294*f; }
        else 
        {*r = 0.647058823529+0.501960784314*f;
          *g = 1.34901960784-0.878431372549*f;
          *b = 1.61960784314-1.50588235294*f;  }
      }
      else
      {
        if (f < 0.65625) 
        {*r = 0.803921568627+0.250980392157*f;
          *g = 1.42745098039-1.00392156863*f;
          *b = 1.69803921569-1.63137254902*f;  }
        else if (f < 0.6875) 
        { *r = 0.96862745098+0.0*f;
          *g = 1.50980392157-1.12941176471*f;
          *b = 1.61568627451-1.50588235294*f; }
        else if (f < 0.71875) 
        { *r = 0.96862745098+0.0*f;
          *g = 1.59607843137-1.25490196078*f;
          *b = 1.70196078431-1.63137254902*f; }
        else 
        {*r = 1.23921568627-0.376470588235*f;
          *g = 1.6862745098-1.38039215686*f;
          *b = 1.61176470588-1.50588235294*f;  }
      }
    }
    else
    {
      if (f < 0.875)
      {
        if (f < 0.78125) 
        { *r = 1.23921568627-0.376470588235*f;
          *g = 1.78039215686-1.50588235294*f;
          *b = 1.61176470588-1.50588235294*f; }
        else if (f < 0.8125) 
        {*r = 1.43529411765-0.627450980392*f;
          *g = 1.87843137255-1.63137254902*f;
          *b = 1.61176470588-1.50588235294*f; }
        else if (f < 0.84375) 
        { *r = 1.63921568627-0.878431372549*f;
          *g = 1.98039215686-1.7568627451*f;
          *b = 1.50980392157-1.38039215686*f; }
        else 
        {*r = 1.63921568627-0.878431372549*f;
          *g = 2.0862745098-1.88235294118*f;
          *b = 1.50980392157-1.38039215686*f;  }
      }
      else
      {
        if (f < 0.900625) 
        { *r = 1.85882352941-1.12941176471*f;
          *g = 2.19607843137-2.00784313725*f;
          *b = 1.50980392157-1.38039215686*f; }
        else if (f < 0.9375) 
        {*r = 1.97254901961-1.25490196078*f;
          *g = 4.2431372549-4.26666666667*f;
          *b = 1.39607843137-1.25490196078*f;  }
        else if (f < 0.96875) 
        {*r = 2.09019607843-1.38039215686*f;
          *g = 2.83137254902-2.76078431373*f;
          *b = 1.27843137255-1.12941176471*f;  }
        else 
        {*r = 2.21176470588-1.50588235294*f;
          *g = 4.53333333333-4.51764705882*f;
          *b = 1.27843137255-1.12941176471*f;  }
      }
    }
  }
}

// Retourne une couleur complementaire
// IN: r,g,b
// OUT: ro,go,bo complementaire
void complementColor(float r, float g, float b,
                     float& ro, float& go, float& bo)
{
  // passage HSL; H entre 0 et 360,S et L entre 0 et 1.
  float Cmax = std::max(r,std::max(g,b));
  float Cmin = std::min(r,std::min(g,b));
  float delta = Cmax-Cmin;
  float H = 0.;
  if (delta == 0.) H = 0.;
  else if (Cmax == r) H = std::fmod((g-b)/delta, (float)6.) * 60.;
  else if (Cmax == g) H = (2.+(b-r)/delta) * 60.;
  else H = (4.+(r-g)/delta) * 60.;
  float L = (Cmax+Cmin)*0.5;
  float S = 0.;
  if (delta == 0.) S = 0.;
  else S = delta/(1.-std::abs(Cmax+Cmin-1.));
  //printf("HSL %f %f %f\n",H,S,L);
  
  // decalage teinte 
  H = H+180.;
  if (H > 360.) H = H-360.;
  if (L+0.2 > 1.) L = 1.-L;
  else L = L+0.2;

  // retour RGB
  float C = (1.-std::abs(2.*L-1.))*S;
  float X = C*(1.-std::abs(std::fmod(H/60.,2.) -1.));
  float m = L - 0.5*C;
  float r1,g1,b1;
  if (H >= 300.)      { r1 = C;  g1 = 0.; b1 = X; }
  else if (H >= 240.) { r1 = X;  g1 = 0.; b1 = C; }
  else if (H >= 180.) { r1 = 0.; g1 = X;  b1 = C; }
  else if (H >= 120.) { r1 = 0.; g1 = C;  b1 = X; }
  else if (H >= 60.)  { r1 = X;  g1 = C;  b1 = 0.; }
  else                { r1 = C;  g1 = X;  b1 = 0.; }

  ro = r1+m;
  go = g1+m;
  bo = b1+m;
}
//============================================================================
/*
  Tri-color colormap with RGB interpolation
*/
//============================================================================
void col3RGB(Data* d, double f, float* r, float* g, float* b)
{
  double r1 = d->ptrState->colormapR1;
  double g1 = d->ptrState->colormapG1;
  double b1 = d->ptrState->colormapB1;
  double r2 = d->ptrState->colormapR2;
  double g2 = d->ptrState->colormapG2;
  double b2 = d->ptrState->colormapB2;
  double r3 = d->ptrState->colormapR3;
  double g3 = d->ptrState->colormapG3;
  double b3 = d->ptrState->colormapB3;
  
  if (f < 0.5)
  {
    f = 2.*f;
    *r = (1.-f)*r1+f*r3;
    *g = (1.-f)*g1+f*g3;
    *b = (1.-f)*b1+f*b3;
  }
  else
  {
    
    f = 2.*f-1.;
    *r = (1.-f)*r3+f*r2;
    *g = (1.-f)*g3+f*g2;
    *b = (1.-f)*b3+f*b2; 
  }  
}

//============================================================================
/*
  Tri-color colormap with HSV interpolation
*/
//============================================================================
void col3HSV(Data* d, double f, float* r, float* g, float* b)
{
  double r1 = d->ptrState->colormapR1;
  double g1 = d->ptrState->colormapG1;
  double b1 = d->ptrState->colormapB1;
  double r2 = d->ptrState->colormapR2;
  double g2 = d->ptrState->colormapG2;
  double b2 = d->ptrState->colormapB2;
  double r3 = d->ptrState->colormapR3;
  double g3 = d->ptrState->colormapG3;
  double b3 = d->ptrState->colormapB3;
  double h1,s1,v1,h2,v2,s2,h3,v3,s3,h,s,v,ro,go,bo;
  double delta, delta1, delta2;
  d->rgb2hsv(r1,g1,b1,h1,s1,v1);
  d->rgb2hsv(r2,g2,b2,h2,s2,v2);
  d->rgb2hsv(r3,g3,b3,h3,s3,v3);
  
  if (f < 0.5)
  {
    delta = fabs(h3-h1);
    delta1 = fabs(h3-h1-360.);
    delta2 = fabs(h3-h1+360.);
    if (delta1 < delta) h3 += -360.;
    else if (delta2 < delta) h3 += 360.;
    f = 2*f;
    h = (1.-f)*h1+f*h3;
    s = (1.-f)*s1+f*s3;
    v = (1.-f)*v1+f*v3;
    if (h < 0) h += 360.;
    else if (h > 360.) h += -360.;
    d->hsv2rgb(h,s,v,ro,go,bo);
    *r = (float)ro;
    *g = (float)go;
    *b = (float)bo;
  }
  else
  {
    delta = fabs(h2-h3);
    delta1 = fabs(h2-h3-360.);
    delta2 = fabs(h2-h3+360.);
    if (delta1 < delta) h2 += -360.;
    else if (delta2 < delta) h2 += 360.;
    f = 2*f-1.;
    h = (1.-f)*h3+f*h2;
    s = (1.-f)*s3+f*s2;
    v = (1.-f)*v3+f*v2;
    if (h < 0) h += 360.;
    else if (h > 360.) h += -360.;
    d->hsv2rgb(h,s,v,ro,go,bo);
    *r = (float)ro;
    *g = (float)go;
    *b = (float)bo; 
  }
}
