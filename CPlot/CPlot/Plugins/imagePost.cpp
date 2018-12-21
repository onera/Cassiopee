/*    
    Copyright 2013-2019 Onera.

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
// supersample X factor: moyenne des pixels
// IN: im1 : image 1 RGB factor * w x factor * h (deja alloue)
// OUT: im2: image 2 RGB w x h (deja alloue)
//=============================================================================
void Data::superSample(int w, int h, char* im1, char* im2, int factor)
{
  int w3=3*w; int w6=factor*3*w;
  int i2, j2; 
  //int rmoy, gmoy, bmoy;
  uint8_t r1, g1, b1;
  int r, g, b;

  double mul = 1./(factor*factor*1.);
  double co = 1./255.;
  double r2, g2, b2;

  double Y, U, V, Y1, U1, V1;

  // moyenne sur YUV
  for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
    {
      i2 = factor*i; j2 = factor*j;
      r = 0; g = 0; b = 0;
      Y = 0.; U = 0.; V = 0.;
      for (int q = 0; q < factor; q++)
        for (int p = 0; p < factor; p++)
        {
          r1 = (uint8_t)im1[3*(i2+p) + (j2+q)*w6];
          g1 = (uint8_t)im1[3*(i2+p) + (j2+q)*w6+1];
          b1 = (uint8_t)im1[3*(i2+p) + (j2+q)*w6+2];
          
          r2 = (double)r1 * co;
          g2 = (double)g1 * co;
          b2 = (double)b1 * co;
        
          Y1 = 0.299*r2 + 0.587*g2 + 0.114*b2;
          U1 = -0.14713*r2 - 0.28886*g2 + 0.436*b2;
          V1 = 0.615*r2 - 0.51499*g2 - 0.10001*b2;
          
          //Y = std::max(Y, Y1); // on garde le max de brightness des subpixels
          Y += Y1;
          U += U1;
          V += V1;
        }
      
      Y = Y*mul; 
      //Y = std::min(1.5*Y, 1.); // BOOST
      U = U*mul; V = V*mul;

      r2 = Y + 1.13983*V;
      g2 = Y -0.39465*U - 0.58060*V;
      b2 = Y + 2.03211*U;

      r2 = std::min(r2, 1.); 
      r2 = std::max(r2, 0.); 
      g2 = std::min(g2, 1.);
      g2 = std::max(g2, 0.);
      b2 = std::min(b2, 1.); 
      b2 = std::max(b2, 0.); 

      r2 = r2*255.; g2 = g2*255.; b2 = b2*255.;
      r = (int)r2; g = (int)g2; b = (int)b2;
      
      r1 = (uint8_t)r;  g1 = (uint8_t)g;  b1 = (uint8_t)b;  
      
      //printf("%d %d %d ; ", r1,g1,b1);
      im2[3*i+j*w3] = r1;
      im2[3*i+j*w3+1] = g1;
      im2[3*i+j*w3+2] = b1;
    }

  // Filtre lineaire
  /*
  double wg, w1, w2, wa;
  for (int j = 1; j < h-1; j++)
    for (int i = 1; i < w-1; i++)
    {
      i2 = factor*i; j2 = factor*j;
      r = 0.; g = 0.; b = 0.; wa = 0;
      for (int q = -1; q < 3; q++)
        for (int p = -1; p < 3; p++)
        {
          r1 = (uint8_t)im1[3*(i2+p) + (j2+q)*w6];
          g1 = (uint8_t)im1[3*(i2+p) + (j2+q)*w6+1];
          b1 = (uint8_t)im1[3*(i2+p) + (j2+q)*w6+2];
          w1 = 1./(0.5*abs(p-0.5)+1);
          w2 = 1./(0.5*abs(q-0.5)+1);
          wg = w1 * w2; wa += wg;
          r += r1*wg; g += g1*wg; b += b1*wg;
        }
      rmoy = r / wa; gmoy = g / wa; bmoy = b / wa;
      r1 = (uint8_t)rmoy; g1 = (uint8_t)gmoy; b1 = (uint8_t)bmoy;
      //printf("%d %d %d\n", r,g,b);
      im2[3*i+j*w3] = r1;
      im2[3*i+j*w3+1] = g1;
      im2[3*i+j*w3+2] = b1;
    }
  */
}

//=============================================================================
// Gaussian blur sur RGB
// IN: im1: w x h: is also modified
// IN: r: nbre d'iteration de lissage
// IN: eps: facteur de lissage
// OUT: im2: w x h: deja alloue
//=============================================================================
void Data::gaussianBlur(int w, int h, char* im1, char* im2, int r, double eps)
{
  int w3 = 3*w;
  uint8_t r1, r2, r3, r4, r5, rmoy;
  uint8_t g1, g2, g3, g4, g5, gmoy, b1, b2, b3, b4, b5, bmoy;
  int it = 0;

  while (it < r)
  {
    for (int i = 0; i < w*h*3; i++) im2[i] = im1[i];
    
    for (int j = 1; j < h-1; j++)
      for (int i = 1; i < w-1; i++)
      {
        r1 = (uint8_t)im2[3*i + j*w3];
        r2 = (uint8_t)im2[3*(i-1) + j*w3];
        r3 = (uint8_t)im2[3*(i+1) + j*w3];
        r4 = (uint8_t)im2[3*i + (j-1)*w3];
        r5 = (uint8_t)im2[3*i + (j+1)*w3];
        rmoy = (uint8_t)(r1 + eps*(r2+r3+r4+r5-4*r1));
        im1[3*i+j*w3] = (char)rmoy;

        g1 = (uint8_t)im2[3*i + j*w3+1];
        g2 = (uint8_t)im2[3*(i-1) + j*w3+1];
        g3 = (uint8_t)im2[3*(i+1) + j*w3+1];
        g4 = (uint8_t)im2[3*i + (j-1)*w3+1];
        g5 = (uint8_t)im2[3*i + (j+1)*w3+1];
        gmoy = (uint8_t)(g1 + eps*(g2+g3+g4+g5-4*g1));
        im1[3*i+j*w3+1] = (char)gmoy;

        b1 = (uint8_t)im2[3*i + j*w3+2];
        b2 = (uint8_t)im2[3*(i-1) + j*w3+2];
        b3 = (uint8_t)im2[3*(i+1) + j*w3+2];
        b4 = (uint8_t)im2[3*i + (j-1)*w3+2];
        b5 = (uint8_t)im2[3*i + (j+1)*w3+2];
        bmoy = (uint8_t)(b1 + eps*(b2+b3+b4+b5-4.*b1));
        im1[3*i+j*w3+2] = (char)bmoy;
      }
    it++;
  }
  for (int i = 0; i < w*h*3; i++) im2[i] = im1[i];
}

//=============================================================================
// MIX RGB de 2 images
// IN: im1: w x h
// IN: im2: w x h
// OUT: im1: w x h
//=============================================================================
void Data::mixImages(int w, int h, char* im1, char* im2, 
		     double alpha, double beta)
{
  int w3 = 3*w;
  int ind, ind1, ind2;
  uint8_t r1, r2;
  uint8_t g1, g2, b1, b2;
  for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
    {
      ind = 3*i + j*w3;
      ind1 = ind+1;
      ind2 = ind+2;
      r1 = (uint8_t)im1[ind];
      g1 = (uint8_t)im1[ind1];
      b1 = (uint8_t)im1[ind2];
      r2 = (uint8_t)im2[ind];
      g2 = (uint8_t)im2[ind1];
      b2 = (uint8_t)im2[ind2];
      im1[ind] = (uint8_t)(alpha*r1 + beta*r2);
      im1[ind1] = (uint8_t)(alpha*g1 + beta*g2);
      im1[ind2] = (uint8_t)(alpha*b1 + beta*b2);
    }
}

//=============================================================================
// Sharpen image (unsharp mask)
// IN: im1: w x h
// IN: amount: Sharpen amount
// IN: radius: Blur radius
// IN: threshold: thresold of sharpening
// OUT: im2 : w x h deja alloue
// Pour l'instant, n'a pas encore montre son efficacite
//=============================================================================
void Data::sharpenImage(int w, int h, char* im1, char* im2, double amount,
                        int radius, int threshold)
{
  uint8_t v1, v2;
  int val;
  char* im3 = (char*)malloc(w*h*3*sizeof(char));
  // Save image
  for (int i = 0; i < w*h*3; i++) { im3[i] = im1[i]; }

  // Blur image
  gaussianBlur(w, h, im1, im2, radius, 0.1);

  // substract
  for (int i = 0; i < w*h*3; i++)
  {
    v1 = (uint8_t)im3[i];
    v2 = (uint8_t)im2[i];
    val = v1-v2; val = std::abs(val); 
    im1[i] = (uint8_t)val;
  }

  // output  
  for (int i = 0; i < w*h*3; i++)
  {
    v1 = (uint8_t)im1[i];
    if (v1 > threshold) { 
      double vv = amount*(double)v1;
      val = (int)vv; 
      v1 = (uint8_t)im3[i];
      val += (int)v1; 
      if (val > 255) val = 255;
      im2[i] = (uint8_t)val; }
    else im2[i] = im3[i];
  }
  
  free(im3);
}
