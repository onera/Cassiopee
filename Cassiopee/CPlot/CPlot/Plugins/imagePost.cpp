/*    
    Copyright 2013-2025 Onera.

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

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

//=============================================================================
// supersample X factor: moyenne des pixels
// IN: im1 : image 1 RGB factor * w x factor * h (deja alloue)
// OUT: im2: image 2 RGB w x h (deja alloue)
//=============================================================================
void Data::superSample(E_Int w, E_Int h, char* im1, char* im2, E_Int factor)
{
  E_Int w3=3*w; E_Int w6=factor*3*w;
  E_Int i2, j2; 
  //int rmoy, gmoy, bmoy;
  uint8_t r1, g1, b1;
  int r, g, b;

  double mul = 1./(factor*factor*1.);
  double co = 1./255.;
  double r2, g2, b2;

  double Y, U, V, Y1, U1, V1;

  // moyenne sur YUV
  for (E_Int j = 0; j < h; j++)
    for (E_Int i = 0; i < w; i++)
    {
      i2 = factor*i; j2 = factor*j;
      r = 0; g = 0; b = 0;
      Y = 0.; U = 0.; V = 0.;
      for (E_Int q = 0; q < factor; q++)
        for (E_Int p = 0; p < factor; p++)
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
      
      r1 = (uint8_t)r; g1 = (uint8_t)g; b1 = (uint8_t)b;  
      
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
void Data::gaussianBlur(E_Int w, E_Int h, char* im1, char* im2, E_Int r, double eps)
{
  E_Int w3 = 3*w;
  uint8_t r1, r2, r3, r4, r5, rmoy;
  uint8_t g1, g2, g3, g4, g5, gmoy, b1, b2, b3, b4, b5, bmoy;
  E_Int it = 0;

  while (it < r)
  {
    for (E_Int i = 0; i < w*h*3; i++) im2[i] = im1[i];
    
    for (E_Int j = 1; j < h-1; j++)
      for (E_Int i = 1; i < w-1; i++)
      {
        r1 = (uint8_t)im2[3*i + j*w3];
        r2 = (uint8_t)im2[3*(i-1) + j*w3];
        r3 = (uint8_t)im2[3*(i+1) + j*w3];
        r4 = (uint8_t)im2[3*i + (j-1)*w3];
        r5 = (uint8_t)im2[3*i + (j+1)*w3];
        rmoy = (uint8_t)(r1 + eps*(r2+r3+r4+r5-4.*r1));
        im1[3*i+j*w3] = (char)rmoy;

        g1 = (uint8_t)im2[3*i + j*w3+1];
        g2 = (uint8_t)im2[3*(i-1) + j*w3+1];
        g3 = (uint8_t)im2[3*(i+1) + j*w3+1];
        g4 = (uint8_t)im2[3*i + (j-1)*w3+1];
        g5 = (uint8_t)im2[3*i + (j+1)*w3+1];
        gmoy = (uint8_t)(g1 + eps*(g2+g3+g4+g5-4.*g1));
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
  for (E_Int i = 0; i < w*h*3; i++) im2[i] = im1[i];
}

//=============================================================================
// MIX RGB de 2 images
// IN: im1: w x h
// IN: im2: w x h
// OUT: im1: w x h
//=============================================================================
void Data::mixImages(E_Int w, E_Int h, char* im1, char* im2, 
	double alpha, double beta)
{
  E_Int w3 = 3*w;
  E_Int ind, ind1, ind2;
  uint8_t r1, r2;
  uint8_t g1, g2, b1, b2;
  for (E_Int j = 0; j < h; j++)
    for (E_Int i = 0; i < w; i++)
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
void Data::sharpenImage(E_Int w, E_Int h, char* im1, char* im2, double amount,
                        E_Int radius, E_Int threshold)
{
  uint8_t v1, v2;
  int val; double vv;
  char* im3 = (char*)malloc(w*h*3*sizeof(char));
  // Save image
  for (E_Int i = 0; i < w*h*3; i++) { im3[i] = im1[i]; }

  // Blur image
  gaussianBlur(w, h, im1, im2, radius, 0.1);

  // substract
  for (E_Int i = 0; i < w*h*3; i++)
  {
    v1 = (uint8_t)im3[i];
    v2 = (uint8_t)im2[i];
    val = v1-v2; val = std::abs(val); 
    im1[i] = (uint8_t)val;
  }

  // output  
  for (E_Int i = 0; i < w*h*3; i++)
  {
    v1 = (uint8_t)im1[i];
    if (v1 > threshold) 
    { 
      vv = amount*(double)v1;
      val = (int)vv; 
      v1 = (uint8_t)im3[i];
      val += (int)v1; 
      if (val > 255) val = 255;
      im2[i] = (uint8_t)val; 
    }
    else im2[i] = im3[i];
  }
  
  free(im3);
}

// return rgb for x, y coordinates in (w,h)
void getRGB(E_Int w, E_Int h, char* im, double x, double y, double& r, double& g, double& b)
{
  E_Int i, j, ind, ip1, jp1;
  double co = 1./255.;
  double cfx, cfx1, cfy, cfy1;
  double rl, gl, bl;
  r = 0.0; g = 0.0; b = 0.0;

  i = E_Int(x);
  i = K_FUNC::E_max(i, 0); i = K_FUNC::E_min(i, w-1);
  j = E_Int(y);
  j = K_FUNC::E_max(j, 0); j = K_FUNC::E_min(j, h-1);
  ind = 3*i + j*3*w;
  cfx = x-i; cfx1 = 1.-cfx;
  cfy = y-j; cfy1 = 1.-cfy;

  rl = (uint8_t)im[ind];
  gl = (uint8_t)im[ind+1];
  bl = (uint8_t)im[ind+2];
  rl = rl*co; gl = gl*co; bl = bl*co;
  r += cfx1*cfy1*rl; g += cfx1*cfy1*gl; b += cfx1*cfy1*bl;

  ip1 = i+1;
  ip1 = K_FUNC::E_max(ip1, 0); ip1 = K_FUNC::E_min(ip1, w-1);
  jp1 = j+1;
  jp1 = K_FUNC::E_max(jp1, 0); jp1 = K_FUNC::E_min(jp1, h-1);
  
  ind = 3*ip1 + j*3*w;
  rl = (uint8_t)im[ind];
  gl = (uint8_t)im[ind+1];
  bl = (uint8_t)im[ind+2];
  rl = rl*co; gl = gl*co; bl = bl*co;
  r += cfx*cfy1*rl; g += cfx*cfy1*gl; b += cfx*cfy1*bl;

  ind = 3*ip1 + jp1*3*w;
  rl = (uint8_t)im[ind];
  gl = (uint8_t)im[ind+1];
  bl = (uint8_t)im[ind+2];
  rl = rl*co; gl = gl*co; bl = bl*co;
  r += cfx*cfy*rl; g += cfx*cfy*gl; b += cfx*cfy*bl;

  ind = 3*i + jp1*3*w;
  rl = (uint8_t)im[ind];
  gl = (uint8_t)im[ind+1];
  bl = (uint8_t)im[ind+2];
  rl = rl*co; gl = gl*co; bl = bl*co;
  r += cfx1*cfy*rl; g += cfx1*cfy*gl; b += cfx1*cfy*bl;
}

// return luma for x,y coordinates in (w,h)
double getLuma(E_Int w, E_Int h, char* im, double x, double y)
{
  double r, g, b, luma;
  getRGB(w, h, im, x, y, r, g, b);
  luma = 0.299*r + 0.587*g + 0.114*b;
  luma = std::sqrt(luma);
  return luma;
}

//============================================================================
// directional local blur (FXAA)
//============================================================================
// IN: im1: w x h
// OUT: im2 : w x h deja alloue
// from T. Lottes - nvidia
//=============================================================================
void Data::localBlur(E_Int w, E_Int h, char* im1, char* im2)
{
  E_Int ind;
  E_Int w3 = 3*w;
  double r, g, b;
  double lumC, lumNW, lumNE, lumSW, lumSE, threshold;
  double lumN, lumS, lumE, lumW;
  double rangeMin, rangeMax, range;
  double uvx, uvy, uv1x, uv1y, uv2x, uv2y;
  double lumEnd1, lumEnd2;
  double offsetx, offsety;
  double gradient1, gradient2, distance1, distance2, distanceFinal;
  double lumNS, lumEW, lumSWNW, lumSESW;
  double lumSENE, lumNENW;
  double lum1, lum2, edgeHorizontal, edgeVertical;
  bool isHorizontal, is1Steepest, reached1, reached2, reachedBoth;
  bool correctVariation, isDirection1;
  double gradientScaled, stepLength, lumMean;
  double subPixelOffsetFinal, finalUvx, finalUvy;
  double subPixelOffset2, subPixelOffset1, lumAverage;
  double finalOffset, edgeLength, pixelOffset;

#define EDGE_THRESHOLD_MIN 0.0312
#define EDGE_THRESHOLD_MAX 0.125
#define SUBPIXEL_QUALITY 0.75
  static double QUALITY[12] = {1., 1., 1., 1., 1., 1.5, 2., 2., 2., 2., 4., 8.};
  E_Int ITERATIONS = 12;

  for (E_Int j = 1; j < h-1; j++)
    for (E_Int i = 1; i < w-1; i++)
    {
      uvx = i*1.; uvy = j*1.;
      lumC = getLuma(w, h, im1, uvx, uvy);
      lumW = getLuma(w, h, im1, uvx-1, uvy);
      lumE = getLuma(w, h, im1, uvx+1, uvy);
      lumN = getLuma(w, h, im1, uvx, uvy+1);
      lumS = getLuma(w, h, im1, uvx, uvy-1);

      lumNE = getLuma(w, h, im1, uvx+1, uvy+1);
      lumSE = getLuma(w, h, im1, uvx+1, uvy-1);
      lumNW = getLuma(w, h, im1, uvx-1, uvy+1);
      lumSW = getLuma(w, h, im1, uvx-1, uvy-1);

      rangeMin = std::min(lumC, lumN);
      rangeMin = std::min(rangeMin, lumS);
      rangeMin = std::min(rangeMin, lumE);
      rangeMin = std::min(rangeMin, lumW);

      rangeMax = std::max(lumC, lumN);
      rangeMax = std::max(rangeMax, lumS);
      rangeMax = std::max(rangeMax, lumE);
      rangeMax = std::max(rangeMax, lumW);
      range = rangeMax - rangeMin; // delta luma

      threshold = rangeMax*EDGE_THRESHOLD_MAX;
      threshold = std::max(EDGE_THRESHOLD_MIN, threshold);
      if (range < threshold)
      {
        // no treatment
        ind = 3*i + j*w3;
        im2[ind] = im1[ind];
        im2[ind+1] = im1[ind+1];
        im2[ind+2] = im1[ind+2];
      }
      else
      {
        // edge horiz/vert detection
        lumNS = lumN + lumS;
        lumEW = lumE + lumW;

        lumSWNW = lumSW + lumNW;
        lumSESW = lumSE + lumSW;
        lumSENE = lumSE + lumNE;
        lumNENW = lumNE + lumNW;

        edgeHorizontal = std::abs(-2.0*lumE + lumSWNW) + std::abs(-2.0*lumC + lumNS )*2.0 + std::abs(-2.0*lumW + lumSENE);
        edgeVertical = std::abs(-2.0*lumN + lumNENW) + std::abs(-2.0*lumC + lumEW)*2.0 + std::abs(-2.0*lumS + lumSESW);
        
        isHorizontal = (edgeHorizontal >= edgeVertical);

        // luma values in ortho direction
        lum1 = isHorizontal ? lumS : lumW;
        lum2 = isHorizontal ? lumN : lumE;

        // Compute gradients
        gradient1 = lum1 - lumC;
        gradient2 = lum2 - lumC;

        is1Steepest = std::abs(gradient1) >= std::abs(gradient2);
        gradientScaled = 0.25*std::max(std::abs(gradient1),std::abs(gradient2));
        stepLength = 1.;

        // Average luma in ortho edge direction
        lumMean = 0.0;

        if (is1Steepest)
        {
          stepLength = -stepLength;
          lumMean = 0.5*(lum1 + lumC);
        } 
        else
        {
          lumMean = 0.5*(lum2 + lumC);
        }

        // Shift UV in the correct direction by half a pixel
        if (isHorizontal)
        {
          uvy += stepLength*0.5;
        } 
        else 
        {
          uvx += stepLength*0.5;
        }

        // Compute offset in ortho edge direction
        offsetx = isHorizontal ? 1.0 : 0.0;
        offsety = isHorizontal ? 0.0 : 1.0;
        
        uv1x = uvx - offsetx;
        uv1y = uvy - offsety;
        uv2x = uvx + offsetx;
        uv2y = uvy + offsety;
        
        lumEnd1 = getLuma(w, h, im1, uv1x, uv1y);
        lumEnd2 = getLuma(w, h, im1, uv2x, uv2y);
        lumEnd1 -= lumMean;
        lumEnd2 -= lumMean;

        // check edge end
        reached1 = std::abs(lumEnd1) >= gradientScaled;
        reached2 = std::abs(lumEnd2) >= gradientScaled;
        reachedBoth = reached1 && reached2;

        if (!reached1)
        {
          uv1x -= offsetx;
          uv1y -= offsety;
        }
        if (!reached2)
        {
          uv2x += offsetx;
          uv2y += offsety;
        }   

        if (!reachedBoth)
        {
          for (E_Int i = 0; i < ITERATIONS; i++)
          {
            if (!reached1)
            {
              lumEnd1 = getLuma(w, h, im1, uv1x, uv1y);
              lumEnd1 -= lumMean;
            }
            if (!reached2)
            {
              lumEnd2 = getLuma(w, h, im1, uv2x, uv2y);
              lumEnd2 -= lumMean;
            }
            reached1 = std::abs(lumEnd1) >= gradientScaled;
            reached2 = std::abs(lumEnd2) >= gradientScaled;
            reachedBoth = reached1 && reached2;

            if (!reached1)
            {
              uv1x -= offsetx * QUALITY[i];
              uv1y -= offsety * QUALITY[i];
            }
            if (!reached2)
            {
              uv2x += offsetx * QUALITY[i];
              uv2y += offsety * QUALITY[i];
            }

            if (reachedBoth) { break; }
          }
        }

        // Compute the distances to each extremity of the edge
        distance1 = isHorizontal ? (uvx - uv1x) : (uvy - uv1y);
        distance2 = isHorizontal ? (uv2x - uvx) : (uv2y - uvy);
        isDirection1 = distance1 < distance2;
        distanceFinal = std::min(distance1, distance2);
        edgeLength = (distance1 + distance2);
        pixelOffset = -distanceFinal / edgeLength + 0.5; // entre 0. et 0.5
        //pixelOffset = 0.25;

        // If the luma at center is smaller than at its neighbour, the delta luma at each end should be positive (same variation).
        // (in the direction of the closer side of the edge.)
        correctVariation = ((isDirection1 ? lumEnd1 : lumEnd2) < 0.0) != (lumC < lumMean);
        finalOffset = correctVariation ? pixelOffset : 0.0;
        
        // Weighted average of the luma over the 3x3 neighborhood
        lumAverage = (1.0/12.0) * (2.0 * (lumNS + lumEW) + lumSWNW + lumSENE);
        subPixelOffset1 = std::abs(lumAverage - lumC)/range;
        subPixelOffset1 = std::max(subPixelOffset1, 0.);
        subPixelOffset1 = std::min(subPixelOffset1, 1.);
      
        subPixelOffset2 = (-2.0*subPixelOffset1 + 3.0) * subPixelOffset1*subPixelOffset1;
        subPixelOffsetFinal = subPixelOffset2 * subPixelOffset2 * SUBPIXEL_QUALITY;

        finalOffset = std::max(finalOffset, subPixelOffsetFinal);
        
        // Compute the final UV
        uvx = i*1.; uvy = j*1.;
        finalUvx = uvx; finalUvy = uvy;

        //finalOffset = 0.5; // DBX
        //finalUvx += finalOffset * stepLength;
        //finalUvy += finalOffset * stepLength;
        //finalOffset = 0.5*finalOffset;
      
        if (isHorizontal)
        {
          finalUvy += finalOffset * stepLength;
        } 
        else 
        {
          finalUvx += finalOffset * stepLength;
        }

        // Read the color at the new UV coordinates, and use it.
        getRGB(w, h, im1, finalUvx, finalUvy, r, g, b);
        ind = 3*i + j*w3;
        im2[ind] = 255.*r;
        im2[ind+1] = 255.*g;
        im2[ind+2] = 255.*b;
      }
    }
}

//============================================================================
// specific post-processing applied to interlaced color buffer (3), return out
// do darken color with depth, blur with depth (dof)
//============================================================================
void Data::specPostProcess(char* in, E_Int ni, E_Int nj, float* depth, char* out)
{
  uint8_t r, g, b;
  printf("Specific post process\n");
  float dmin, dmax;
  // compute min/max of depth
  dmin = 1.e30; dmax = -1.e30;
  for (E_Int i = 0; i < ni*nj; i++)
  {
    dmin = std::min(dmin, depth[i]);
    // pix with 4.e7 are no body
    if (depth[i] < 1.e7) dmax = std::max(dmax, depth[i]);
  }
  
  printf("dmin=%g dmax=%g\n", dmin, dmax);
  E_Float dx = _view.xeye - _view.xcam;
  E_Float dy = _view.yeye - _view.ycam;
  E_Float dz = _view.zeye - _view.zcam;
  printf("dist=%g\n", std::sqrt(dx*dx+dy*dy+dz*dz));
  dmax = 0.3; // hard value for create for all 360 views
  
  // normalize depth
  for (E_Int i = 0; i < ni*nj; i++)
  {
    depth[i] = (depth[i]-dmin)/(dmax-dmin);
  }
 
  // darken far pixels
  E_Float percentage = 0.5;
  E_Float ramp1 = 0.35; E_Float ramp2 = 0.6;
  E_Float p, q, s;
  for (E_Int i = 0; i < ni*nj; i++)
  {
    if (depth[i] <= ramp1)
    {
      out[3*i] = in[3*i];
      out[3*i+1] = in[3*i+1];
      out[3*i+2] = in[3*i+2];
    }
    else if (depth[i] <= ramp2)
    {
      p = (percentage-1.)/(ramp2-ramp1);
      q = 1.-p*ramp1;
      s = p*depth[i]+q;
      r = (uint8_t)in[3*i];
      g = (uint8_t)in[3*i+1];
      b = (uint8_t)in[3*i+2];
      r = (uint8_t)(r*s);
      g = (uint8_t)(g*s);
      b = (uint8_t)(b*s);
      out[3*i] = r;
      out[3*i+1] = g;
      out[3*i+2] = b;
    }
    else
    {
      r = (uint8_t)in[3*i];
      g = (uint8_t)in[3*i+1];
      b = (uint8_t)in[3*i+2];
      r = (uint8_t)(r*percentage);
      g = (uint8_t)(g*percentage);
      b = (uint8_t)(b*percentage);
      out[3*i] = r;
      out[3*i+1] = g;
      out[3*i+2] = b;
    }
  }

  for (E_Int i = 0; i < ni*nj; i++)
  {
    in[3*i] = out[3*i];
    in[3*i+1] = out[3*i+1];
    in[3*i+2] = out[3*i+2];
  }
  
  // DBX - output depth
  /*
  for (E_Int i = 0; i < ni*nj; i++)
  {
    r = (uint8_t)(depth[i]*255.);
    out[3*i] = r;
    out[3*i+1] = r;
    out[3*i+2] = r;
    }*/

  // dof
  E_Float blurSigma = 0.8;
  E_Float sigma, sigma2, c;
  E_Int ind;
  E_Int n = 5; // max coc
  
  // blur
  E_Float rc, gc, bc;
  for (E_Int j = n; j < nj-n; j++)
  for (E_Int i = n; i < ni-n; i++)
  {
    ind = i+j*ni;
    if (depth[ind] <= ramp1) sigma = 0.;
    else if (depth[ind] <= ramp2)
    {
      p = (blurSigma)/(ramp2-ramp1);
      q = -p*ramp1;
      sigma = p*depth[ind]+q;
      sigma = 0.;
    }
    else sigma = blurSigma;

    sigma2 = sigma*sigma;

    rc = 0.; gc = 0.; bc = 0.;
    for (E_Int ki = -n; ki <= n; ki++)
    for (E_Int kj = -n; kj <= n; kj++)
    {
      if (sigma2 < 1.e-12)
      {
        if (ki == 0 && kj == 0) c = 1.;
        else c = 0.;
      }
      else c = exp( -(ki*ki+kj*kj)/(2.*sigma2) ) / (sigma2*2.*M_PI);
      r = (uint8_t)in[3*(ind+ki+kj*ni)];
      g = (uint8_t)in[3*(ind+ki+kj*ni)+1];
      b = (uint8_t)in[3*(ind+ki+kj*ni)+2];
      rc += (r*c);
      gc += (g*c);
      bc += (b*c);
    }

    out[3*ind] = (uint8_t)rc;
    out[3*ind+1] = (uint8_t)gc;
    out[3*ind+2] = (uint8_t)bc;
  }

}
