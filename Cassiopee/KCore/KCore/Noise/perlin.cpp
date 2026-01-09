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

#include "Noise/noise.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define N 0x1000
#define NP 12   // 2^N
#define NM 0xfff

#define s_curve(t) (t * t * (3. - 2. * t))
#define lerp(t, a, b) (a + t * (b - a))
#define setup(i, b0, b1, r0, r1)\
        t = vec[i] + N;\
        b0 = ((int)t) & data.BM;\
        b1 = (b0+1) & data.BM;\
        r0 = t - (int)t;\
        r1 = r0 - 1.;
#define at2(rx, ry) (rx * q[0] + ry * q[1])
#define at3(rx, ry, rz) (rx * q[0] + ry * q[1] + rz * q[2])

//=============================================================================
void K_NOISE::normalize2(double* v)
{
  double s, si;
  s = sqrt(v[0]*v[0] + v[1]*v[1]);
  si = 1./s;
  v[0] = v[0] * si;
  v[1] = v[1] * si;
}

//=============================================================================
void K_NOISE::normalize3(double* v)
{
  double s, si;
  s = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  si = 1./s;
  v[0] = v[0] * si;
  v[1] = v[1] * si;
  v[2] = v[2] * si;
}

//=============================================================================
void K_NOISE::initPerlinNoise(int frequency, K_NOISE::PDS& data)
{
  data.start = 1;
  data.B = frequency;
  data.BM = frequency-1;
}

//=============================================================================
void K_NOISE::initNoise(K_NOISE::PDS& data)
{
  int i, j, k;
  
  srand(30757);
  int B = data.B;

  for (i = 0; i < B; i++)
  {
    data.p[i] = i;
    data.g1[i] = (double)((rand() % (B + B)) - B) / B;
    
    for (j = 0; j < 2; j++)
      data.g2[i][j] = (double)((rand() % (B + B)) - B) / B;
    normalize2(data.g2[i]);

    for (j = 0; j < 3; j++)
      data.g3[i][j] = (double)((rand() % (B + B)) - B) / B;
    normalize3(data.g3[i]);
  }

  while (--i)
  {
    k = data.p[i];
    data.p[i] = data.p[j = rand() % B];
    data.p[j] = k;
  }

  for (i = 0; i < B + 2; i++)
  {
    data.p[B + i] = data.p[i];
    data.g1[B + i] = data.g1[i];
    for (j = 0; j < 2; j++)
      data.g2[B + i][j] = data.g2[i][j];
    for (j = 0; j < 3; j++)
      data.g3[B + i][j] = data.g3[i][j];
  }
}

//=============================================================================
double K_NOISE::noise1(double arg, K_NOISE::PDS& data)
{
  int bx0, bx1;
  double rx0, rx1, sx, t, u, v, vec[1];
  
  vec[0] = arg;
  if (data.start)
  {
    data.start = 0; initNoise(data);
  }
  
  setup(0, bx0, bx1, rx0, rx1);
  
  sx = s_curve(rx0);
  u = rx0 * data.g1[data.p[bx0]];
  v = rx1 * data.g1[data.p[bx1]];

  return(lerp(sx, u, v));
}

//=============================================================================
double K_NOISE::noise2(double vec[2], K_NOISE::PDS& data)
{
  int bx0, bx1, by0, by1, b00, b10, b01, b11;
  double rx0, rx1, ry0, ry1, *q, sx, sy, a, b, t, u, v;
  int i, j;
  
  if (data.start)
  {
    data.start = 0; initNoise(data);
  }

  setup(0, bx0, bx1, rx0, rx1);
  setup(1, by0, by1, ry0, ry1);
  
  i = data.p[bx0];
  j = data.p[bx1];
  
  b00 = data.p[i + by0];
  b10 = data.p[j + by0];
  b01 = data.p[i + by1];
  b11 = data.p[j + by1];

  sx = s_curve(rx0);
  sy = s_curve(ry0);
  
  q = data.g2[b00]; u = at2(rx0, ry0);
  q = data.g2[b10]; v = at2(rx1, ry0);
  a = lerp(sx, u, v);
  
  q = data.g2[b01]; u = at2(rx0, ry1);
  q = data.g2[b11]; v = at2(rx1, ry1);
  b = lerp(sx, u, v);
  
  return lerp(sy, a, b);
}

//=============================================================================
double K_NOISE::noise3(double vec[3], K_NOISE::PDS& data)
{
  int bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11;
  double rx0, rx1, ry0, ry1, rz0, rz1, *q, sy, sz, a, b, c, d, t, u, v;
  int i, j;
  
  if (data.start)
  {
    data.start = 0; initNoise(data);
  }

  setup(0, bx0, bx1, rx0, rx1);
  setup(1, by0, by1, ry0, ry1);
  setup(2, bz0, bz1, rz0, rz1);
  
  i = data.p[bx0];
  j = data.p[bx1];

  b00 = data.p[i + by0];
  b10 = data.p[j + by0];
  b01 = data.p[i + by1];
  b11 = data.p[j + by1];

  t  = s_curve(rx0);
  sy = s_curve(ry0);
  sz = s_curve(rz0);
  
  q = data.g3[b00 + bz0]; u = at3(rx0, ry0, rz0);
  q = data.g3[b10 + bz0]; v = at3(rx1, ry0, rz0);
  a = lerp(t, u, v);
  
  q = data.g3[b01 + bz0]; u = at3(rx0, ry1, rz0);
  q = data.g3[b11 + bz0]; v = at3(rx1, ry1, rz0);
  b = lerp(t, u, v);
  
  c = lerp(sy, a, b);

  q = data.g3[b00 + bz1]; u = at3(rx0, ry0, rz1);
  q = data.g3[b10 + bz1]; v = at3(rx1, ry0, rz1);
  a = lerp(t, u, v);

  q = data.g3[b01 + bz1]; u = at3(rx0, ry1, rz1);
  q = data.g3[b11 + bz1]; v = at3(rx1, ry1, rz1);
  b = lerp(t, u, v);
  
  d = lerp(sy, a, b);
  
  return lerp(sz, c, d);
}

//=============================================================================
double K_NOISE::perlinNoise1D(double x, double alpha, double beta, int n, 
                              K_NOISE::PDS& data)
{
  int i;
  double val, sum=0;
  double p, scale=1;
  
  p = x;
  for (i = 0; i < n; i++)
  {
    val = noise1(p, data);
    sum += val / scale;
    scale *= alpha;
    p *= beta;
  }
  return (sum);
}

//=============================================================================
double K_NOISE::perlinNoise2D(double x, double y, double alpha, double beta, 
                              int n, K_NOISE::PDS& data)
{
  int i;
  double val, sum=0;
  double p[2], scale=1;
  
  p[0] = x;
  p[1] = y;
  for (i = 0; i < n; i++)
  {
    val = noise2(p, data);
    sum += val / scale;
    scale *= alpha;
    p[0] *= beta;
    p[1] *= beta;
  }
  return(sum);
}

//=============================================================================
double K_NOISE::perlinNoise3D(double x, double y, double z, double alpha,
                              double beta, int n, K_NOISE::PDS& data)
{
  int i;
  double val, sum=0;
  double p[3], scale=1.;
  
  p[0] = x;
  p[1] = y;
  p[2] = z;
  for (i = 0; i < n; i++)
  {
    val = noise3(p, data);
    sum += val / scale;
    scale *= alpha;
    p[0] *= beta;
    p[1] *= beta;
    p[2] *= beta;
  }
  return (sum);
}
