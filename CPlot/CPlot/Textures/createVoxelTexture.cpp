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
void addNoise(unsigned char* array, int n, 
              double x, double y, double z, double r);

//=============================================================================
// Cree une texture 3D pour le voxel buffer
//=============================================================================
int Data::createVoxelTexture()
{
#ifdef __SHADERS__
  if (glewIsSupported("GL_EXT_texture3D") != 0)
  {
    glGenTextures(1, &_texVoxelBuffer);

    GLenum target = GL_TEXTURE_3D;
    GLenum filter = GL_LINEAR;
    GLenum address = GL_CLAMP_TO_BORDER;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(target, _texVoxelBuffer);

    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, filter);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, filter);

    glTexParameteri(target, GL_TEXTURE_WRAP_S, address);
    glTexParameteri(target, GL_TEXTURE_WRAP_T, address);
    glTexParameteri(target, GL_TEXTURE_WRAP_R, address);
    
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  }
#endif
  return 1;
}

//=============================================================================
/*
  Calcul le voxelArray pour la zone.
  On insere les points + bruit
  On insere les centres
  Puis des points ponderes sur les faces.

  D'autres algo pourraient etre testes:
  - le type d'insertion (pts, elts...)
  - le bruit ajoute (0, perlin absolu, perlin suivant les normales)
  - le champ (ou pas).
*/
//=============================================================================
void Data::voxelize(UnstructZone& zn, UnstructZone& z)
{
  int n = _voxelBufferSize;
  int n2 = n*n;

  unsigned char* voxelArray = new unsigned char[n2*n];

  double cx, cy, cz;
  double xf, yf, zf;
  int tx, ty, tz;
  double xmin = zn.xmin;
  double ymin = zn.ymin;
  double zmin = zn.zmin;
  double dx = n*1./(zn.xmax-zn.xmin+1.e-12);
  double dy = n*1./(zn.ymax-zn.ymin+1.e-12);
  double dz = n*1./(zn.zmax-zn.zmin+1.e-12);
  double centerx=0, centery=0, centerz=0;
  
  memset(voxelArray, 0, n2*n*sizeof(unsigned char));

  // Pyroclatic debug noise
  //centerx = n*0.5; centery = n*0.5; centerz = n*0.5;
  //addNoise(voxelArray, n, centerx, centery, centerz, z.shaderParam1*0.05);
  //zn._voxelArray = voxelArray;
  //return;
  // end debug

  // Ajout des noeuds
  for (E_Int i = 0; i < z.npts; i++)
  {
    // Coord du pts
    cx = z.x[i]; cy = z.y[i]; cz = z.z[i];
    // Coord dans le repere de la texture
    xf = (cx-xmin)*dx; yf = (cy-ymin)*dy; zf = (cz-zmin)*dz;
    centerx = centerx+xf; centery = centery+yf; centerz = centerz+zf;
    tx = int(xf); ty = int(yf); tz = int(zf);
    voxelArray[tx + n*ty + n2*tz] = 255;
    addNoise(voxelArray, n, xf, yf, zf, z.shaderParam1*0.01);
  } 
 
  // Center noise
  centerx = centerx / z.npts; centery = centery / z.npts; centerz = centerz / z.npts;
  //addNoise(voxelArray, n, centerx, centery, centerz, z.shaderParam1*0.2);

  // Ajout de pts supplementaires par elements
  int ne = z.ne;
  int ind;
  int eltSize = z.eltSize;
  double eltSizei = 1./eltSize;
  double eltSize1i = 1./(eltSize+1.);
  int* connect = z.connect;
  
  // Centre de l'element
  for (E_Int i = 0; i < ne; i++)
  {
    // Centre de l'element
    cx = 0; cy = 0; cz = 0;
    
    for (E_Int np = 0; np < eltSize; np++)
    {
      ind = connect[i+ne*np]-1;
      cx = cx + z.x[ind]; cy = cy + z.y[ind]; cz = cz + z.z[ind]; 
    }
    cx = cx * eltSizei; cy = cy * eltSizei; cz = cz * eltSizei;  
    // Coord dans le repere de la texture
    xf = (cx-xmin)*dx; yf = (cy-ymin)*dy; zf = (cz-zmin)*dz;
    tx = int(xf); ty = int(yf); tz = int(zf);
    voxelArray[tx + n*ty + n2*tz] = 255;
    //addNoise(voxelArray, n, xf, yf, zf, z.shaderParam1 / z.npts);
  }

  // eltSize point sur la face
  double* w = new double[eltSize];
  for (E_Int i = 0; i < ne; i++)
  {
    for (E_Int j = 0; j < eltSize; j++)
    {
      for (E_Int k = 0; k < eltSize; k++) w[k] = eltSize1i;
      w[j] = 2.*eltSize1i;

      cx = 0; cy = 0; cz = 0;
      for (E_Int np = 0; np < eltSize; np++)
      {
        ind = connect[i+ne*np]-1;
        cx = cx + w[np]*z.x[ind]; 
        cy = cy + w[np]*z.y[ind]; 
        cz = cz + w[np]*z.z[ind]; 
      }
      // Coord dans le repere de la texture
      xf = (cx-xmin)*dx; yf = (cy-ymin)*dy; zf = (cz-zmin)*dz;
      tx = int(xf); ty = int(yf); tz = int(zf);
      voxelArray[tx + n*ty + n2*tz] = 255;
      //addNoise(voxelArray, n, xf, yf, zf, z.shaderParam1 / z.npts);
    }
  }
  delete [] w;

  zn._voxelArray = voxelArray;
}

//=============================================================================
void Data::voxelize(StructZone& zn, StructZone& z)
{
  int n = _voxelBufferSize;
  int n2 = n*n;

  unsigned char* voxelArray = new unsigned char[n2*n];

  double cx, cy, cz;
  double xf, yf, zf;
  int tx, ty, tz;
  double xmin = zn.xmin;
  double ymin = zn.ymin;
  double zmin = zn.zmin;
  double dx = n*1./(zn.xmax-zn.xmin+1.e-12);
  double dy = n*1./(zn.ymax-zn.ymin+1.e-12);
  double dz = n*1./(zn.zmax-zn.zmin+1.e-12);
  double centerx=0, centery=0, centerz=0;
  
  memset(voxelArray, 0, n2*n*sizeof(unsigned char));

  // Pyroclatic debug noise
  //centerx = n*0.5; centery = n*0.5; centerz = n*0.5;
  //addNoise(voxelArray, n, centerx, centery, centerz, z.shaderParam1*0.05);
  //zn._voxelArray = voxelArray;
  //return;
  // end debug

  // Ajout des noeuds
  for (E_Int i = 0; i < z.npts; i++)
  {
    // Coord du pts
    cx = z.x[i]; cy = z.y[i]; cz = z.z[i];
    // Coord dans le repere de la texture
    xf = (cx-xmin)*dx; yf = (cy-ymin)*dy; zf = (cz-zmin)*dz;
    centerx = centerx+xf; centery = centery+yf; centerz = centerz+zf;
    tx = int(xf); ty = int(yf); tz = int(zf);
    voxelArray[tx + n*ty + n2*tz] = 255;
    addNoise(voxelArray, n, xf, yf, zf, z.shaderParam1*0.01);
  } 
 
  // Center noise
  centerx = centerx / z.npts; centery = centery / z.npts; centerz = centerz / z.npts;
  //addNoise(voxelArray, n, centerx, centery, centerz, z.shaderParam1*0.2);

  zn._voxelArray = voxelArray;
}

//=============================================================================
// Ajoute un bruit de perlin a array, centre en (x,y,z) et de rayon r
// x,y,z sont des coord dans le repere de la texture (0-n-1)
// n est la taille array (n*n*n)
//=============================================================================
void addNoise(unsigned char* array, int n, 
              double x, double y, double z, double r)
{
  float frequency = 3.0f / n;
  K_NOISE::PDS data;
  K_NOISE::initPerlinNoise(n, data);

  double d;
  double dxf, dyf, dzf;
  float off;
  double onen = 1./(n-1.);
  int n2 = n*n;

  double delta = 10*r*n;
  int imin = K_FUNC::E_max(int(x-delta),0);
  int imax = K_FUNC::E_min(int(x+delta),n);
  int jmin = K_FUNC::E_max(int(y-delta),0);
  int jmax = K_FUNC::E_min(int(y+delta),n);
  int kmin = K_FUNC::E_max(int(z-delta),0);
  int kmax = K_FUNC::E_min(int(z+delta),n);

  for (int i = imin; i < imax; i++)
  {
    for (int j = jmin; j < jmax; j++)
    {
      for (int k = kmin; k < kmax; k++)
      {
        dxf = x-i; dyf = y-j; dzf = z-k;
        off = fabsf(K_NOISE::perlinNoise3D(i*frequency,
                                           j*frequency,
                                           k*frequency,
                                           5, 6, 3, data)); // 5, 6, 3
        
        d = sqrtf(dxf*dxf+dyf*dyf+dzf*dzf)*onen;
        //off = 0; // Annule le bruit
        if (fabs(d-off) < r) array[i + n*j + n2*k] = 255;
      }
    }
  }
}
