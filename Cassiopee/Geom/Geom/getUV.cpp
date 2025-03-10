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

// Surface uv map

# include <string.h>
# include "geom.h"
# include "xatlas/xatlas.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// used internally
static void RandomColor(uint8_t *color)
{
  for (int i = 0; i < 3; i++)
    color[i] = uint8_t((rand() % 255 + 192) * 0.5f);
}

static void SetPixel(uint8_t *dest, int destWidth, int x, int y, const uint8_t *color)
{
  uint8_t *pixel = &dest[x * 3 + y * (destWidth * 3)];
  pixel[0] = color[0];
  pixel[1] = color[1];
  pixel[2] = color[2];
}

// https://github.com/miloyip/line/blob/master/line_bresenham.c
// License: public domain.
static void RasterizeLine(uint8_t *dest, int destWidth, const int *p1, const int *p2, const uint8_t *color)
{
  const int dx = abs(p2[0] - p1[0]), sx = p1[0] < p2[0] ? 1 : -1;
  const int dy = abs(p2[1] - p1[1]), sy = p1[1] < p2[1] ? 1 : -1;
  int err = (dx > dy ? dx : -dy) / 2;
  int current[2];
  current[0] = p1[0];
  current[1] = p1[1];
  while (SetPixel(dest, destWidth, current[0], current[1], color), current[0] != p2[0] || current[1] != p2[1])
  {
    const int e2 = err;
    if (e2 > -dx) { err -= dy; current[0] += sx; }
    if (e2 < dy) { err += dx; current[1] += sy; }
  }
}

/*
https://github.com/ssloy/tinyrenderer/wiki/Lesson-2:-Triangle-rasterization-and-back-face-culling
Copyright Dmitry V. Sokolov

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
static void RasterizeTriangle(uint8_t *dest, int destWidth, const int *t0, const int *t1, const int *t2, const uint8_t *color)
{
  if (t0[1] > t1[1]) std::swap(t0, t1);
  if (t0[1] > t2[1]) std::swap(t0, t2);
  if (t1[1] > t2[1]) std::swap(t1, t2);
  int total_height = t2[1] - t0[1];
  for (int i = 0; i < total_height; i++) {
    bool second_half = i > t1[1] - t0[1] || t1[1] == t0[1];
    int segment_height = second_half ? t2[1] - t1[1] : t1[1] - t0[1];
    float alpha = (float)i / total_height;
    float beta = (float)(i - (second_half ? t1[1] - t0[1] : 0)) / segment_height;
    int A[2], B[2];
    for (int j = 0; j < 2; j++) {
      A[j] = int(t0[j] + (t2[j] - t0[j]) * alpha);
      B[j] = int(second_half ? t1[j] + (t2[j] - t1[j]) * beta : t0[j] + (t1[j] - t0[j]) * beta);
    }
    if (A[0] > B[0]) std::swap(A, B);
    for (int j = A[0]; j <= B[0]; j++)
      SetPixel(dest, destWidth, j, t0[1] + i, color);
  }
}

// Rasterize triangle with a different color for each vertex (field)
static void RasterizeTriangle2(uint8_t *dest, int destWidth, const int *t0, const int *t1, const int *t2, 
                               uint8_t *colorv0, uint8_t *colorv1, uint8_t* colorv2)
{
  //printf("rastering triangle : u0=%d v0=%d\n", t0[0], t0[1]);
  //printf("rastering triangle : u1=%d v1=%d\n", t1[0], t1[1]);
  //printf("rastering triangle : u2=%d v2=%d\n", t2[0], t2[1]);
  
  if (t0[1] > t1[1]) { std::swap(t0, t1); std::swap(colorv0, colorv1); }
  if (t0[1] > t2[1]) { std::swap(t0, t2); std::swap(colorv0, colorv2); }
  if (t1[1] > t2[1]) { std::swap(t1, t2); std::swap(colorv1, colorv2); }

  //printf("after rastering triangle : u0=%d v0=%d\n", t0[0], t0[1]);
  //printf("after rastering triangle : u1=%d v1=%d\n", t1[0], t1[1]);
  //printf("after rastering triangle : u2=%d v2=%d\n", t2[0], t2[1]);

  int total_height = t2[1] - t0[1];
  int A[2], B[2];
  uint8_t color[3];
  E_Float u,v,r0,r1,r2,u0,u1,u2,v0,v1,v2,dd,a0,a1,a2;
  for (int i = 0; i < total_height; i++)
  {
    bool second_half = i > t1[1] - t0[1] || t1[1] == t0[1];
    int segment_height = second_half ? t2[1] - t1[1] : t1[1] - t0[1];
    float alpha = (float)i / total_height;
    float beta = (float)(i - (second_half ? t1[1] - t0[1] : 0)) / segment_height;
    for (int j = 0; j < 2; j++) 
    {
      A[j] = int(t0[j] + (t2[j] - t0[j]) * alpha);
      B[j] = int(second_half ? t1[j] + (t2[j] - t1[j]) * beta : t0[j] + (t1[j] - t0[j]) * beta);
    }
    if (A[0] > B[0]) std::swap(A, B);
    for (int j = A[0]; j <= B[0]; j++)
    {
      // global pixel uv : t0[1] + i, j
      u = j;
      v = t0[1]+i;
      u0 = t0[0]; u1 = t1[0]; u2 = t2[0];
      v0 = t0[1]; v1 = t1[1]; v2 = t2[1];
      //printf("u0=%f, v0=%f\n", u0,v0);
      //printf("u1=%f, v1=%f\n", u1,v1);
      //printf("u2=%f, v2=%f\n", u2,v2);
      
      dd = (u2-u0)*(v1-v0)-(v2-v0)*(u1-u0);
      
      r0 = colorv0[0]; r1 = colorv1[0]; r2 = colorv2[0];  
      a0 = (r2-r0)*(v1-v0)-(v2-v0)*(r1-r0);
      a0 = a0 / dd;
      if (v1 != v0)
      {
        a1 = (r1-r0) - a0*(u1-u0);
        a1 = a1 / (v1-v0);
      }
      else if (v2 != v0)
      {
        a1 = (r2-r0) - a0*(u2-u0);
        a1 = a1 / (v2-v0);
      }
      else
      {
        a1 = (r2-r1) - a0*(u2-u1);
        a1 = a1 / (v2-v1);
      }
      a2 = r0 - a0*u0 - a1*v0;
      color[0] = int(a0*u + a1*v + a2);
      //printf("u=%f v=%f r=%d, r0=%d r1=%d r2=%d dd=%f \n ",u,v,color[0],colorv0[0],colorv1[0],colorv2[0],dd);

      r0 = colorv0[1]; r1 = colorv1[1]; r2 = colorv2[1];      
      a0 = (r2-r0)*(v1-v0)-(v2-v0)*(r1-r0);
      a0 = a0 / dd;

      if (v1 != v0)
      {
        a1 = (r1-r0) - a0*(u1-u0);
        a1 = a1 / (v1-v0);
      }
      else if (v2 != v0)
      {
        a1 = (r2-r0) - a0*(u2-u0);
        a1 = a1 / (v2-v0);
      }
      else
      {
        a1 = (r2-r1) - a0*(u2-u1);
        a1 = a1 / (v2-v1);
      }

      a2 = r0 - a0*u0 - a1*v0;
      color[1] = int(a0*u + a1*v + a2);

      r0 = colorv0[2]; r1 = colorv1[2]; r2 = colorv2[2];      
      a0 = (r2-r0)*(v1-v0)-(v2-v0)*(r1-r0);
      a0 = a0 / dd;
      
      if (v1 != v0)
      {
        a1 = (r1-r0) - a0*(u1-u0);
        a1 = a1 / (v1-v0);
      }
      else if (v2 != v0)
      {
        a1 = (r2-r0) - a0*(u2-u0);
        a1 = a1 / (v2-v0);
      }
      else
      {
        a1 = (r2-r1) - a0*(u2-u1);
        a1 = a1 / (v2-v1);
      }

      a2 = r0 - a0*u0 - a1*v0;
      color[2] = int(a0*u + a1*v + a2);

      // u = (j-A[0])/(t1[0]-t0[0])
      // v = i*1./(t1[1]-t0[1])
      // r = u*r1 + v*r2 + (1-u-v)*r0
      //color[0] = colorv0[0];
      //color[1] = colorv0[1];
      //color[2] = colorv0[2];
      SetPixel(dest, destWidth, j, t0[1] + i, color);
    }
  }
}

//  public-domain code by Darel Rex Finley, 2007
// http://alienryderflex.com/polygon_fill/
static void RasterizePolygon(uint8_t *dest, int destWidth, int vertices[][2], const int vertexCount, const uint8_t *color)
{
  int IMAGE_TOP = INT_MAX, IMAGE_BOT = 0, IMAGE_LEFT = INT_MAX, IMAGE_RIGHT = 0;
  for (int i = 0; i < vertexCount; i++) 
  {
    const int *vertex = vertices[i];
    IMAGE_TOP = vertex[1] < IMAGE_TOP ? vertex[1] : IMAGE_TOP;
    IMAGE_BOT = vertex[1] > IMAGE_BOT ? vertex[1] : IMAGE_BOT;
    IMAGE_LEFT = vertex[0] < IMAGE_LEFT ? vertex[0] : IMAGE_LEFT;
    IMAGE_RIGHT = vertex[0] > IMAGE_RIGHT ? vertex[0] : IMAGE_RIGHT;
  }
  int nodes, nodeX[255], pixelX, pixelY, i, j, swap;
  //  Loop through the rows of the image.
  for (pixelY=IMAGE_TOP; pixelY<IMAGE_BOT; pixelY++) 
  {
    //  Build a list of nodes.
    nodes=0; j=vertexCount-1;
    for (i=0; i<vertexCount; i++) 
    {
      if ((vertices[i][1]<(double)pixelY && vertices[j][1]>=(double)pixelY) || (vertices[j][1]<(double)pixelY && vertices[i][1]>=(double)pixelY))
      {
        nodeX[nodes++]=(int) (vertices[i][0]+(pixelY-vertices[i][1])/(vertices[j][1]-vertices[i][1])*(vertices[j][0]-vertices[i][0]));
      }
      j=i;
    }
    //  Sort the nodes, via a simple Bubble sort.
    i=0;
    while (i<nodes-1) 
    {
      if (nodeX[i]>nodeX[i+1]) 
      {
        swap=nodeX[i]; nodeX[i]=nodeX[i+1]; nodeX[i+1]=swap; if (i) i--; 
      }
      else 
      {
        i++;
      }
    }
    //  Fill the pixels between node pairs.
    for (i=0; i<nodes; i+=2) 
    {
      if (nodeX[i  ]>=IMAGE_RIGHT) break;
      if (nodeX[i+1]> IMAGE_LEFT ) 
      {
        if (nodeX[i  ]< IMAGE_LEFT) nodeX[i  ]=IMAGE_LEFT ;
        if (nodeX[i+1]> IMAGE_RIGHT) nodeX[i+1]=IMAGE_RIGHT;
        for (pixelX=nodeX[i]; pixelX<nodeX[i+1]; pixelX++)
          SetPixel(dest, destWidth, pixelX, pixelY, color);
      }
    }
  }
}

// ============================================================================
/* Get UV mapping of a surface (use xatlas) */
// ============================================================================
PyObject* K_GEOM::getUV(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float normalDeviationWeight;
  E_Float texResolution; 
  PyObject* fields; // champ a sortir en rgb
  if (!PYPARSETUPLE_(args, O_ RR_ O_, &array, &normalDeviationWeight, 
      &texResolution, &fields)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray2(array, varString, f, im, jm, km, cn, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getUV: invalid array.");
    return NULL;
  }

  E_Int api = f->getApi();

  if (res == 1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                      "getUV: ony for TRI or QUAD array.");
    return NULL;  
  }

  if (strcmp(eltType, "TRI") != 0 && strcmp(eltType, "QUAD") != 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                      "getUV: ony for TRI or QUAD array.");
    return NULL;  
  }

  // check coords
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getUV: can't find coordinates in array.");
    return NULL;  
  }
  posx++; posy++; posz++;

  // UV must have been already added
  E_Int posu = K_ARRAY::isNamePresent("_u_", varString);
  E_Int posv = K_ARRAY::isNamePresent("_v_", varString);
  if (posu == -1 || posv == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getUV: can't find _u_,_v_ in array.");
    return NULL;  
  }
  posu++; posv++;

  // Check fields (optional image rgb output)
  bool exportFields = false;
  std::vector<char*> fNames;
  if (fields != Py_None)
  {
    for (E_Int i = 0; i < PyList_Size(fields); i++) 
    {
      PyObject* tpl0 = PyList_GetItem(fields, i);
      if (PyString_Check(tpl0)) 
      {
        char* str = PyString_AsString(tpl0);
        fNames.push_back(str);
      }
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(tpl0)) 
      {
        char* str = (char*)PyUnicode_AsUTF8(tpl0);
        fNames.push_back(str);
      }
#endif
    }
  }
  E_Int posvx=0, posvy=0, posvz=0;
  if (fNames.size() == 0) exportFields = false;
  if (fNames.size() >= 1)
  {
    posvx = K_ARRAY::isNamePresent(fNames[0], varString);
    if (posvx >= 0) { exportFields = true; posvy = posvx; posvz = posvx; } 
  }
  if (fNames.size() >= 2 && exportFields)
  {
    posvy = K_ARRAY::isNamePresent(fNames[1], varString);
    if (posvy >= 0) { posvz = posvy; }
    else exportFields = false; 
  }
  if (fNames.size() >= 3 && exportFields)
  {
    posvz = K_ARRAY::isNamePresent(fNames[2], varString);
    if (posvz < 0) exportFields = false;
  }

  // Pass mesh to xatlas
  xatlas::Atlas *atlas = xatlas::Create();

  E_Int nvertex = f->getSize();
  E_Int nelts = cn->getSize();
  E_Int nf = cn->getNfld(); // 3 or 4

  // compact coords
  float* coords = new float [3*nvertex];
  E_Float* px = f->begin(posx);
  E_Float* py = f->begin(posy);
  E_Float* pz = f->begin(posz);

#define SCALEFACTOR 1000.

  for (E_Int i = 0; i < nvertex; i++) 
  {
    coords[3*i] = px[i] * SCALEFACTOR;
    coords[3*i+1] = py[i] * SCALEFACTOR;
    coords[3*i+2] = pz[i] * SCALEFACTOR;
  }  
  // ecrit connect
  uint32_t* connect = new uint32_t [nf*nelts];
  for (E_Int i = 0; i < nelts; i++)
    for (E_Int n = 0; n < nf; n++)
    {
      connect[nf*i+n] = (*cn)(i,n+1)-1;
    }

  //printf("nvert=%d %d %d\n", nvertex, nf, nelts);
  xatlas::MeshDecl meshDecl;
  meshDecl.vertexCount = (uint32_t)nvertex;
  meshDecl.vertexPositionData = coords;
  meshDecl.vertexPositionStride = sizeof(float) * 3;

  meshDecl.indexCount = (uint32_t)(nelts*nf);
  meshDecl.indexData = connect;

  uint8_t* fvc = NULL;
  if (nf == 4)
  {
    meshDecl.faceCount = (uint32_t)nelts;
    // size of each element
    fvc = new uint8_t [nelts];
    for (E_Int i = 0; i < nelts; i++) fvc[i] = 4;
    meshDecl.faceVertexCount = fvc;
  }
  meshDecl.indexFormat = xatlas::IndexFormat::UInt32;
  xatlas::AddMeshError error = xatlas::AddMesh(atlas, meshDecl, (uint32_t)1);
  if (error != xatlas::AddMeshError::Success)
  {
    xatlas::Destroy(atlas);
    delete [] coords;
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getUV: can not create atlas.");
    return NULL;
  }

  // Ajout chartOptions et packOptions
  printf("-> normal=%f texResolution=%f\n", normalDeviationWeight, texResolution);
  xatlas::ChartOptions opt1;
  opt1.normalDeviationWeight = normalDeviationWeight; // poids de la deviation des normales dans la metrique
  //opt1.maxCost = 2.; // si grand, moins de charts
  //opt1.maxIterations = 2; // improve quality
  //opt1.roundnessWeight = 0.01f;
  //opt1.straightnessWeight = 6.0f;
  //opt1.normalSeamWeight = 4.0f; // If > 1000, normal seams are fully respected.
  //opt1.textureSeamWeight = 0.5f;
  
  xatlas::PackOptions opt2;
  //opt2.texelsPerUnit = texelsPerUnit;
  opt2.resolution = texResolution; // 1920

  xatlas::Generate(atlas, opt1, opt2);
  printf("   %d charts\n", atlas->chartCount);
  printf("   %d atlases\n", atlas->atlasCount);
  for (uint32_t i = 0; i < atlas->atlasCount; i++)
    printf("      %d: %0.2f%% utilization\n", i, atlas->utilization[i] * 100.0f);
  printf("   %ux%u resolution\n", atlas->width, atlas->height);

  // output
  // Fill uv in orig model
  E_Float* pu = f->begin(posu);
  E_Float* pv = f->begin(posv);
  const xatlas::Mesh &mesh = atlas->meshes[0];
  E_Int nvertexOut = mesh.vertexCount;
  printf("output nvertex=%d\n", mesh.vertexCount);

  for (E_Int i = 0; i < nvertexOut; i++)
  {
    const xatlas::Vertex &vertex = mesh.vertexArray[i];
    E_Int iref = vertex.xref;
    //float *pos = &coords[iref * 3];
    //printf("v %g %g %g\n", pos[0], pos[1], pos[2]);
    pu[iref] = vertex.uv[0] / atlas->width;
    pv[iref] = vertex.uv[1] / atlas->height;
    //printf("vt %f %f\n", pu[i],pv[i]);
  }

  // cleanup
  delete [] coords;
  delete [] connect;
  delete [] fvc;


  //=========================
  // Export mesh with seams
  //=========================
  PyObject* o;
  if (exportFields)
    o = K_ARRAY::buildArray2(8, "x,y,z,u,v,vx,vy,vz", nvertexOut, nelts, -1, eltType, 0, 1, 1, api);
  else 
    o = K_ARRAY::buildArray2(5, "x,y,z,u,v", nvertexOut, nelts, -1, eltType, 0, 1, 1, api);

  FldArrayF* fo; FldArrayI* cno; 
  K_ARRAY::getFromArray2(o, fo, cno);
  E_Float* puo = fo->begin(4);
  E_Float* pvo = fo->begin(5);
  E_Float* pxo = fo->begin(1);
  E_Float* pyo = fo->begin(2);
  E_Float* pzo = fo->begin(3);
  E_Float* pvxo = NULL; E_Float* pvyo = NULL; E_Float* pvzo = NULL;
  E_Float* pvx = NULL; E_Float* pvy = NULL; E_Float* pvz = NULL;

  if (exportFields)
  {
    pvxo = fo->begin(6); pvx = f->begin(posvx+1);
    pvyo = fo->begin(7); pvy = f->begin(posvy+1);
    pvzo = fo->begin(8); pvz = f->begin(posvz+1);
  }
  for (E_Int i = 0; i < nvertexOut; i++)
  {
    const xatlas::Vertex &vertex = mesh.vertexArray[i];
    E_Int iref = vertex.xref;
    pxo[i] = px[iref];
    pyo[i] = py[iref];
    pzo[i] = pz[iref];
    puo[i] = vertex.uv[0] / atlas->width;
    pvo[i] = vertex.uv[1] / atlas->height;
    if (exportFields)
    {
      pvxo[i] = pvx[iref];
      pvyo[i] = pvy[iref];
      pvzo[i] = pvz[iref];
    }
  }

  uint32_t currentIndex = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    for (E_Int n = 1; n <= nf; n++)
    {
      (*cno)(i,n) = mesh.indexArray[currentIndex] + 1;
      //printf("%d %d - %d %d\n", i, n, nvertexOut, mesh.indexArray[currentIndex]+1);
      //(*cno)(i,n) = 1;
      currentIndex++;
    }
  }

  PyObject* tpl = PyList_New(0);
  PyList_Append(tpl, o); RELEASESHAREDU(o, fo, cno);

  if (atlas->width <= 0 || atlas->height <= 0) return tpl; 

  // si la vitesse est presente, on calcule le max de la norme
  E_Float vmax = 0.;
  E_Float dx, dy, dz;
  if (exportFields)
  { 
    for (E_Int i = 0; i < nvertex; i++)
    {
      dx = pvx[i]; dy = pvy[i]; dz = pvz[i];
      vmax = std::max(vmax, dx*dx+dy*dy+dz*dz);
    }
    vmax = sqrt(vmax);
    vmax = std::max(vmax, 1.e-12);
  }

  //===============================
  // export Texture Image as array
  //===============================
  printf("Rasterizing result...\n");
  std::vector<uint8_t> outputTrisImage, outputChartsImage, outputBumpImage, outputFieldImage;
  const uint32_t imageDataSize = atlas->width * atlas->height * 3;
  outputTrisImage.resize(atlas->atlasCount * imageDataSize);
  outputChartsImage.resize(atlas->atlasCount * imageDataSize);
  outputBumpImage.resize(atlas->atlasCount * imageDataSize);
  //const xatlas::Mesh &mesh = atlas->meshes[0];

  // Rasterize mesh triangles/polygons in Image
  const uint8_t white[] = { 255, 255, 255 };
  auto faceCount = (const uint32_t)nelts;
  uint32_t faceFirstIndex = 0;
  int verts[255][2];
  uint8_t color[3];
  uint8_t color0[3]; uint8_t color1[3]; uint8_t color2[3];
  
  for (uint32_t f = 0; f < faceCount; f++)
  {
    int32_t atlasIndex = -1;
    const uint32_t faceVertexCount = nf;
    for (uint32_t j = 0; j < faceVertexCount; j++) 
    {
      const xatlas::Vertex &v = mesh.vertexArray[mesh.indexArray[faceFirstIndex + j]];
      atlasIndex = v.atlasIndex; // The same for every vertex in the face.
      verts[j][0] = int(v.uv[0]);
      verts[j][1] = int(v.uv[1]);
    }
    if (atlasIndex < 0) continue; // Skip faces that weren't atlased.
    RandomColor(color);
    uint8_t *imageData = &outputTrisImage[atlasIndex * imageDataSize];
    if (exportFields)
    {
      // attention, only triangle for now
      const xatlas::Vertex &v0 = mesh.vertexArray[mesh.indexArray[faceFirstIndex + 0]];
      const xatlas::Vertex &v1 = mesh.vertexArray[mesh.indexArray[faceFirstIndex + 1]];
      const xatlas::Vertex &v2 = mesh.vertexArray[mesh.indexArray[faceFirstIndex + 2]];
      E_Int i0 = v0.xref;
      E_Int i1 = v1.xref;
      E_Int i2 = v2.xref;

      color0[0] = int(255*pvx[i0]/vmax);
      color0[1] = int(255*pvy[i0]/vmax);
      color0[2] = int(255*pvz[i0]/vmax);
      color1[0] = int(255*pvx[i1]/vmax); 
      color1[1] = int(255*pvy[i1]/vmax);
      color1[2] = int(255*pvz[i1]/vmax);
      color2[0] = int(255*pvx[i2]/vmax); 
      color2[1] = int(255*pvy[i2]/vmax);
      color2[2] = int(255*pvz[i2]/vmax);

      if (faceVertexCount == 3)
        RasterizeTriangle2(imageData, atlas->width, verts[0], verts[1], verts[2], color0, color1, color2);
    }
    else
    {
      if (faceVertexCount == 3)
        RasterizeTriangle(imageData, atlas->width, verts[0], verts[1], verts[2], color);
      else
        RasterizePolygon(imageData, atlas->width, verts, (int)faceVertexCount, color);
      for (uint32_t j = 0; j < faceVertexCount; j++)
        RasterizeLine(imageData, atlas->width, verts[j], verts[(j + 1) % faceVertexCount], white);
    }
    faceFirstIndex += faceVertexCount;
  }

  //============================================
  // Rasterize mesh triangles in Bumpmap image
  //============================================
  faceFirstIndex = 0;
  for (uint32_t f = 0; f < faceCount; f++)
  {
    int32_t atlasIndex = -1;
    int verts[255][2];
    const uint32_t faceVertexCount = nf;
    for (uint32_t j = 0; j < faceVertexCount; j++) 
    {
      const xatlas::Vertex &v = mesh.vertexArray[mesh.indexArray[faceFirstIndex + j]];
      atlasIndex = v.atlasIndex; // The same for every vertex in the face.
      verts[j][0] = int(v.uv[0]);
      verts[j][1] = int(v.uv[1]);
    }
    if (atlasIndex < 0) continue; // Skip faces that weren't atlased.
    uint8_t color[3];
    color[0] = 126; color[1] = 126; color[2] = 255;
    uint8_t color2[3];
    color2[0] = 126; color2[1] = 126; color2[2] = 250;
    
    uint8_t *imageData = &outputBumpImage[atlasIndex * imageDataSize];
    if (faceVertexCount == 3)
      RasterizeTriangle(imageData, atlas->width, verts[0], verts[1], verts[2], color);
    else
      RasterizePolygon(imageData, atlas->width, verts, (int)faceVertexCount, color);
    for (uint32_t j = 0; j < faceVertexCount; j++)
      RasterizeLine(imageData, atlas->width, verts[j], verts[(j + 1) % faceVertexCount], color2);
    faceFirstIndex += faceVertexCount;
  }

  // Rasterize mesh triangles in field image
  /*
  if (posvx >= 0 && posvy >= 0 && posvz >= 0)
  {
    E_Float* pvx = f->begin(posvx+1);
    E_Float* pvy = f->begin(posvy+1);
    E_Float* pvz = f->begin(posvz+1);
    
    outputFieldImage.resize(atlas->atlasCount * imageDataSize);
    faceFirstIndex = 0;
    for (uint32_t f = 0; f < faceCount; f++)
    {
      int32_t atlasIndex = -1;
      int verts[255][2];
      const uint32_t faceVertexCount = nf;
      for (uint32_t j = 0; j < faceVertexCount; j++) 
      {
        const xatlas::Vertex &v = mesh.vertexArray[mesh.indexArray[faceFirstIndex + j]];
        atlasIndex = v.atlasIndex; // The same for every vertex in the face.
        verts[j][0] = int(v.uv[0]);
        verts[j][1] = int(v.uv[1]);
      }
      if (atlasIndex < 0) continue; // Skip faces that weren't atlased.
      uint8_t color[3];
      color[0] = pvx[];
      
      uint8_t *imageData = &outputFieldImage[atlasIndex * imageDataSize];
      if (faceVertexCount == 3)
        RasterizeTriangle(imageData, atlas->width, verts[0], verts[1], verts[2], color);
      else
        RasterizePolygon(imageData, atlas->width, verts, (int)faceVertexCount, color);
      for (uint32_t j = 0; j < faceVertexCount; j++)
        RasterizeLine(imageData, atlas->width, verts[j], verts[(j + 1) % faceVertexCount], color2);
      faceFirstIndex += faceVertexCount;
    }
  }
  */

  // Rasterize mesh charts.
  /*
  for (uint32_t j = 0; j < mesh.chartCount; j++) 
  {
    const xatlas::Chart *chart = &mesh.chartArray[j];
    uint8_t color[3];
    RandomColor(color);
    for (uint32_t k = 0; k < chart->faceCount; k++)
    {
      const uint32_t face = chart->faceArray[k];
      const uint32_t faceVertexCount = nf;
      faceFirstIndex = 0;
      for (uint32_t l = 0; l < face; l++)
      faceFirstIndex += nvertex;
      int verts[255][2];
      for (uint32_t l = 0; l < faceVertexCount; l++) 
      {
        const xatlas::Vertex &v = mesh.vertexArray[mesh.indexArray[faceFirstIndex + l]];
        verts[l][0] = int(v.uv[0]);
        verts[l][1] = int(v.uv[1]);
      }
      uint8_t *imageData = &outputChartsImage[chart->atlasIndex * imageDataSize];
      if (faceVertexCount == 3)
        RasterizeTriangle(imageData, atlas->width, verts[0], verts[1], verts[2], color);
      else
        RasterizePolygon(imageData, atlas->width, verts, (int)faceVertexCount, color);
      for (uint32_t l = 0; l < faceVertexCount; l++)
        RasterizeLine(imageData, atlas->width, verts[l], verts[(l + 1) % faceVertexCount], white);
    }
  }
  */

  // Export images as a cartesian zone
  for (uint32_t i = 0; i < atlas->atlasCount; i++) 
  {
    E_Int width = atlas->width;
    E_Int height = atlas->height;
    uint8_t *imageData = &outputTrisImage[i * imageDataSize];
    PyObject* o = K_ARRAY::buildArray2(6, "x,y,z,r,g,b", width,height,1, api);
    FldArrayF* fo; K_ARRAY::getFromArray2(o, fo);
    E_Float* pr = fo->begin(4);
    E_Float* pg = fo->begin(5);
    E_Float* pb = fo->begin(6);
    E_Float* px = fo->begin(1);
    E_Float* py = fo->begin(2);
    E_Float* pz = fo->begin(3);
      
    for (E_Int lj = 0; lj < height; lj++)
    for (E_Int li = 0; li < width;  li++)
    {
      E_Int ind = li+lj*width;
      px[ind] = li; py[ind] = lj; pz[ind] = 0.;
      pr[ind] = imageData[3*ind];
      pg[ind] = imageData[3*ind+1];
      pb[ind] = imageData[3*ind+2];
    }
    PyList_Append(tpl, o); RELEASESHAREDS(o, fo);
  }
  
  // Export images as a cartesian zone
  for (uint32_t i = 0; i < atlas->atlasCount; i++) 
  {
    E_Int width = atlas->width;
    E_Int height = atlas->height;
    uint8_t *imageData = &outputBumpImage[i * imageDataSize];
    PyObject* o = K_ARRAY::buildArray2(6, "x,y,z,r,g,b", width,height,1, api);
    FldArrayF* fo; K_ARRAY::getFromArray2(o, fo);
    E_Float* pr = fo->begin(4);
    E_Float* pg = fo->begin(5);
    E_Float* pb = fo->begin(6);
    E_Float* px = fo->begin(1);
    E_Float* py = fo->begin(2);
    E_Float* pz = fo->begin(3);
      
    for (E_Int lj = 0; lj < height; lj++)
    for (E_Int li = 0; li < width;  li++)
    {
      E_Int ind = li+lj*width;
      px[ind] = li; py[ind] = lj; pz[ind] = 0.;
      pr[ind] = imageData[3*ind];
      pg[ind] = imageData[3*ind+1];
      pb[ind] = imageData[3*ind+2];
    }
    PyList_Append(tpl, o); RELEASESHAREDS(o, fo);
  }

  // End
  xatlas::Destroy(atlas);
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
