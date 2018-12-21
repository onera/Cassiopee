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
# include "generator.h"
using namespace std;
using namespace K_FLD;
# include "Tetgen/tetgen.h"
int __Error__ = 0;
jmp_buf __env__;

//=========================================================================
/* Generation de maillage tetra a partir d'un maillage surfacique (tetgen) 
   Le maillage surfacique n'est pas modifie */
//=========================================================================
PyObject* K_GENERATOR::tetgen(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float grading, maxh;
  PyObject* holes; // coord. d'un point a l'interieur
#if defined E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "OddO", &array, &maxh, &grading, &holes)) 
    return NULL;
#else
  if (!PyArg_ParseTuple(args, "OffO", &array, &maxh, &grading, &holes)) 
    return NULL;
#endif
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res <= 0) return NULL;
  if (res == 1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "tetgen: input must be TRI.");
    return NULL;
  }

  if (K_STRING::cmp(eltType, "TRI") != 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "tetgen: input must be TRI.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "tetgen: coordinates not found in array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  
  // Getting hole list
  E_Int nbholes = 0;
  E_Float* pholes = NULL;
  if (PyList_Check(holes) == true)
  {
    E_Int n = PyList_Size(holes);
    for (E_Int i = 0; i < n; i++)
    {
      PyObject* o = PyList_GetItem(holes, i);
      if (PyList_Check(o) == true && PyList_Size(o) == 3)
        nbholes++;
    }
    pholes = new E_Float [nbholes*3];
    for (E_Int i = 0; i < n; i++)
    {
      PyObject* o = PyList_GetItem(holes, i);
      for (E_Int j = 0; j < 3; j++)
        pholes[3*i+j] = PyFloat_AsDouble(PyList_GetItem(o, j));
    }
  }

  PyObject* pout = NULL;

  tetgenio in, out;
  tetgenbehavior b;
  // Maille l'interieur d'un PLC
  b.plc = 1;
  // Ne touche pas le maillage surfacique
  b.nobisect = 1; // -Y (conserve la frontiere ext)
  // b.nobisect_param = 1; // 0, 1 ou 2
  // b.addsteiner_algo = 1; // 1 ou 2
  // Ajoute des pts pour la qualite
  b.quality = 1; // -q (ajoute de points)
  //b.minratio = 1.2; // radius-edge ratio >= 0.612
  //b.diagnose = 1;

  //b.regionattrib = 1; // met des attributs de region

  if (maxh > 0)
  {
    b.fixedvolume = 1;
    b.maxvolume = (maxh*maxh*maxh)*0.25; // -a (max volume constraint)
  }
  else b.varvolume = 1;

  // coplanar test
  b.epsilon = 1.e-10;

  // forced simple
  //b.nobisect = 0;
  //b.quality = 0;
  //b.fixedvolume = 0;
  //b.epsilon = 1.e-8;

  // pour du remaillage coarsening
  // b.refine = 1;
  // b.coarsen = 1;
  // b.coarsen_param = 1; // ??
  // b.coarsen_percent = 0.5; // x % less points

  /*
  char file[256];
  strcpy(file, "example");
  in.load_plc(file, tetgenbehavior::POLY);
  printf("loaded.\n");
  */

  // Remplissage de in a partir de array
  in.mesh_dim = 3;
  in.numberofpointattributes = 0;

  E_Int n = f->getSize();
  E_Float* ptx = f->begin(posx);
  E_Float* pty = f->begin(posy);
  E_Float* ptz = f->begin(posz);
  in.numberofpoints = n;
  in.pointlist = new E_Float [n*3];
  E_Float* ptin = in.pointlist;
  for (E_Int i = 0; i < n; i++)
  {
    ptin[0] = ptx[i]; ptin[1] = pty[i]; ptin[2] = ptz[i]; ptin += 3;
  }
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);
  E_Int ne = cn->getSize();
  in.numberoffacets = ne;
  in.facetlist = new tetgenio::facet[ne];
  for (E_Int i = 0; i < ne; i++)
  {
    tetgenio::facet* f = &(in.facetlist[i]);
    tetgenio::init(f);
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[1];
    tetgenio::polygon* p = &f->polygonlist[0];
    tetgenio::init(p);
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    
    p->vertexlist[0] = cn1[i]-1;
    p->vertexlist[1] = cn2[i]-1;
    p->vertexlist[2] = cn3[i]-1;
  }

  // Ajoute des holes
  // in.numberofholes = 1;
  // in.holelist = new double [3* in.numberofholes];
  // in.holelist[0] = 0.5; // x
  // in.holelist[1] = 0.5; // y
  // in.holelist[2] = 0.5; // z

  in.numberofholes = nbholes;
  in.holelist = new double [3* in.numberofholes];
  for (E_Int i = 0; i < nbholes; i++)
  {
    in.holelist[3*i+0] = pholes[3*i+0];
    in.holelist[3*i+1] = pholes[3*i+1];
    in.holelist[3*i+2] = pholes[3*i+2];
  }
  delete [] pholes;
  

  // Ajoute les regions lists
  //in.numberofregions = 1;
  //in.regionlist = new double [5* in.numberofregions];
  //in.regionlist[0] = 0.5; // x
  //in.regionlist[1] = 0.5; // y 
  //in.regionlist[2] = 0.5; // z
  //in.regionlist[3] = 1; // entier reperant la region (region attrib)
  //in.regionlist[4] = 0.01; // ?

  RELEASESHAREDU(array, f, cn);
  
  // in: surface mesh or previous tetra mesh
  // out: output mesh
  // addin: if != NULL, list of constraint points
  // bgmin: if != NULL, background mesh with size function
  int val;
  val = setjmp(__env__); // enregistre l'environnement pour retour d'erreur
  if (val == 0) // premier passage
    tetrahedralize(&b, &in, &out, NULL, NULL);

  if (val == 1) // erreur
  { 
    PyErr_SetString(PyExc_TypeError,
                    "tetgen: failed.");
    return NULL; 
  }

  /* Build output */
  E_Int np = out.numberofpoints;
  ne = out.numberoftetrahedra;
  printf("Generate %d points and %d elements.\n",np,ne);
  pout = K_ARRAY::buildArray(3, "x,y,z", np, ne, 4, NULL, false, 0);
  E_Float* fp = K_ARRAY::getFieldPtr(pout);
  E_Float* fx = fp; E_Float* fy = fp+np; E_Float* fz = fp+2*np;
  E_Float* pt = out.pointlist;
  for (E_Int i = 0; i < np; i++)
  {
    fx[i] = pt[0]; fy[i] = pt[1]; fz[i] = pt[2]; pt += 3;
  }
  E_Int* cp = K_ARRAY::getConnectPtr(pout);
  E_Int* cp1 = cp; E_Int* cp2 = cp+ne; 
  E_Int* cp3 = cp+2*ne; E_Int* cp4 = cp+3*ne;
  E_Int ind = 0;
  int* tl = out.tetrahedronlist;
  for (E_Int i = 0; i < ne; i++)
  { 
    cp1[i] = tl[ind]+1;
    cp2[i] = tl[ind+1]+1;
    cp3[i] = tl[ind+2]+1;
    cp4[i] = tl[ind+3]+1;
    ind += 4;
  }

  return pout;
}
