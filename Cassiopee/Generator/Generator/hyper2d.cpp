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

// Hyperbolic mesh generators

# include "generator.h"
# include "CompGeom/compGeom.h"
# include <vector>

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Fortran functions */
// ============================================================================
extern "C"
{
 
  void k6hyper2d_(const E_Int& ni, const E_Int& nj,
                  const E_Float* distrib,
                  const E_Float* xi, const E_Float* yi, const E_Float* zi,
                  const E_Int& type,
                  E_Float* xo, E_Float* yo, E_Float* zo,
                  E_Int* IP, E_Float* A, E_Float* B, E_Float* C,
                  E_Float* RHS, E_Float* Z, E_Float* ZA, E_Float* vol,
                  const E_Int& eta_start, const E_Int& eta_end, const E_Float& beta);

  void k6hyper2d2_(const E_Int& ni, const E_Int& nj,
                   const E_Float* distrib,
                   const E_Float* xi, const E_Float* yi, const E_Float* zi,
                   const E_Int& type, const E_Float& alpha,
                   E_Float* xo, E_Float* yo, E_Float* zo,
                   E_Int* IP, E_Float* A, E_Float* B, E_Float* C,
                   E_Float* RHS, E_Float* Z, E_Float* ZA, E_Float* vol);

  void k6hyper2d3_(const E_Int& ni, const E_Int& nj,
                   const E_Float* distrib,
                   const E_Float* xi, const E_Float* yi, const E_Float* zi,
                   const E_Int& type, 
                   const E_Float& alpha1, const E_Float& alpha2,
                   E_Float* xo, E_Float* yo, E_Float* zo,
                   E_Int* IP, E_Float* A, E_Float* B, E_Float* C,
                   E_Float* RHS, E_Float* Z, E_Float* ZA, 
                   E_Float* vol, E_Float* alphad);

  void k6hyper2d4_(const E_Int& ni, const E_Int& nj,
                   const E_Float* distrib,
                   const E_Float* xi, const E_Float* yi, const E_Float* zi,
                   const E_Int& type,
                   E_Float* xo, E_Float* yo, E_Float* zo,
                   E_Int* IP, E_Float* A, E_Float* B, E_Float* C,
                   E_Float* RHS, E_Float* Z, E_Float* ZA, E_Float* vol);
}

// ============================================================================
/* Hyperbolic mesh generation for a mesh. Orthogonal. */
// ============================================================================
PyObject* K_GENERATOR::hyper2DMesh(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayd;
  char* meshType;
  E_Int eta_start, eta_end;
  E_Float beta;

  if (!PYPARSETUPLE_(args, OO_ S_ II_ R_, 
                    &array, &arrayd, &meshType, 
                    &eta_start, &eta_end, &beta))
    return NULL;

  E_Int type = 0;
  if (strcmp(meshType, "C") == 0) type = 0;
  else if (strcmp(meshType, "O") == 0) type = 1;
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "hyper2D: unknown mesh type.");
    return NULL;
  }

  // Check array
  E_Int ni, nj, nk, nid, njd, nkd;
  FldArrayF* f;
  FldArrayF* fd;
  char* varString; char* eltType;
  FldArrayI* cn;
  char* varStringd;
  char* eltTyped;
  FldArrayI* cnd;
  FldArrayF coord1, s, dx, dy, dz;
  FldArrayF coord;
  FldArrayI IP;
  FldArrayF A, B, C, RHS, Z, ZA, vol;
  vector<E_Int> pos;
  vector<E_Int> posd;
  E_Int posx, posy, posz, posxd, posyd, poszd;

  // Extraction des infos sur le maillage
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  // Extraction des infos sur la distribution
  E_Int resd = 
    K_ARRAY::getFromArray3(arrayd, varStringd, fd, nid, njd, nkd, cnd, eltTyped);

  E_Int api = f->getApi();
 
  if (res == 1 && resd == 1)
  {
    char* varString0 = new char[strlen(varString)+strlen(varStringd)+4];
    E_Int res0 = 
      K_ARRAY::getPosition(varString, varStringd, pos, posd, varString0);
    if (res0 == -1)
    {
      RELEASESHAREDS(array, f);
      RELEASESHAREDS(arrayd, fd);
      PyErr_SetString(PyExc_TypeError,
                      "hyper2d: common variables list is empty.");
      return NULL;
    } 
    else if (res0 == 0)
    {
      printf("Warning: hyper2d: some variables are different.");
      printf(" Only common variables are kept.\n");
    }
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    posxd = K_ARRAY::isCoordinateXPresent(varStringd);
    posyd = K_ARRAY::isCoordinateYPresent(varStringd);
    poszd = K_ARRAY::isCoordinateZPresent(varStringd);

    if (posx == -1 || posy == -1 || posz == -1 ||
        posxd == -1 || posyd == -1 || poszd == -1)
    {
      RELEASESHAREDS(array, f);
      RELEASESHAREDS(arrayd, fd);
      PyErr_SetString(PyExc_TypeError,
                      "hyper2d: coordinates not found.");
      return NULL;
    }
    posx++; posy++; posz++;
    posxd++; posyd++; poszd++;

    // check forcing : si la distribution est a peu 
    // pres celle de la grille fournie, on ne remaille pas le profil
    E_Bool forced = false;
    if (nid == ni)
    {
        // abscisse curv
        s.malloc(ni);
        E_Float* fx = f->begin(posx);
        E_Float* fy = f->begin(posy);
        E_Float* fz = f->begin(posz);
        E_Float* fdx = fd->begin(posx);
        E_Float L = 0.;
        E_Float l, lx, ly, lz;
        s[0] = 0.;
        for (E_Int i = 0; i < ni-1; i++)
        {
            lx = fx[i+1]-fx[i];
            ly = fy[i+1]-fy[i];
            lz = fz[i+1]-fz[i];
            l = sqrt(lx*lx+ly*ly+lz*lz);
            s[i+1] = s[i]+l;
            L += l; 
        }
        for (E_Int i = 0; i < ni; i++)
        {
            s[i] = s[i]/L;
        }
        forced = true;
        for (E_Int i = 0; i < ni; i++)
        {
            if (std::abs(s[i]-fdx[i])>1.e-6) forced = false;
        }
    }

    // Make the one D mapping
    coord1.malloc(nid, 3);
    if (forced == false)
    {
        s.malloc(ni);
        dx.malloc(ni); dy.malloc(ni); dz.malloc(ni);
	K_COMPGEOM::onedmap(ni, f->begin(posx), f->begin(posy), f->begin(posz),
			    nid, fd->begin(posxd),
			    coord1.begin(1), coord1.begin(2), coord1.begin(3),
			    s.begin(), dx.begin(), dy.begin(), dz.begin());
    }
    else
    {
        printf("forcing...\n");
        E_Float* fx = f->begin(posx);
        E_Float* fy = f->begin(posy);
        E_Float* fz = f->begin(posz);
        E_Float* c1x = coord1.begin(1);
        E_Float* c1y = coord1.begin(2);
        E_Float* c1z = coord1.begin(3);
        for (E_Int i = 0; i < nid; i++)
        {
            c1x[i] = fx[i]; c1y[i] = fy[i]; c1z[i] = fz[i];
        }
    }

    // Generate the mesh using hyperbolic grid generator
    coord.malloc(nid*njd, 3);

    if (type == 0)
    {
        IP.malloc(2*(nid-2));
        A.malloc(2*2*(nid-2));
        B.malloc(2*2*(nid-2));
        C.malloc(2*2*(nid-2));
        RHS.malloc(2*(nid-2));
    }
    else
    {
        IP.malloc(2*(nid-2));
        A.malloc(2*2*(nid-1));
        B.malloc(2*2*(nid-1));
        C.malloc(2*2*(nid-1));
        RHS.malloc(2*(nid-1));
        Z.malloc(2*2*(nid-1));
        ZA.malloc(2*(nid-2));
    }
    vol.malloc(nid*njd);
    
    k6hyper2d_(nid, njd,
               fd->begin(posyd),
               coord1.begin(1), coord1.begin(2), coord1.begin(3),
               type,
               coord.begin(1), coord.begin(2), coord.begin(3),
               IP.begin(), A.begin(), B.begin(), C.begin(),
               RHS.begin(), Z.begin(), ZA.begin(), vol.begin(),
               eta_start, eta_end, beta);
    
    // Array Creation 
    RELEASESHAREDS(array, f);
    RELEASESHAREDS(arrayd, fd);
    PyObject* tpl = 
      K_ARRAY::buildArray3(coord, varString0, nid, njd, 1, api);
    delete [] varString0;
    return tpl;
  }
  else if (res == 2 || resd == 2)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resd, arrayd, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "hyper2d: not used for unstructured arrays.");
    return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "hyper2d: unknown type of array.");
    return NULL;
  }
}

// ============================================================================
/* Hyperbolic mesh generation for a mesh. Constant alpha angle. */
// ============================================================================
PyObject* K_GENERATOR::hyper2D2Mesh(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayd;
  char* meshType;
  E_Float alpha;
  if ( !PYPARSETUPLE_(args, OO_ S_ R_, &array, &arrayd,
                      &meshType, &alpha ) )
  {
    return NULL;
  }

  //check arguments
  E_Int type = 0;
  if (strcmp(meshType, "C") == 0)
    type = 0;
  else if (strcmp(meshType, "O") == 0)
    type = 1;
  else
  {
     PyErr_SetString(PyExc_TypeError,
                     "hyper2D2: unknown mesh type.");
    return NULL;
  }

  // Check array
  E_Int ni, nj, nk, nid, njd, nkd;
  FldArrayF* f; FldArrayF* fd;
  char* eltType;
  FldArrayI* cn;
  char* varStringd; char* varString;
  char* eltTyped;
  FldArrayI* cnd;
  FldArrayF coord1, s, dx, dy, dz;
  FldArrayF coord;
  FldArrayI IP;
  FldArrayF A, B, C, RHS, Z, ZA, vol;
  vector<E_Int> pos;
  vector<E_Int> posd;
  E_Int posx, posy, posz, posxd, posyd, poszd;
  
  //extraction des infos sur le maillage
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  //extraction des infos sur la distribution
  E_Int resd = 
    K_ARRAY::getFromArray3(arrayd, varStringd, fd, 
                           nid, njd, nkd, cnd, eltTyped);

  E_Int api = f->getApi();
 
  if (res == 1 && resd == 1)
  {
    char* varString0 = new char [strlen(varString)+strlen(varStringd)+4];
    E_Int res0 = 
      K_ARRAY::getPosition(varString, varStringd, pos, posd, varString0);
    if (res0 == -1)
    {
      RELEASESHAREDS(array, f);
      RELEASESHAREDS(arrayd, fd);
      PyErr_SetString(PyExc_TypeError,
                      "hyper2d2: common variables list is empty.");
      return NULL;
    } 
    else if (res0 == 0)
    {
      printf(" Warning: hyper2d2: some variables are different.");
      printf(" Only common variables are kept.\n");
    }
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    posxd = K_ARRAY::isCoordinateXPresent(varStringd);
    posyd = K_ARRAY::isCoordinateYPresent(varStringd);
    poszd = K_ARRAY::isCoordinateZPresent(varStringd);

    if (posx == -1 || posy == -1 || posz == -1 ||
        posxd == -1 || posyd == -1 || poszd == -1)
    {
      RELEASESHAREDS(array, f);
      RELEASESHAREDS(arrayd, fd);
      PyErr_SetString(PyExc_TypeError,
                      "hyper2d2: coordinates not found.");
      return NULL;
    }
    posx++; posy++; posz++;
    posxd++; posyd++; poszd++;
      
    // Make the one D mapping
    coord1.malloc(nid, 3);
    s.malloc(ni);
    dx.malloc(ni); dy.malloc(ni); dz.malloc(ni);
    K_COMPGEOM::onedmap( ni,  f->begin(posx), f->begin(posy), f->begin(posz),
			 nid, fd->begin(posxd),
			 coord1.begin(1), coord1.begin(2), coord1.begin(3),
			 s.begin(), dx.begin(), dy.begin(), dz.begin());
    
    // Generate the mesh using hyperbolic grid generator
    coord.malloc((nid+1)*(njd+1), 3);
    IP.malloc(2*(nid-1));
    A.malloc(2*2*nid);
    B.malloc(2*2*nid);
    C.malloc(2*2*nid);
    RHS.malloc(2*nid);
    Z.malloc(2*2*(nid-1));
    ZA.malloc(2*(nid-2));
    vol.malloc(nid*njd);
    
    k6hyper2d2_(nid, njd, 
                fd->begin(posyd),
                coord1.begin(1), coord1.begin(2), coord1.begin(3),
                type, alpha,
                coord.begin(1), coord.begin(2), coord.begin(3),
                IP.begin(), A.begin(), B.begin(), C.begin(),
                RHS.begin(), Z.begin(), ZA.begin(), vol.begin());
    
    coord.reAllocMat(nid*njd, 3);
    // Array Creation 
    RELEASESHAREDS(array, f);
    RELEASESHAREDS(arrayd, fd);
    PyObject* tpl = 
      K_ARRAY::buildArray3(coord, varString0, nid, njd, 1, api);
    delete [] varString0;
    return tpl;
  }
  else if (res == 2 || resd == 2)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resd, arrayd, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "hyper2d2: not used for unstructured arrays.");
    return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "hyper2d2: unrecognised type of array.");
    return NULL;
  }
}

// ============================================================================
/* Hyperbolic mesh generation for a mesh. Given alpha boundary angles. */
// ============================================================================
PyObject* K_GENERATOR::hyper2D3Mesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* arrayd;
  char* meshType;
  E_Float alpha1, alpha2;

  if (!PYPARSETUPLE_(args, OO_ S_ II_, &array, &arrayd,
                     &meshType, &alpha1, &alpha2))
  {
    return NULL;
  }

  E_Int type = 0;
  if (strcmp(meshType, "C") == 0)
    type = 0;
  else if (strcmp(meshType, "O") == 0)
    type = 1;
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "hyper2D3: unknown mesh type.");
    return NULL;
  }

  // Check array
  E_Int ni, nj, nk, nid, njd, nkd;
  FldArrayF* f;
  FldArrayF* fd;
  char* eltType;
  FldArrayI* cn;
  char* varString;
  char* varStringd;
  char* eltTyped;
  FldArrayI* cnd;
  FldArrayF coord1, s, dx, dy, dz;
  FldArrayF coord;
  FldArrayI IP;
  FldArrayF A, B, C, RHS, Z, ZA, vol, alphad;
  vector<E_Int> pos;
  vector<E_Int> posd;
  E_Int posx, posy, posz, posxd, posyd, poszd;
    
  //extraction des infos sur le maillage
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  //extraction des infos sur la distribution
  E_Int resd = 
    K_ARRAY::getFromArray3(arrayd, varStringd, fd, 
                           nid, njd, nkd, cnd, eltTyped);

  E_Int api = f->getApi();
 
  if ( res == 1 && resd == 1)
  {
    char* varString0 = new char [strlen(varString)+strlen(varStringd)+4];
    E_Int res0 = 
      K_ARRAY::getPosition(varString, varStringd, pos, posd, varString0);
    if (res0 == -1)
    {
      RELEASESHAREDS(array, f);
      RELEASESHAREDS(arrayd, fd);
      PyErr_SetString(PyExc_TypeError,
                      "hyper2d3: common variables list is empty.");
      return NULL;
    } 
    else if (res0 == 0)
    {
      printf(" Warning: hyper2d3: some variables are different.");
      printf(" Only common variables are kept.\n");
    }
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    posxd = K_ARRAY::isCoordinateXPresent(varStringd);
    posyd = K_ARRAY::isCoordinateYPresent(varStringd);
    poszd = K_ARRAY::isCoordinateZPresent(varStringd);

    if (posx == -1 || posy == -1 || posz == -1 ||
        posxd == -1 || posyd == -1 || poszd == -1)
    {
      RELEASESHAREDS(array, f);
      RELEASESHAREDS(arrayd, fd);
      PyErr_SetString(PyExc_TypeError,
                      "hyper2d3: coordinates not found.");
      return NULL;
    }
    posx++; posy++; posz++;
    posxd++; posyd++; poszd++;
      
    // Make the one D mapping
    coord1.malloc(nid, 3);
    s.malloc(ni);
    dx.malloc(ni);
    dy.malloc(ni);
    dz.malloc(ni);
    
    K_COMPGEOM::onedmap( ni, f->begin(posx), f->begin(posy), f->begin(posz),
                nid, fd->begin(posxd),
                coord1.begin(1), coord1.begin(2), coord1.begin(3),
                s.begin(), dx.begin(), dy.begin(), dz.begin());
    
    // Generate the mesh using hyperbolic grid generator
    coord.malloc((nid+1)*(njd+1), 3);
    IP.malloc(2*(nid-1));
    A.malloc(2*2*nid);
    B.malloc(2*2*nid);
    C.malloc(2*2*nid);
    RHS.malloc(2*nid);
    Z.malloc(2*2*(nid-1));
    ZA.malloc(2*(nid-2));
    vol.malloc(nid*njd);
    alphad.malloc(nid*njd);
    
    k6hyper2d3_(nid, njd, fd->begin(posyd),
                coord1.begin(1), coord1.begin(2), coord1.begin(3),
                type, alpha1, alpha2,
                coord.begin(1), coord.begin(2), coord.begin(3),
                IP.begin(), A.begin(), B.begin(), C.begin(),
                RHS.begin(), Z.begin(), ZA.begin(), vol.begin(),
                alphad.begin());

    coord.reAllocMat(nid*njd, 3);
    // Array Creation 
    RELEASESHAREDS(array, f);
    RELEASESHAREDS(arrayd, fd);
    PyObject* tpl = 
      K_ARRAY::buildArray3(coord, varString0, nid, njd, 1, api);
    delete [] varString0;
    return tpl;
  }
  else if (res == 2 || resd == 2)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resd, arrayd, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "hyper2d3: not used for unstructured arrays.");
    return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "hyper2d3: unrecognised type of array.");
    return NULL;
  }
}
// ============================================================================
/* Hyperbolic mesh generation for a mesh. Orthogonal. 
   New version. */
// ============================================================================
PyObject* K_GENERATOR::hyper2D4Mesh(PyObject* self,
                                    PyObject* args)
{
  PyObject* array;
  PyObject* arrayd;
  char* meshType;
  if (!PYPARSETUPLE_(args, OO_ S_, &array, &arrayd, &meshType))
  {
    return NULL;
  }
  
  E_Int type = 0;
  if (strcmp(meshType, "C") == 0)
    type = 0;
  else if (strcmp(meshType, "O") == 0)
    type = 1;
  else
  {
    PyErr_SetString(PyExc_TypeError,
                     "hyper2D4: unknown mesh type.");
    return NULL;
  }
 
  // Check array
  E_Int ni, nj, nk, nid, njd, nkd;
  FldArrayF* f;
  FldArrayF* fd;
  char* varString; char* eltType;
  FldArrayI* cn;
  char* varStringd; char* eltTyped;
  FldArrayI* cnd;
  FldArrayF coord1, s, dx, dy, dz;
  FldArrayF coord;
  FldArrayI IP;
  FldArrayF A, B, C, RHS, Z, ZA, vol;
  vector<E_Int> pos;
  vector<E_Int> posd;
  E_Int posx, posy, posz, posxd, posyd, poszd;

  //extraction des infos sur le maillage
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  //extraction des infos sur la distribution
  E_Int resd = 
    K_ARRAY::getFromArray3(arrayd, varStringd, fd, 
                           nid, njd, nkd, cnd, eltTyped);

  E_Int api = f->getApi();
 
  if ( res == 1 && resd == 1)
  {
    char* varString0 = new char [strlen(varString)+strlen(varStringd)+4];
    E_Int res0 = 
      K_ARRAY::getPosition(varString, varStringd, pos, posd, varString0);
    if (res0 == -1)
    {
      RELEASESHAREDS(array, f);
      RELEASESHAREDS(arrayd, fd);
      PyErr_SetString(PyExc_TypeError,
                      "hyper2d4: common variables list is empty.");
      return NULL;
    } 
    else if (res0 == 0)
    {
      printf(" Warning: hyper2d4: some variables are different.");
      printf(" Only common variables are kept.\n");
    }
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    posxd = K_ARRAY::isCoordinateXPresent(varStringd);
    posyd = K_ARRAY::isCoordinateYPresent(varStringd);
    poszd = K_ARRAY::isCoordinateZPresent(varStringd);

    if (posx == -1 || posy == -1 || posz == -1 ||
        posxd == -1 || posyd == -1 || poszd == -1)
    {
      RELEASESHAREDS(array, f);
      RELEASESHAREDS(arrayd, fd);
      PyErr_SetString(PyExc_TypeError,
                      "hyper2d4: coordinates not found.");
      return NULL;
    }
    posx++; posy++; posz++;
    posxd++; posyd++; poszd++; 
   
    // Make the one D mapping
    coord1.malloc(nid, 3);
    s.malloc(ni);
    dx.malloc(ni);
    dy.malloc(ni);
    dz.malloc(ni);

    K_COMPGEOM::onedmap(ni, f->begin(posx), f->begin(posy), f->begin(posz),
			nid, fd->begin(posxd),
			coord1.begin(1), coord1.begin(2), coord1.begin(3),
			s.begin(), dx.begin(), dy.begin(), dz.begin());
  
    // Generate the mesh using hyperbolic grid generator
    coord.malloc((nid+1)*(njd+1), 3);
    IP.malloc(2*(nid-1));
    A.malloc(2*2*nid);
    B.malloc(2*2*nid);
    C.malloc(2*2*nid);
    RHS.malloc(2*nid);
    Z.malloc(2*2*(nid-1));
    ZA.malloc(2*(nid-2));
    vol.malloc(nid*njd);
  
    k6hyper2d4_(nid, njd, 
                fd->begin(posyd),
                coord1.begin(1), coord1.begin(2), coord1.begin(3),
                type,
                coord.begin(1), coord.begin(2), coord.begin(3),
                IP.begin(), A.begin(), B.begin(), C.begin(),
                RHS.begin(), Z.begin(), ZA.begin(), vol.begin());
   
    coord.reAllocMat(nid*njd, 3);
    // Array Creation 
    RELEASESHAREDS(array, f);
    RELEASESHAREDS(arrayd, fd);
    PyObject* tpl = 
      K_ARRAY::buildArray3(coord, varString0, nid, njd, 1, api);
    delete [] varString0;
    return tpl;
  }
  else if (res == 2 || resd == 2)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resd, arrayd, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "hyper2d4: not used for unstructured arrays.");
    return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "hyper2d4: unrecognised type of array.");
    return NULL;
  }
}    
