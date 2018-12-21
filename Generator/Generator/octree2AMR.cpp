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

#include "generator.h"

using namespace K_FLD;
using namespace std;
using namespace K_CONST;

//=============================================================================
/* Generates a set of HEXA grids of AMR levels from an octree HEXA mesh */
//=============================================================================
PyObject* K_GENERATOR::octree2AMR(PyObject* self, PyObject* args)
{
  PyObject* octree;
  E_Int vmin;
#ifdef E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "Ol", &octree, &vmin)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Oi", &octree, &vmin)) return NULL;
#endif
  if (vmin < 1) 
  { 
    PyErr_SetString(PyExc_TypeError, 
                    "octree2AMR: vmin must be greater than 0.");
    return NULL;
  }
  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f;
  char* varString; char* eltType;
  K_FLD::FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray(octree, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  E_Int dim = 3;
  if (res != 2) 
  {
    if (res == 1) RELEASESHAREDS(octree, f); 
    PyErr_SetString(PyExc_TypeError, 
                    "octree2AMR: array must be unstructured.");
    return NULL;   
  }
  if (strcmp(eltType,"HEXA") == 0) dim = 3;
  else if (strcmp(eltType,"QUAD") == 0) dim = 2;
  else 
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError, 
                    "octree2AMR: array must be HEXA or QUAD.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "octree2AMR: coordinates not found in array.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  vector<FldArrayF*> vectOfCartGrids;
  /* Etape 1 : detection du niveau le plus fin*/ 
  //calcul de dh
  E_Int nelts = cn->getSize(); E_Int nvert = cn->getNfld();
  E_Int* cnt1 = cn->begin(1); E_Int* cnt2 = cn->begin(2);
  FldArrayF dht(nelts); E_Float dhmin = K_CONST::E_MAX_FLOAT;E_Float dhmax = 0.;
  E_Float* xo = f->begin(posx);   E_Float* yo = f->begin(posy); E_Float* zo = f->begin(posz); 
  for (E_Int et = 0; et < nelts; et++)
  {
    E_Int ind1 = cnt1[et]-1; E_Int ind2 = cnt2[et]-1;
    dht[et] = xo[ind2]-xo[ind1]; dhmin = K_FUNC::E_min(dhmin, dht[et]);
    dhmax = K_FUNC::E_max(dhmax, dht[et]);
  }  
  //niveau le plus fin 
  ni = vmin; nj = vmin; nk = vmin; if ( dim == 2 ) nk = 1;
  E_Int ninj = ni*nj; E_Int ninjnk = ninj*nk;
  E_Int ind1, ind2, ind;
  E_Float xmin, ymin, xmax, ymax, zmin, zmax, dh;
  E_Float vmini = 1./(vmin-1);
  FldArrayF indic(nelts); indic.setAllValuesAtNull();
  E_Float* indict = indic.begin();
  for (E_Int et = 0; et < nelts; et++)
  {
    if (K_FUNC::fEqualZero(dht[et]-dhmin) == true) 
    {
      indict[et] = -1.;
      FldArrayF* coord = new FldArrayF(ninjnk,3); 
      E_Float* xt = coord->begin(1); E_Float* yt = coord->begin(2); E_Float* zt = coord->begin(3);
      ind1 = cnt1[et]-1; ind2 = cnt2[et]-1;
      xmin = xo[ind1]; xmax = xo[ind2]; ymin = yo[ind1]; zmin = zo[ind1];
      dh = (xmax-xmin)*vmini;
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            ind = i + j*ni + k*ninj; 
            xt[ind] = xmin + i*dh; yt[ind] = ymin + j*dh; zt[ind] = zmin + k*dh;
          }
      vectOfCartGrids.push_back(coord);       
    }
  }
  /*---------------------*/
  /* niveau l quelconque */
  /*---------------------*/
  E_Int npts; 
  E_Float dhloc = dhmin;
  FldArrayF f2(f->getSize(),3); FldArrayI cn2 = *cn; FldArrayF indico = indic;
  f2.setOneField(*f,posx,1); f2.setOneField(*f,posy,2);f2.setOneField(*f,posz,3);// copy 
  E_Int level = 1;
  FldArrayF fo; FldArrayI cno;
  FldArrayF indicout(nelts,1); indicout.setAllValuesAtNull();
  while ( dhloc < dhmax ) 
  {
    npts = f2.getSize(); nelts = cn2.getSize();
    FldArrayIS dejaVu(nelts); dejaVu.setAllValuesAt(0);
    fo.malloc(npts, 3); cno.malloc(nelts, nvert); 
    E_Float* indicto = indico.begin();
    // calcul de la bbox de l octree
    E_Float* xt = f2.begin(1); E_Float* yt = f2.begin(2); E_Float* zt = f2.begin(3); 
    xmin = K_CONST::E_MAX_FLOAT; xmax =-K_CONST::E_MAX_FLOAT;
    ymin = K_CONST::E_MAX_FLOAT; ymax =-K_CONST::E_MAX_FLOAT;
    zmin = K_CONST::E_MAX_FLOAT; zmax =-K_CONST::E_MAX_FLOAT;
    FldArrayF bbox(nelts,6); K_COMPGEOM::boundingBoxOfUnstrCells(cn2, xt, yt, zt, bbox);   
    E_Float* xmint = bbox.begin(1); E_Float* xmaxt = bbox.begin(4);
    E_Float* ymint = bbox.begin(2); E_Float* ymaxt = bbox.begin(5);
    E_Float* zmint = bbox.begin(3); E_Float* zmaxt = bbox.begin(6);
    for (E_Int et = 0; et < nelts; et++)
    {
      xmin = K_FUNC::E_min(xmin,xmint[et]); xmax = K_FUNC::E_max(xmax,xmaxt[et]);
      ymin = K_FUNC::E_min(ymin,ymint[et]); ymax = K_FUNC::E_max(ymax,ymaxt[et]);
      zmin = K_FUNC::E_min(zmin,zmint[et]); zmax = K_FUNC::E_max(zmax,zmaxt[et]);
    }

    //fusion des elements de niveau l-1
    E_Int no = 0; E_Int eto = 0; 
    for (E_Int et = 0; et < nelts; et++)
    {
      if (dejaVu[et] == 0 && indico[et] == -1.) 
      {
        mergeOctreeElement(et, npts, indico[et], cn2, xmin, ymin, zmin, xmax, ymax, zmax,
			   xt, yt, zt, dht.begin(1), indicto,
			   cno, fo, indicout, no, eto, dejaVu);
      }
    }
    indico.setAllValuesAt(-1.);
    // elements restants
    for (E_Int et = 0; et < nelts; et++)
    {
      if (dejaVu[et] == 0) 
      {
        for (E_Int nov = 1; nov <= nvert; nov++)
        {
          ind = cn2(et, nov)-1;
          fo(no,1) = xt[ind]; fo(no,2) = yt[ind]; fo(no,3) = zt[ind];
          no++; cno(eto,nov) = no; 
          if ( no == fo.getSize() ) fo.reAllocMat(fo.getSize()+npts,3);
        }
        if ( K_FUNC::fEqualZero(dht[et]-dhloc) == false) indico[eto] = 0.;
        eto++; 
        if (eto == cno.getSize()){cno.reAllocMat(cno.getSize()+nelts,nvert); indico.resize(indico.getSize()+nelts);}
      }
    }
    fo.reAllocMat(no,3); cno.reAllocMat(eto,nvert); indico.resize(eto);

    // construction des grilles cartesiennes a partir du niveau l
    E_Int neltso = eto; E_Int* cno1 = cno.begin(1); E_Int* cno2 = cno.begin(2);
    xo = fo.begin(1); yo = fo.begin(2); zo = fo.begin(3);
    for (E_Int eto = 0; eto < neltso; eto++)
    {
      if ( indico[eto] == -1. )
      {
        FldArrayF* coord = new FldArrayF(ninjnk,3); 
        E_Float* xt = coord->begin(1); E_Float* yt = coord->begin(2); E_Float* zt = coord->begin(3);
        ind1 = cno1[eto]-1; ind2 = cno2[eto]-1;
        xmin = xo[ind1]; xmax = xo[ind2]; ymin = yo[ind1]; zmin = zo[ind1];
        dh = (xmax-xmin)*vmini;
        for (E_Int k = 0; k < nk; k++)
          for (E_Int j = 0; j < nj; j++)
            for (E_Int i = 0; i < ni; i++)
            {
              ind = i + j*ni + k*ninj; 
              xt[ind] = xmin + i*dh; yt[ind] = ymin + j*dh; zt[ind] = zmin + k*dh;
            }
        vectOfCartGrids.push_back(coord);       
      }
    }

    cn2 = cno; f2 = fo;
    //calcul de dh
    dht.resize(neltso); dhmin = K_CONST::E_MAX_FLOAT; dhmax = 0.;
    for (E_Int eto = 0; eto < neltso; eto++)
    {
      ind1 = cno1[eto]-1; ind2 = cno2[eto]-1;
      dht[eto] = xo[ind2]-xo[ind1]; dhmin = K_FUNC::E_min(dhmin, dht[eto]);
      dhmax = K_FUNC::E_max(dhmax, dht[eto]);
    } 
    dhloc = dhmin; // pas d espace du niveau courant
    level++;
  }// while dhloc < dhmax 

  //dernier  niveau
  E_Int* cno1 = cno.begin(1); E_Int* cno2 = cno.begin(2);
  for (E_Int eto = 0; eto < cno.getSize(); eto++)
  {
    FldArrayF* coord = new FldArrayF(ninjnk,3); 
    E_Float* xt = coord->begin(1); E_Float* yt = coord->begin(2); E_Float* zt = coord->begin(3);
    ind1 = cno1[eto]-1; ind2 = cno2[eto]-1;
    xmin = xo[ind1]; xmax = xo[ind2]; ymin = yo[ind1]; zmin = zo[ind1];
    dh = (xmax-xmin)*vmini;
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j*ni + k*ninj; 
          xt[ind] = xmin + i*dh; yt[ind] = ymin + j*dh; zt[ind] = zmin + k*dh;
        }
    vectOfCartGrids.push_back(coord);       
  }
  
  /* sortie */
  RELEASESHAREDU(octree, f, cn); 
  PyObject* l = PyList_New(0);
  E_Int nzones = vectOfCartGrids.size();
  for (E_Int v = 0; v < nzones; v++)
  {
    PyObject* tpl = K_ARRAY::buildArray(*vectOfCartGrids[v], "x,y,z", ni, nj, nk);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  // nettoyage
  for (E_Int v = 0; v < nzones; v++) delete vectOfCartGrids[v];
  return l;
}
//=============================================================================

