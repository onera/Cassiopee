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

// Adapt a 27-tree

#include "generator.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Adapt the octree */
//=============================================================================
PyObject* K_GENERATOR::adaptOctree3(PyObject* self, PyObject* args)
{
  PyObject *octree, *indica;
  if (!PyArg_ParseTuple(args, "OO", &octree, &indica)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(octree, varString, f, ni, nj, nk, cn, 
                                    eltType, true);
  if ( res != 2 ) 
  {
    if ( res == 1 ) RELEASESHAREDS(octree, f); 
    PyErr_SetString(PyExc_TypeError,
                    "adaptOctree3: array must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType, "HEXA") != 0 && strcmp(eltType, "QUAD") != 0 )
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "adaptOctree3: the octree must be HEXA or QUAD.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "adaptOctree3: the octree must contain coordinates.");
    return NULL;  
  }
  posx++; posy++; posz++;
 
  // Check indicator
  E_Int nii, nji, nki;
  FldArrayF* fi;
  FldArrayI* cni;
  char* varStringi;
  char* eltTypei;
  E_Int resi = K_ARRAY::getFromArray(indica, varStringi, fi, 
                                     nii, nji, nki, cni, eltTypei, true);
  if ( resi != 1 && resi != 2 ) 
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "adaptOctree3: indicator array must be structured.");
    return NULL;
  }
  E_Int posi = K_ARRAY::isNamePresent("indicator", varStringi);
  if ( posi == -1 ) 
  { 
    RELEASESHAREDB(resi, indica, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: adaptOctree3: no refinement indicator in the octree. Nothing done."); 
    return octree;
  }
  posi++;
  if ( fi->getSize() != cn->getSize() ) 
  {
    RELEASESHAREDB(resi, indica, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: adaptOctree3: refinement indicator size must be equal to the number of elements. Nothing done."); 
    return octree;
  }
  /*-----------------------------*/
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  E_Float* indict = fi->begin(posi);
  E_Int nelts = cn->getSize(); E_Int nvert = cn->getNfld();
  E_Int npts = f->getSize(); 
  FldArrayF* fo = new FldArrayF(npts, 3);
  FldArrayI* cno = new FldArrayI(nelts, nvert);
  FldArrayF* indicout = new FldArrayF(nelts,1); indicout->setAllValuesAtNull();
  E_Int no = 0; E_Int eto = 0;
  E_Float indic;

  // calcul de la bbox de l'octree
  E_Float xmin = K_CONST::E_MAX_FLOAT; E_Float xmax =-K_CONST::E_MAX_FLOAT;
  E_Float ymin = K_CONST::E_MAX_FLOAT; E_Float ymax =-K_CONST::E_MAX_FLOAT;
  E_Float zmin = K_CONST::E_MAX_FLOAT; E_Float zmax =-K_CONST::E_MAX_FLOAT;
  FldArrayF bbox(nelts,6);
  K_COMPGEOM::boundingBoxOfUnstrCells(*cn, xt, yt, zt, bbox);
  E_Float* xmint = bbox.begin(1); E_Float* xmaxt = bbox.begin(4);
  E_Float* ymint = bbox.begin(2); E_Float* ymaxt = bbox.begin(5);
  E_Float* zmint = bbox.begin(3); E_Float* zmaxt = bbox.begin(6);
  for (E_Int et = 0; et < nelts; et++)
  {
    xmin = K_FUNC::E_min(xmin,xmint[et]); xmax = K_FUNC::E_max(xmax,xmaxt[et]);
    ymin = K_FUNC::E_min(ymin,ymint[et]); ymax = K_FUNC::E_max(ymax,ymaxt[et]);
    zmin = K_FUNC::E_min(zmin,zmint[et]); zmax = K_FUNC::E_max(zmax,zmaxt[et]);
  }
  bbox.malloc(0,1);

  // calcul du pas sur chaque grille 
  FldArrayF dht(nelts);
  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
  E_Int ind1, ind2;
  for (E_Int et = 0; et < nelts; et++)
  {
    ind1 = cn1[et]-1; ind2 = cn2[et]-1;
    dht[et] = xt[ind2]-xt[ind1];
  }
  vector< vector<E_Int> > cEEN(nelts);
  E_Int ok = getNeighbourElts(npts, xt, yt, zt, *cn, cEEN); 
  if ( ok == 0 ) 
  {
    RELEASESHAREDB(resi, indica, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: adaptOctree3: cannot compute connectivity. Check connectivity."); 
    return octree;
  }
  // modification de l indicateur : si une cellule de niveau l a un indicateur a 1
  // et une cellule adjacente est de niveau l+1, indicateur < 1 : mettre indicateur(l+1) a 1
  // pour un saut de 2*h au pire si raffinement entre deux elements adjacents (reequilibrage)
  modifyIndicator(nelts, dht, cEEN, posi, indict);
  //short deraf = 0; // pour savoir s'il faut reequilibrer
  //raffinement / deraffinement
  FldArrayIS dejaVu(nelts); dejaVu.setAllValuesAtNull();
  for (E_Int et = 0; et < nelts; et++)
  {
    if ( dejaVu[et] == 0 ) 
    {
      indic = indict[et]; 
      if ( indic > 0. ) 
      {
        splitElement27(et, npts, indic, *cn, xt, yt, zt, *cno, *fo, 
                       *indicout, no, eto);
        dejaVu[et] = 1; 
      }
      else if ( indic < 0. ) 
      {                     
        //si merge : dejaVu de et et ses voisins mis a 1
        ok = mergeOctreeElement27(et, npts, indic, *cn, xmin, ymin, zmin, xmax, ymax, zmax,
                                  xt, yt, zt, dht.begin(1), indict,
                                  *cno, *fo, *indicout, no, eto, dejaVu);
        // if (ok == 1) deraf = 1;
      }   
    }
  }
  // elements restants : a conserver tels quels
  for (E_Int et = 0; et < nelts; et++)
  {
    if ( dejaVu[et] == 0 ) 
    {
      for (E_Int nov = 1; nov <= nvert; nov++)
      {
        E_Int ind = (*cn)(et, nov)-1;
        (*fo)(no,1) = xt[ind]; (*fo)(no,2) = yt[ind]; (*fo)(no,3) = zt[ind];
        no++; (*cno)(eto, nov) = no; 
        if ( no == fo->getSize() ) 
        {
          fo->reAllocMat(fo->getSize()+npts, 3);
        }
      }
      (*indicout)[eto]=0.;
      eto++;
      if ( eto == cno->getSize() )
      {
        cno->reAllocMat(cno->getSize()+nelts,nvert); 
        indicout->reAllocMat(eto,1);
      }
      dejaVu[et] = 1;
    }
  }

  // Nettoyage
  for (E_Int v = 0; v < nelts; v++) cEEN[v].clear();
  cEEN.clear();

  // Realloc
  fo->reAllocMat(no,3); cno->reAllocMat(eto,nvert); indicout->reAllocMat(eto,1);
  K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, eltType, *fo, *cno);

  // Sortie
  RELEASESHAREDU(octree, f, cn);  RELEASESHAREDB(resi, indica, fi, cni);
  PyObject* l = PyList_New(0); 
  PyObject* tpl = K_ARRAY::buildArray(*fo, "x,y,z", *cno, -1, eltType);
  delete fo;
  PyList_Append(l, tpl); Py_DECREF(tpl);
  char eltType2[8];
  if (strcmp(eltType, "QUAD")  == 0) strcpy(eltType2, "QUAD*");
  else strcpy(eltType2, "HEXA*");
  PyObject* tpl2 = K_ARRAY::buildArray(*indicout, "indicator", *cno, 
                                       -1, eltType2);
  delete cno; delete indicout;
  PyList_Append(l, tpl2); Py_DECREF(tpl2);
  return l;
}
//=============================================================================
/* Split un element en 9/27 elements fils */
//=============================================================================
void K_GENERATOR::splitElement27(
  E_Int et, E_Int npts, E_Float indic,
  FldArrayI& cn, E_Float* xt, E_Float* yt, E_Float* zt, 
  FldArrayI& cno, FldArrayF& fo, FldArrayF& indicout, 
  E_Int& no, E_Int& eto)
{
  E_Int nelts = cno.getSize(); E_Int nvert = cno.getNfld();
  E_Int nptso = fo.getSize(); 
  if ( no >= nptso-217) 
  {nptso = nptso+217; fo.reAllocMat(nptso,3);}
  if ( eto >= cno.getSize()-28) 
  {nelts = nelts+28; cno.reAllocMat(nelts,nvert); indicout.reAllocMat(nelts,1);}

  E_Float* xo = fo.begin(1); E_Float* yo = fo.begin(2); E_Float* zo = fo.begin(3);
  E_Float* indico = indicout.begin(1);
  if ( nvert == 4 ) // QUAD
  {
    E_Int indmin = cn(et,1)-1; E_Int indmax = cn(et,3)-1; 
    E_Float z = zt[indmin];
    E_Float xmin = xt[indmin]; E_Float ymin = yt[indmin];
    E_Float xmax = xt[indmax]; E_Float ymax = yt[indmax];
    E_Float dhs3 = (xmax-xmin)/3.;

    // 1er quad
    xo[no] = xmin; yo[no] = ymin; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;
    
    // 2eme quad
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 3eme quad
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 4e quad
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;
    
    // 5e quad
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 6e quad
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 7e quad
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymax; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;
    
    // 8e quad
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 9e quad
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = z; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymax; zo[no] = z; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = z; no++; cno(eto,4) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

  }
  else // HEXA
  {
    E_Int indmin = cn(et,1)-1;  E_Int indmax = cn(et,7)-1;
    E_Float xmin = xt[indmin]; E_Float ymin = yt[indmin]; E_Float zmin = zt[indmin];
    E_Float xmax = xt[indmax]; E_Float ymax = yt[indmax]; E_Float zmax = zt[indmax];
    E_Float dhs3 = (xmax-xmin)/3.;

    // 1er HEXA
    xo[no] = xmin; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,8) = no;
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 2eme HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 3eme HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 4eme HEXA
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 5eme HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 6e HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 7e HEXA
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;
    
    // 8e HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 9e HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 2e LAYER
    //----------
    // 10e HEXA
    xo[no] = xmin; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no;
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 11e HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 12e HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 13e HEXA
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 14e HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 15e HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 16e HEXA
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;
    
    // 17e HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 18e HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin+dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin+dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    //3e LAYER
    //--------
    // 19e HEXA
    xo[no] = xmin; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,8) = no;
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 20e HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 21e HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 22e HEXA
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 23e HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 24e HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin+dhs3; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 25e HEXA
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;
    
    // 26e HEXA
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin+dhs3; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;

    // 27e HEXA
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,1) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmin+2.*dhs3; no++; cno(eto,2) = no; 
    xo[no] = xmax; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,3) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmin+2.*dhs3; no++; cno(eto,4) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,5) = no; 
    xo[no] = xmax; yo[no] = ymin+2.*dhs3; zo[no] = zmax; no++; cno(eto,6) = no; 
    xo[no] = xmax; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,7) = no; 
    xo[no] = xmin+2.*dhs3; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,8) = no; 
    indico[eto] = K_FUNC::E_max(0., indic-1.); eto++;
  }
}

//=============================================================================
/* Fusionne un element avec 8/26 voisins valides */
//=============================================================================
E_Int
K_GENERATOR::mergeOctreeElement27(E_Int et, E_Int npts, E_Float indic,
                                  K_FLD::FldArrayI& cn, 
                                  E_Float xs, E_Float ys, E_Float zs,
                                  E_Float xe, E_Float ye, E_Float ze,
                                  E_Float* xt, E_Float* yt, E_Float* zt, 
                                  E_Float* dht, E_Float* indict, 
                                  K_FLD::FldArrayI& cno, K_FLD::FldArrayF& fo, 
                                  K_FLD::FldArrayF& indicout,
                                  E_Int& no, E_Int& eto, 
                                  K_FLD::FldArrayIS& dejaVu)
{
  E_Float eps = 1.e-10;
  // Recherche des elements voisins de meme niveau
  E_Int nvert = cno.getNfld(); E_Int nelts = cno.getSize();
  E_Int nptso = fo.getSize();
  if ( no >= nptso-217) 
  {nptso = nptso+217; fo.reAllocMat(nptso,3);}
  if ( eto >= cno.getSize()-28) 
  {nelts = nelts+28; cno.reAllocMat(nelts,nvert); indicout.reAllocMat(nelts,1); }

  E_Float* xo = fo.begin(1); E_Float* yo = fo.begin(2); E_Float* zo = fo.begin(3);
  E_Float* indico = indicout.begin(1);
  E_Float dh = dht[et];
  E_Float xmin, ymin, xmax, ymax, zmin, zmax;

  std::vector<E_Int> candidats;
  getValidNgbrsForMerge(et, indict, dht, xs, ys, zs, xe, ye, ze, 
                        xt, yt, zt, dejaVu, cn, candidats);
  E_Int ncandidats = candidats.size();

  // 1- verifier si on a au moins un nb de candidats suffisant (3 pour le 2D).
  // 2- sinon : 8 QUAD candidats : forment-ils un carre ? Immediat 
  // 3- plus de nmax candidats : impossible, sinon ca sent le bug qqpart...
  E_Int nmax = 8; if (nvert == 8) nmax = 26;
  E_Int indA1, indB1, indC1, indE1;
  E_Int indA2, indB2, indC2, indE2;
  
  if (ncandidats < nmax ) 
  {
    E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
    indA1 = cn1[et]-1; indB1 = cn2[et]-1; 
    return -1; //cas 1
  }
  else if ( ncandidats == nmax ) //cas 2
  {
    if ( nvert == 4 ) // QUAD
    {
      E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3); //E_Int* cn4 = cn.begin(4);
      indA1 = cn1[et]-1; indB1 = cn2[et]-1; indC1 = cn3[et]-1;
      xmin = xt[indA1]; xmax = xt[indB1]; 
      ymin = yt[indA1]; ymax = yt[indC1];
      zmin = zt[indA1]; zmax = zmin;
      for (E_Int v = 0; v < ncandidats; v++)
      {
        E_Int et2 = candidats[v]; 
        indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1; //indD2 = cn4[et2]-1;
        xmin = K_FUNC::E_min(xmin, xt[indA2]); xmax = K_FUNC::E_max(xmax, xt[indB2]); 
        ymin = K_FUNC::E_min(ymin, yt[indA2]); ymax = K_FUNC::E_max(ymax, yt[indC2]); 
      }
      if (K_FUNC::fEqualZero(3.*dh-(xmax-xmin),eps) == false || 
          K_FUNC::fEqualZero(3.*dh-(ymax-ymin),eps) == false ) return -2;

      // construction du quad
      xo[no] = xmin; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,1) = no; 
      xo[no] = xmax; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,2) = no; 
      xo[no] = xmax; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,3) = no; 
      xo[no] = xmin; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,4) = no; 
      indico[eto] = K_FUNC::E_min(0.,indic+1);
      eto++;
      for (E_Int v = 0; v < ncandidats; v++)
      {E_Int et2 = candidats[v]; dejaVu[et2] = 1;}
      dejaVu[et] = 1;
      return 1;
    }
    else // HEXA
    {
      E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2); 
      E_Int* cn3 = cn.begin(3); //E_Int* cn4 = cn.begin(4);
      E_Int* cn5 = cn.begin(5); //E_Int* cn6 = cn.begin(6); 
      //E_Int* cn7 = cn.begin(7); E_Int* cn8 = cn.begin(8);
      indA1 = cn1[et]-1; indB1 = cn2[et]-1; indC1 = cn3[et]-1; //indD1 = cn4[et]-1;
      indE1 = cn5[et]-1; //indF1 = cn6[et]-1; indG1 = cn7[et]-1; indH1 = cn8[et]-1;
      xmin = xt[indA1]; xmax = xt[indB1]; 
      ymin = yt[indA1]; ymax = yt[indC1];
      zmin = zt[indA1]; zmax = zt[indE1];
    
      for (E_Int v = 0; v < ncandidats; v++)
      {
        E_Int et2 = candidats[v]; 
        indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1; indE2 = cn5[et2]-1;
        xmin = K_FUNC::E_min(xmin, xt[indA2]); xmax = K_FUNC::E_max(xmax, xt[indB2]); 
        ymin = K_FUNC::E_min(ymin, yt[indA2]); ymax = K_FUNC::E_max(ymax, yt[indC2]);
        zmin = K_FUNC::E_min(zmin, zt[indA2]); zmax = K_FUNC::E_max(zmax, zt[indE2]);
      }
      
      if (K_FUNC::fEqualZero(3.*dh-(xmax-xmin),eps) == false || 
          K_FUNC::fEqualZero(3.*dh-(ymax-ymin),eps) == false ||
          K_FUNC::fEqualZero(3.*dh-(zmax-zmin),eps) == false) return -3;
  
      // construction de l hexa
      xo[no] = xmin; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,1) = no; 
      xo[no] = xmax; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,2) = no; 
      xo[no] = xmax; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,3) = no; 
      xo[no] = xmin; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,4) = no; 
      xo[no] = xmin; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,5) = no; 
      xo[no] = xmax; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,6) = no; 
      xo[no] = xmax; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,7) = no; 
      xo[no] = xmin; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,8) = no; 
      indico[eto] = K_FUNC::E_min(0.,indic+1);
      eto++;
      for (E_Int v = 0; v < ncandidats; v++)
      {E_Int et2 = candidats[v]; dejaVu[et2] = 1;}
      dejaVu[et] = 1;
      return 1;
    }
  }
  printf("Error : mergeElement : too many candidates (%d) to merge with element %d. Check the octree.\n", ncandidats, et+1);
  
  return -1;
}
//=============================================================================
