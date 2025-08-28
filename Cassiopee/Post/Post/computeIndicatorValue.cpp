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
# include "post.h"
# include "Nuga/include/BbTree.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Projection du champ indicateur sur le maillage octree non structure a 
   partir du champ indicateur sur le maillage structure de depart. Retourne 
   la valeur absolue maximale du champ */
//=============================================================================
PyObject* K_POST::computeIndicatorValue(PyObject* self, PyObject* args)
{
  PyObject *octree, *zones, *fieldA; 
  if (!PyArg_ParseTuple(args, "OOO", &octree, &zones, &fieldA)) return NULL;

  // Verif octree HEXA/QUAD
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(octree, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeIndicatorValue: array must be unstructured.");
    RELEASESHAREDB(res, octree, f, cn); return NULL;   
  }

  E_Int posxo = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posyo = K_ARRAY::isCoordinateYPresent(varString);
  E_Int poszo = K_ARRAY::isCoordinateZPresent(varString);
  if (posxo == -1 || posyo == -1 || poszo == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeIndicatorValue: coordinates not found in array.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  posxo++; posyo++; poszo++;

  // coordonnees des grilles 
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit;
  vector<E_Int> njt; 
  vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = true;
  E_Bool skipDiffVars = true;
  K_ARRAY::getFromArrays(
    zones, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = structF.size();
  if (nzones == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeIndicatorValue: no valid zone detected.");
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);      
    RELEASESHAREDU(octree, f, cn); return NULL; 
  }
  // valeur du champ pour l indicateur
  vector<char*> structVarString2;
  vector<char*> unstrVarString2;
  vector<FldArrayF*> structF2;
  vector<FldArrayF*> unstrF2;
  vector<E_Int> nit2;
  vector<E_Int> njt2; 
  vector<E_Int> nkt2;
  vector<FldArrayI*> cnt2;
  vector<char*> eltTypet2;
  vector<PyObject*> objst2, objut2;
  skipNoCoord = false;

  K_ARRAY::getFromArrays(
    fieldA, resl, structVarString2, unstrVarString2,
    structF2, unstrF2, nit2, njt2, nkt2, cnt2, eltTypet2, objst2, objut2, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones2 = structF2.size();
  if (nzones2 != nzones) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeIndicatorValue: no valid array for indicator field detected.");
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);      
    for (E_Int nos = 0; nos < nzones2; nos++)
      RELEASESHAREDS(objst2[nos], structF2[nos]);      
    RELEASESHAREDU(octree, f, cn); return NULL; 
  }
  for (E_Int noz = 0; noz < nzones2; noz++)
  {
    if (structF2[noz]->getNfld() != 1) 
    {
      printf("Warning: computeIndicatorValue: several fields defined. First one is selected.\n");
    }
  }
  // Verification de posxi, posyi, poszi dans arrays
  E_Int posx=1, posy=2, posz=3;
  for (E_Int i = 0; i < nzones; i++)
  {
    E_Int posxi = K_ARRAY::isCoordinateXPresent(structVarString[i]); posxi++;
    E_Int posyi = K_ARRAY::isCoordinateYPresent(structVarString[i]); posyi++;
    E_Int poszi = K_ARRAY::isCoordinateZPresent(structVarString[i]); poszi++;
    if (posxi != posx || posyi != posy || poszi != posz) 
    {
      PyErr_SetString(PyExc_TypeError,
                    "computeIndicatorValue: invalid coordinates in zones.");
      for (E_Int nos = 0; nos < nzones; nos++)
        RELEASESHAREDS(objst[nos], structF[nos]);      
      for (E_Int nos = 0; nos < nzones2; nos++)
        RELEASESHAREDS(objst2[nos], structF2[nos]);      
      RELEASESHAREDU(octree, f, cn); return NULL;
    }
  }
  // Construction du bbtree 
  E_Float minB[3];  E_Float maxB[3];
  typedef K_SEARCH::BoundingBox<3>  BBox3DType;
  FldArrayF bboxes(nzones,6);
  vector<BBox3DType*> boxes(nzones);// liste des bbox
  for (E_Int v = 0; v < nzones; v++)
  {
    K_COMPGEOM::boundingBox(nit[v], njt[v], nkt[v], posx, posy, posz, *structF[v], 
                            bboxes(v,1), bboxes(v,2), bboxes(v,3),bboxes(v,4),bboxes(v,5),bboxes(v,6));
  }
  E_Float* xmin = bboxes.begin(1); E_Float* xmax = bboxes.begin(4);
  E_Float* ymin = bboxes.begin(2); E_Float* ymax = bboxes.begin(5);
  E_Float* zmin = bboxes.begin(3); E_Float* zmax = bboxes.begin(6);
  for (E_Int v = 0; v < nzones; v++)
  {
    minB[0] = xmin[v]; minB[1] = ymin[v]; minB[2] = zmin[v];
    maxB[0] = xmax[v]; maxB[1] = ymax[v]; maxB[2] = zmax[v]; 
    boxes[v] = new BBox3DType(minB, maxB);
  }
  K_SEARCH::BbTree3D bbtree(boxes);

  //======================
  E_Int nptso = f->getSize();
  FldArrayF* field = new FldArrayF(nptso); field->setAllValuesAtNull();
  E_Float* fp = field->begin();
  E_Int nvert = cn->getNfld();
  E_Int indmin, indmax;
  // Determination des zones intersectant l elt
  E_Float* xo = f->begin(posxo);
  E_Float* yo = f->begin(posyo);
  E_Float* zo = f->begin(poszo);
  vector<E_Int> indicesBB;
  E_Int nindicesBB;
  E_Float* xmint = bboxes.begin(1);
  E_Float* ymint = bboxes.begin(2);
  E_Float* zmint = bboxes.begin(3);
  E_Float* xmaxt = bboxes.begin(4);

  E_Int found = 0; 
  E_Int ind, ic, jc, kc, novert; E_Float valmax;
  E_Int icmin, icmax, jcmin, jcmax, kcmin, kcmax;
  E_Int nelts = cn->getSize();
  if (nvert == 8) // HEXA
  {
    E_Int* cn1 = cn->begin(1);
    E_Int* cn7 = cn->begin(7);
    for (E_Int et = 0; et <nelts; et++)
    {
      indmin = cn1[et]-1; indmax = cn7[et]-1;
      minB[0] = xo[indmin]; minB[1] = yo[indmin]; minB[2] = zo[indmin];
      maxB[0] = xo[indmax]; maxB[1] = yo[indmax]; maxB[2] = zo[indmax];
      indicesBB.clear();
      bbtree.getOverlappingBoxes(minB, maxB, indicesBB);
      nindicesBB = indicesBB.size();
      for (E_Int nov = 0; nov < nindicesBB; nov++)
      {
        found = 0; 
        E_Int noblk = indicesBB[nov];
        //npts = structF[noblk]->getSize();
        // test si la grille contient 8 pts coincidents avec et
        E_Float* ft = structF2[noblk]->begin();
        E_Int nib = nit[noblk]; E_Int njb = njt[noblk]; E_Int nkb = nkt[noblk];
        E_Int nic = K_FUNC::E_max(1,nib-1);
        E_Float dh = (xmaxt[noblk]-xmint[noblk])/nic;
        E_Float dhi = K_CONST::ONE/dh;
        // pt1 de l element
        icmin = 10000000; icmax = -1; jcmin = icmin; jcmax =-1; kcmin=icmin; kcmax=-1;
        for (novert = 1; novert <= nvert; novert++)
        {
          ind = (*cn)(et,novert)-1;          
          ic = E_Int((xo[ind]-xmint[noblk])*dhi);
          jc = E_Int((yo[ind]-ymint[noblk])*dhi);
          kc = E_Int((zo[ind]-zmint[noblk])*dhi);
          if ( ic < 0 || ic >= nib || jc < 0 || jc >= njb || kc < 0 || kc >= nkb ) goto next;
          icmin  = K_FUNC::E_min(ic,icmin); icmax  = K_FUNC::E_max(ic,icmax);
          jcmin  = K_FUNC::E_min(jc,jcmin); jcmax  = K_FUNC::E_max(jc,jcmax);
          kcmin  = K_FUNC::E_min(kc,kcmin); kcmax  = K_FUNC::E_max(kc,kcmax);
          found++; 
        }
    
        if (found >= 8) // les sommets sont interieurs de la grille cartesienne
        {
          valmax = 0.;
          for (kc = kcmin; kc <= kcmax; kc++)
            for (jc = jcmin; jc <= jcmax; jc++)
              for (ic = icmin; ic <= icmax; ic++)
              {
                ind = ic + jc * nib + kc*nib*njb;
                valmax = K_FUNC::E_max(valmax,ft[ind]);
              }
          for (novert = 1; novert <= nvert; novert++)
          {
            ind = (*cn)(et,novert)-1;       
            fp[ind] = K_FUNC::E_max(fp[ind],valmax);
          }
        }
        next:;
      }// boucle sur les zones intersectantes

    }// boucle sur les elts
  }
  else if ( nvert == 4 ) // QUAD
  {
    E_Int* cn1 = cn->begin(1);
    E_Int* cn3 = cn->begin(3);
    for (E_Int et = 0; et <nelts; et++)
    {
      indmin = cn1[et]-1; indmax = cn3[et]-1;
      minB[0] = xo[indmin]; minB[1] = yo[indmin]; minB[2] = zo[indmin];
      maxB[0] = xo[indmax]; maxB[1] = yo[indmax]; maxB[2] = zo[indmax];
      indicesBB.clear();
      bbtree.getOverlappingBoxes(minB, maxB, indicesBB);
      nindicesBB = indicesBB.size();
      for (E_Int nov = 0; nov < nindicesBB; nov++)
      {
        found = 0; E_Int noblk = indicesBB[nov];
        //npts = structF[noblk]->getSize();
        // test si la grille contient 4 pts coincidents avec et
        E_Float* ft = structF2[noblk]->begin();
        E_Int nib = nit[noblk]; E_Int njb = njt[noblk];
        E_Int nic = K_FUNC::E_max(1,nib-1);
        E_Float dh = (xmaxt[noblk]-xmint[noblk])/nic;
        E_Float dhi = K_CONST::ONE/dh;
        // pt1 de l element
        icmin = 10000000; icmax = -1; jcmin = icmin; jcmax =-1; 
        for (novert = 1; novert <= nvert; novert++)
        {
          ind = (*cn)(et,novert)-1;          
          ic = E_Int((xo[ind]-xmint[noblk])*dhi);
          jc = E_Int((yo[ind]-ymint[noblk])*dhi);
          if ( ic < 0 || ic >= nib || jc < 0 || jc >= njb  ) goto next2;
          icmin = K_FUNC::E_min(ic,icmin); icmax = K_FUNC::E_max(ic,icmax);
          jcmin = K_FUNC::E_min(jc,jcmin); jcmax = K_FUNC::E_max(jc,jcmax);
          found++; 
        }
    
        if ( found >= 4 ) // les sommets sont interieurs de la grille cartesienne
        {
          valmax = 0.;
          for (jc = jcmin; jc <= jcmax; jc++)
            for (ic = icmin; ic <= icmax; ic++)
            {
              ind = ic + jc * nib;
              valmax = K_FUNC::E_max(valmax,ft[ind]);
            }
          for (novert = 1; novert <= nvert; novert++)
          {
            ind = (*cn)(et,novert)-1;       
            fp[ind] = K_FUNC::E_max(fp[ind],valmax);
          }
        }
        next2:;
      }// boucle sur les zones intersectantes
      
    }// boucle sur les elts
  }

  vector<char*> vars;
  K_ARRAY::extractVars(structVarString2[0], vars); 
  PyObject* tpl = K_ARRAY::buildArray(*field, vars[0], *cn, -1, eltType, true);
  
  for (E_Int v = 0; v < nzones; v++) delete boxes[v];
  delete field;  RELEASESHAREDB(res, octree, f, cn); 
  for (size_t i = 0; i < vars.size(); i++) delete[] vars[i];
  vars.clear();
  
  for (E_Int nos = 0; nos < nzones; nos++)
    RELEASESHAREDS(objst[nos], structF[nos]);  
  for (E_Int nos = 0; nos < nzones2; nos++)
    RELEASESHAREDS(objst2[nos], structF2[nos]);
    
  return tpl;
}
