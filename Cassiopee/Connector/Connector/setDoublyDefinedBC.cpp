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
# include "connector.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* CL doublement definies : */
//=============================================================================
PyObject* K_CONNECTOR::setDoublyDefinedBC(PyObject* self, PyObject* args)
{
  PyObject *a1, *celln1, *arrays, *cellns, *range1;
  E_Int depth;
  if (!PYPARSETUPLE_(args, OOOO_ O_ I_,
                    &a1, &celln1, &arrays, &cellns, &range1, &depth))
  {
      return NULL;
  }
  
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setDoublyDefinedBC: third argument must be a list.");
    return NULL;
  }
  if ( PyList_Check(cellns) == 0 )
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setDoublyDefinedBC: fourth argument must be a list.");
    return NULL;
  }
  //Range des parois doublement definies
  if ( PyList_Check(range1) == 0 )
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setDoublyDefinedBC: fifth argument must be a list.");
    return NULL;
  }
  if( depth != 1 && depth != 2 ) 
  {
    printf("Warning: setDoublyDefinedBC: depth must be equal to 1 or 2. Set to 2.\n"); 
    depth = 2;
  }
  // Check array
  E_Int im, jm, km;
  FldArrayF* f;
  FldArrayI* cn0;
  char* varString;
  char* eltType0;
  E_Int res = K_ARRAY::getFromArray(a1, varString, f, im, jm, km, 
                                    cn0, eltType0, true); 
  if (res != 1) 
  {
    RELEASESHAREDB(res, a1, f, cn0);
    PyErr_SetString(PyExc_TypeError, 
                    "setDoublyDefinedBC: 1st argument not valid.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setDoublyDefinedBC: 1st arg must contain coordinates.");
    RELEASESHAREDS(a1, f);
    return NULL;
  }
  posx++; posy++; posz++;

  // Check array
  E_Int imc, jmc, kmc;
  FldArrayF* fc;
  char* varStringc;
  res = K_ARRAY::getFromArray(celln1, varStringc, fc, imc, jmc, kmc, 
                              cn0, eltType0, true); 
  if (res != 1) 
  {
    RELEASESHAREDS(a1, f); RELEASESHAREDB(res, celln1, fc, cn0);
    PyErr_SetString(PyExc_TypeError, 
                    "setDoublyDefinedBC: 2nd argument not valid.");
    return NULL;
  }
  E_Int posc = K_ARRAY::isCellNatureField2Present(varStringc);
  if (posc == -1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "setDoublyDefinedBC: 2nd arg must contain cellN variable.");
    RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); return NULL;
  }
  posc++;

  //check dimensions
  if ((km == 1) && (kmc == 0))
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setDoublyDefinedBC: implemented only in 3D.");
    RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); return NULL;
  }    
  res = checkDimensions(im, jm, km, imc, jmc, kmc);
  if (res == 0) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setDoublyDefinedBC: 1st and 2nd arguments not compliant.");
    RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); return NULL;
  }

  // Extract infos from coord arrays 
  // seulement arrays structures avec coordonnees ici 
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<K_FLD::FldArrayF*> structF;
  vector<K_FLD::FldArrayF*> unstrF;
  vector<E_Int> nit; 
  vector<E_Int> njt; 
  vector<E_Int> nkt;
  vector<K_FLD::FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = true;
  E_Bool skipDiffVars = true;
  res = K_ARRAY::getFromArrays(
    arrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, 
    skipUnstructured, true);
  E_Int ns = structF.size();
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setDoublyDefinedBC: 3rd argument is not valid."); 
    RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); 
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objs[is], structF[is]);
    return NULL;
  }  
  if (ns == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "setDoublyDefinedBC: no valid donor zone found."); 
    RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); 
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objs[is], structF[is]);
    return NULL;
  }
  // Extract infos from celln arrays : seulement structures
  vector<E_Int> resc;
  vector<char*> structVarStringc;
  vector<char*> unstrVarStringc;
  vector<K_FLD::FldArrayF*> structFc;
  vector<K_FLD::FldArrayF*> unstrFc;
  vector<E_Int> nitc; 
  vector<E_Int> njtc; 
  vector<E_Int> nktc;
  vector<K_FLD::FldArrayI*> cntc;
  vector<char*> eltTypec;
  vector<PyObject*> objsc, objuc;
  skipNoCoord = false;
  skipStructured = false;
  skipUnstructured = true;
  skipDiffVars = true;
  res = K_ARRAY::getFromArrays(
    cellns, resc, structVarStringc, unstrVarStringc,
    structFc, unstrFc, nitc, njtc, nktc, cntc, eltTypec, objsc, objuc, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true); 
  E_Int nsc = structFc.size();
  if ( res == -1 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "setDoublyDefinedBC: 4th argument is not valid."); 
    RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); 

    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objs[is], structF[is]);
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsc[is], structFc[is]);
    return NULL;
  }

  E_Int nzones = nit.size(); E_Int ncellns = nitc.size();
  if ( nzones != ncellns ) 
  {
    RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); 
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objs[is], structF[is]);
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsc[is], structFc[is]);
    PyErr_SetString(PyExc_TypeError, 
                    "setDoublyDefinedBC: 3rd and 4th arguments must be lists of same size.");
    return NULL;
  }
  //check dimensions
  for (E_Int v = 0; v < nzones; v++)
  {
    res = checkDimensions(nit[v], njt[v], nkt[v], nitc[v], njtc[v], nktc[v]);
    if (res == 0) 
    {
      RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); 
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objs[is], structF[is]);
      for (E_Int is = 0; is < nsc; is++)
        RELEASESHAREDS(objsc[is], structFc[is]);
      PyErr_SetString(PyExc_TypeError, 
                      "setDoublyDefinedBC: 3rd and 4th arguments not compliant.");
      return NULL;  
    }
  }
  // Check: range de la BC
  E_Int range[6];
  E_Int sizeRange = PyList_Size(range1);
  for (int i = 0; i <  sizeRange; i++)
  {
    PyObject* tpl = PyList_GetItem(range1, i);
    if (PyLong_Check(tpl) == 0 && PyInt_Check(tpl) == 0)
    {
      RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); 
      for (E_Int is = 0; is < nsc; is++)
        RELEASESHAREDS(objsc[is], structFc[is]);
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objs[is], structF[is]);
      PyErr_SetString(PyExc_TypeError,
                      "setDoublyDefinedBC: range value must be an integer.");
      return NULL;
    }
    range[i] = PyLong_AsLong(tpl);
    if ( range[i] < 1 ) 
    {
      RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc); 
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objs[is], structF[is]);
      for (E_Int is = 0; is < nsc; is++)
        RELEASESHAREDS(objsc[is], structFc[is]);
      PyErr_SetString(PyExc_TypeError,
                      "setDoublyDefinedBC: range value must be positive.");
      return NULL;
    }    
  }
  E_Int posxt = K_ARRAY::isCoordinateXPresent(structVarString[0]);
  E_Int posyt = K_ARRAY::isCoordinateYPresent(structVarString[0]);
  E_Int poszt = K_ARRAY::isCoordinateZPresent(structVarString[0]);
  E_Int posct = K_ARRAY::isCellNatureField2Present(structVarStringc[0]);
  posxt++; posyt++; poszt++; posct++;
  // maillage en noeuds->centres a interpoler
  FldArrayF nodes(im*jm*km,3); FldArrayF centers; 
  nodes.setOneField(*f, posx, 1);
  nodes.setOneField(*f, posy, 2);
  nodes.setOneField(*f, posz, 3);
  FldArrayF* cellnout = new FldArrayF(fc->getSize());
  cellnout->setOneField(*fc, posc, 1);
  res = K_LOC::node2centerStruct(nodes, im, jm, km,  -1, 1, centers);

  E_Int ok = modifyCellNForDoublyDefined(imc, jmc, kmc, depth, range, 
                                         centers.begin(1), centers.begin(2), centers.begin(3), cellnout->begin(),
                                         posxt, posyt, poszt, posct, 
                                         nit, njt, nkt, structF, nitc, njtc, nktc, structFc);
  if ( ok == -1 ) 
  {
    delete cellnout;
    RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc);      
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsc[is], structFc[is]);
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objs[is], structF[is]);
    return NULL;
  }
  // build output array
  PyObject* tpl = K_ARRAY::buildArray(*cellnout, varStringc, imc, jmc, kmc);
  delete cellnout;
  RELEASESHAREDS(a1, f); RELEASESHAREDS(celln1, fc);      
  for (E_Int is = 0; is < nsc; is++)
    RELEASESHAREDS(objsc[is], structFc[is]);
  for (E_Int is = 0; is < ns; is++)
    RELEASESHAREDS(objs[is], structF[is]);
  return tpl;
}
 
//=============================================================================
/* Verifie que la dimension en centres et la dim en noeuds sont coherents */
//=============================================================================
E_Int K_CONNECTOR::checkDimensions(E_Int im, E_Int jm, E_Int km, 
                                   E_Int imc, E_Int jmc, E_Int kmc)
{
  if ( imc == K_FUNC::E_max(im-1,1) && 
       jmc == K_FUNC::E_max(jm-1,1) && 
       kmc == K_FUNC::E_max(km-1,1) ) return 1;
  else return 0;
}
//=============================================================================
/* Modification du celln sur les parois doubly defined lorsque le 
   point n est pas interpolable
   IN : xc, yc, zc : maillage en centres de taille imc,jmc,kmc
*/
//=============================================================================
E_Int K_CONNECTOR::modifyCellNForDoublyDefined(
  E_Int imc, E_Int jmc, E_Int kmc, E_Int depth, E_Int* range, 
  E_Float* xc, E_Float* yc, E_Float* zc, E_Float* celln,
  E_Int posxt, E_Int posyt, E_Int poszt, E_Int posct, 
  vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt, 
  vector<FldArrayF*>& structF, vector<E_Int>& nitc, vector<E_Int>& njtc, 
  vector<E_Int>& nktc, vector<FldArrayF*>& structFc)
{
  // Interpolation type
  //K_INTERP::InterpAdt::InterpolationType interpType=K_INTERP::InterpAdt::O2CF;
  E_Int nindi = 1; E_Int ncf = 8; 
  FldArrayI indi(nindi);// indice de la cellule d'interp
  FldArrayF cf(ncf);// coefs d'interp
  FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);

  E_Int nzones = structF.size();
  E_Int isBuilt;

  //creation des kmeshes et interpDatas
  vector<K_INTERP::InterpAdt*> listOfInterpDatas;
  vector<FldArrayF*> listOfExtCenters;
  vector<E_Int> niet; vector<E_Int> njet; vector<E_Int> nket;
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    E_Int ni = nit[zone]; E_Int nj = njt[zone]; E_Int nk = nkt[zone];
    FldArrayF& field = *structF[zone];
    FldArrayF coord(ni*nj*nk,3); 
    coord.setOneField(field, posxt, 1);
    coord.setOneField(field, posyt, 2);
    coord.setOneField(field, poszt, 3);

    E_Int nie = ni+1; E_Int nje = nj+1; E_Int nke = nk+1;
    FldArrayF* fec = new FldArrayF(nie*nje*nke,3);
    K_LOC::node2ExtCenterStruct(ni, nj, nk,  coord,
                                nie, nje, nke, *fec);
    listOfExtCenters.push_back(fec);
    niet.push_back(nie); njet.push_back(nje); nket.push_back(nke); 
    K_INTERP::InterpAdt* interpData = new K_INTERP::InterpAdt(fec->getSize(), 
                                                              fec->begin(1),fec->begin(2),fec->begin(3),
                                                              &nie, &nje, &nke, isBuilt);  
    if (isBuilt == 0 ) 
    {
      for (size_t noi = 0; noi < listOfInterpDatas.size(); noi++)
        delete listOfInterpDatas[noi];
      PyErr_SetString(PyExc_TypeError,
                      "setDoublyDefinedBC: 2D structured donor zone must be with z=constant.");
      return -1;
    }
    listOfInterpDatas.push_back(interpData);
  }
  //range de la condition aux limites en centres
  E_Int ic1 = K_FUNC::E_max(1,range[0]-1)-1;
  E_Int ic2 = K_FUNC::E_max(1,range[1]-1)-1;
  E_Int jc1 = K_FUNC::E_max(1,range[2]-1)-1;
  E_Int jc2 = K_FUNC::E_max(1,range[3]-1)-1;
  E_Int kc1 = K_FUNC::E_max(1,range[4]-1)-1;
  E_Int kc2 = K_FUNC::E_max(1,range[5]-1)-1;

  E_Int imcjmc = imc*jmc;
  E_Int ind1, ni2, nj2, ni2nj2, ind2, inc;
  E_Float x, y, z, prod;
  short found=0;
  switch (depth) 
  {
    case 1:
      for (E_Int k1 = kc1; k1 <= kc2; k1++)
        for (E_Int j1 = jc1; j1 <= jc2; j1++)
          for (E_Int i1 = ic1; i1 <= ic2; i1++)
          {
            ind1 = i1 + j1 * imc + k1 * imcjmc;
            x = xc[ind1]; y = yc[ind1]; z = zc[ind1];
            found = 0;
            
            for (E_Int v = 0; v < nzones; v++)
            {
              K_INTERP::InterpAdt* interpData = listOfInterpDatas[v];
              ni2 = nit[v]-1; nj2 = njt[v]-1; ni2nj2 = ni2*nj2;
              FldArrayF& celln2 = *structFc[v];
              E_Float* celln2t = celln2.begin(posct); 
              E_Float voli = 0.; E_Int type = 0; E_Int noblk = 0;
              short tmp = K_INTERP::getInterpolationCell(x, y, z, interpData,
                                                          listOfExtCenters[v], &niet[v], &njet[v], &nket[v], NULL,
                                                          1, 2, 3, 0, voli, indi, cf, tmpIndi, tmpCf, type, noblk);
              if ( tmp < 1 ) found = 0;
              else 
              {
                found = 1;
                E_Int indExt = indi[0];
                E_Int extrapB = 0;
                FldArrayI indTab;
                K_LOC::fromExtCenters2StdCenters(niet[v], njet[v], nket[v], indExt, type, indTab, extrapB);
                prod = 1.;
                for (E_Int kc = 2*type; kc < 3*type; kc++)
                  for (E_Int jc = type; jc < 2*type; jc++)
                    for (E_Int ic = 0; ic < type; ic++)
                    {
                      E_Int indloc = indTab[ic] + indTab[jc]* ni2 + indTab[kc]*ni2nj2;
                      prod *= celln2t[indloc];
                    }
//                 if ( prod == 0 || prod == 4096 ) found = 0;
                if ( prod != 1. ) found = 0;
                else 
                {
                  found = 1; 
                  if ( celln[ind1] != 0. ) celln[ind1] = 2.;
                  goto next;
                }//interpolable 
              }
            } //parcours de ts les blocs
            
            // pas en recouvrement avec un autre bloc
            if ( found == 0 )  celln[ind1] = K_FUNC::E_min(1.,celln[ind1]); 
            next:;
          }
      break;
      
    case 2:
      if ( ic1 == ic2 && ic1 == 0 ) inc = 1;
      else if ( jc1 == jc2 && jc1 == 0 ) inc = imc;
      else if ( kc1 == kc2 && kc1 == 0 && kmc != 1) inc = imcjmc;
      else if ( ic1 == ic2 && ic1 > 0 ) inc = -1;
      else if ( jc1 == jc2 && jc1  > 0 ) inc =-imc;
      else if ( kc1 == kc2 && kc2  > 0 ) inc =-imcjmc;
      else 
      {
        PyErr_SetString(PyExc_TypeError,
                        "setDoublyDefinedBC: invalid increment.");
        return -1;
      }

      for (E_Int k1 = kc1; k1 <= kc2; k1++)
        for (E_Int j1 = jc1; j1 <= jc2; j1++)
          for (E_Int i1 = ic1; i1 <= ic2; i1++)
          {
            // 1- point depth = 1 interpolable
            found = 0; ind1 = i1 + j1 * imc + k1 * imcjmc;
            ind2 = ind1 + inc;
            x = xc[ind1]; y = yc[ind1]; z = zc[ind1];
            for (E_Int v = 0; v < nzones; v++)
            {
              K_INTERP::InterpAdt* interpData = listOfInterpDatas[v];
              ni2 = nit[v]-1; nj2 = njt[v]-1; ni2nj2 = ni2*nj2;
              FldArrayF& celln2 = *structFc[v];
              E_Float* celln2t = celln2.begin(posct); 
              E_Float voli = 0.; E_Int type = 0; E_Int noblk = 0;
              void* niev = (void*)&(niet[v]);  
              void* njev = (void*)&(njet[v]);  
              void* nkev = (void*)&(nket[v]);  
              short tmp = K_INTERP::getInterpolationCell(x, y, z, interpData,
                                                         listOfExtCenters[v], niev, njev, nkev, NULL,
                                                          1, 2, 3, 0, voli, indi, cf, tmpIndi, tmpCf, type, noblk);
              if ( tmp < 1 ) found = 0;
              else 
              {
                found = 1;
                E_Int indExt = indi[0];
                E_Int extrapB = 0;
                FldArrayI indTab;
                K_LOC::fromExtCenters2StdCenters(niet[v], njet[v], nket[v], indExt, type, indTab, extrapB);
                prod = 1.;
                for (E_Int kc = 2*type; kc < 3*type; kc++)
                  for (E_Int jc = type; jc < 2*type; jc++)
                    for (E_Int ic = 0; ic < type; ic++)
                    {
                      E_Int indloc = indTab[ic] + indTab[jc]* ni2 + indTab[kc]*ni2nj2;
                      prod *= celln2t[indloc];
                    }
                
                //if ( prod == 0 || prod == 4096 ) found = 0;
                if ( prod != 1. ) found = 0;
                else 
                {
                  found = 1;  
                  if ( celln[ind1] != 0. ) celln[ind1] = 2.;
                  goto depth2;}//interpolable 
              }
            } //parcours de ts les blocs

            // pas en recouvrement avec un autre bloc
            if ( found == 0 ) 
            { 
              celln[ind1] = K_FUNC::E_min(1.,celln[ind1]); 
              celln[ind2] = K_FUNC::E_min(1.,celln[ind2]);
              goto next2;
            }

            //2- point depth = 2 interpolable
            depth2:;
            found = 0; 
            x = xc[ind2]; y = yc[ind2]; z = zc[ind2];
            for (E_Int v = 0; v < nzones; v++)
            {
              K_INTERP::InterpAdt* interpData = listOfInterpDatas[v];
              ni2 = nit[v]-1; nj2 = njt[v]-1; ni2nj2 = ni2*nj2;
              FldArrayF& celln2 = *structFc[v];
              E_Float* celln2t = celln2.begin(posct); 
              E_Float voli = 0.; E_Int type = 0; E_Int noblk = 0;
              short tmp = K_INTERP::getInterpolationCell(x, y, z, interpData,
                                                          listOfExtCenters[v], &niet[v], &njet[v], &nket[v], NULL,
                                                          1, 2, 3, 0, voli, indi, cf, tmpIndi, tmpCf, type, noblk);
              if ( tmp < 1 ) found = 0;
              else 
              {
                found = 1;
                E_Int indExt = indi[0];
                E_Int extrapB = 0;
                FldArrayI indTab;
                K_LOC::fromExtCenters2StdCenters(niet[v], njet[v], nket[v], indExt, type, indTab, extrapB);

                prod = 1.;
                for (E_Int kc = 2*type; kc < 3*type; kc++)
                  for (E_Int jc = type; jc < 2*type; jc++)
                    for (E_Int ic = 0; ic < type; ic++)
                    {
                      E_Int indloc = indTab[ic] + indTab[jc]* ni2 + indTab[kc]*ni2nj2;
                      prod *= celln2t[indloc];
                    } 
                //if ( prod == 0 || prod == 4096 ) found = 0;
                if ( prod != 1. ) found = 0;
                else 
                {
                  found = 1; 
                  if (celln[ind2] != 0.) celln[ind2] = 2.; 
                  goto next2;
                }//depth = 2 interpolable 
              }
            } //parcours de ts les blocs
            if ( found == 0 ) 
            {
              celln[ind1] = K_FUNC::E_min(1.,celln[ind1]); 
              celln[ind2] = K_FUNC::E_min(1.,celln[ind2]);
            }
            next2:;
          }      
      break;
  }
  //delete kmeshes, extkmeshes, interpdatas
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    delete listOfExtCenters[zone];
    delete listOfInterpDatas[zone];  
  }
  return 1;
}

