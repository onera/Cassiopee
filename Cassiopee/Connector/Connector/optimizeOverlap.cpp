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
# include "Interp/Interp.h"
# include "connector.h"
using namespace K_FLD;
using namespace std;
#include "stub.h"

//============================================================================
/* Optimise le recouvrement */
//============================================================================
PyObject* K_CONNECTOR::optimizeOverlap(PyObject* self, PyObject* args)
{
  PyObject *coordArray1, *centerArray1;//info de z1
  PyObject *coordArray2, *centerArray2;//info2 de z2
  PyObject *hook1, *hook2;// hook on adt for z1 and z2
  E_Int priorite1;
  E_Int priorite2;
  E_Int isDWO;
  if (!PYPARSETUPLE_(args, OOOO_ III_ OO_,
                    &coordArray1, &centerArray1, &coordArray2, &centerArray2,
                    &priorite1, &priorite2, &isDWO, &hook1, &hook2))
  {
      return NULL;
  }
  E_Int prio1 = 0; E_Int prio2 = 0;
  if (priorite1 < 0)
    printf("Warning: optimizeOverlap: 1st priority must be positive. Set to 0.\n");
  else prio1 = E_Int(priorite1);
  if (priorite2 < 0)
    printf("Warning: optimizeOverlap: 2nd priority must be positive. Set to 0.\n");
  else prio2 = E_Int(priorite2);

  // Check: coordonnees en noeuds de z1
  E_Int im1, jm1, km1;
  FldArrayF* f1;
  FldArrayI* cn1;
  char* varString1;
  char* eltType;
  E_Int res = K_ARRAY::getFromArray3(coordArray1, varString1, f1,
                                     im1, jm1, km1, cn1, eltType);
  if (res != 1)
  {
    RELEASESHAREDB(res, coordArray1, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 1st arg must be structured.");
    return NULL;
  }
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 1st arg must contain coordinates.");
    RELEASESHAREDS(coordArray1, f1); return NULL;
  }
  posx1++; posy1++; posz1++;

  if ( im1 < 2 || jm1 < 2 || km1 < 2 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 1st arg must be a 3D zone.");
    RELEASESHAREDS(coordArray1, f1); return NULL;
  }
  // Check: coordonnees en centres de z1
  E_Int imc1, jmc1, kmc1;
  FldArrayF* fc1; FldArrayI* cnc1;
  char* varStringc1;
  res = K_ARRAY::getFromArray3(centerArray1, varStringc1,fc1,
                               imc1, jmc1, kmc1, cnc1, eltType);
  if (res != 1)
  {
    RELEASESHAREDB(res, centerArray1, fc1, cnc1);
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 2nd arg must be structured.");
    return NULL;
  }
  E_Int posxc1 = K_ARRAY::isCoordinateXPresent(varStringc1);
  E_Int posyc1 = K_ARRAY::isCoordinateYPresent(varStringc1);
  E_Int poszc1 = K_ARRAY::isCoordinateZPresent(varStringc1);
  E_Int poscellN1 = K_ARRAY::isCellNatureField2Present(varStringc1);
  E_Int posvol1 = K_ARRAY::isNamePresent("vol",varStringc1);

  if (posxc1 == -1 || posyc1 == -1 || poszc1 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 2nd arg must contain coordinates.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1); return NULL;
  }
  if (poscellN1 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 2nd arg must contain cellN variable.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1); return NULL;
  }
  if (posvol1 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 2nd arg must contain vol variable.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1); return NULL;
  }
  posxc1++; posyc1++; poszc1++; poscellN1++; posvol1++;

  // Check: coordonnees en noeuds de z2
  E_Int im2, jm2, km2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  res = K_ARRAY::getFromArray3(coordArray2, varString2, f2,
                               im2, jm2, km2, cn2, eltType2);
  if (res != 1)
  {
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDB(res, coordArray2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 3rd arg must be structured.");
    return NULL;
  }
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 3rd arg must contain coordinates.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDS(coordArray2, f2);
    return NULL;
  }
  posx2++; posy2++; posz2++;

  if ( im2 < 2 || jm2 < 2 || km2 < 2 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 3rd arg must be a 3D zone.");
    RELEASESHAREDS(coordArray1, f1);
    RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDS(coordArray2, f2);
    return NULL;
  }

  // Check: coordonnees en centres de z2
  E_Int imc2, jmc2, kmc2;
  FldArrayF* fc2; FldArrayI* cnc2;
  char* varStringc2;
  res = K_ARRAY::getFromArray3(centerArray2, varStringc2, fc2,
                               imc2, jmc2, kmc2, cnc2, eltType);
  if (res != 1)
  {
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDS(coordArray2, f2); RELEASESHAREDB(res, centerArray2, fc2, cnc2);
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 4th arg must be structured.");
    return NULL;
  }
  E_Int posxc2 = K_ARRAY::isCoordinateXPresent(varStringc2);
  E_Int posyc2 = K_ARRAY::isCoordinateYPresent(varStringc2);
  E_Int poszc2 = K_ARRAY::isCoordinateZPresent(varStringc2);
  E_Int poscellN2 = K_ARRAY::isCellNatureField2Present(varStringc2);
  E_Int posvol2 = K_ARRAY::isNamePresent("vol",varStringc2);
  if (posxc2 == -1 || posyc2 == -1 || poszc2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 4th arg must contain coordinates.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDS(coordArray2, f2); RELEASESHAREDS(centerArray2, fc2);
    return NULL;
  }
  if (poscellN2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 4th arg must contain cellN variable.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDS(coordArray2, f2); RELEASESHAREDS(centerArray2, fc2);
    return NULL;
  }
  if (posvol2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: 4th arg must contain vol variable.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDS(coordArray2, f2); RELEASESHAREDS(centerArray2, fc2);
    return NULL;
  }
  posxc2++; posyc2++; poszc2++; poscellN2++; posvol2++;


  E_Int isDW = E_Int(isDWO);
  // Recuperation des ADT
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  void** packet1 = (void**) PyCObject_AsVoidPtr(hook1);
  void** packet2 = (void**) PyCObject_AsVoidPtr(hook2);
#else
  void** packet1 = (void**) PyCapsule_GetPointer(hook1, NULL);
  void** packet2 = (void**) PyCapsule_GetPointer(hook2, NULL);
#endif
  E_Int* type1 = (E_Int*)packet1[0];
  E_Int* type2 = (E_Int*)packet2[0];
  if ( type1[0] != 1 || type2[0] != 1 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: hooks must define an ADT.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDS(coordArray2, f2); RELEASESHAREDS(centerArray2, fc2);
    return NULL;
  }
  if ( type1[1] != 1 || type2[1] != 1 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "optimizeOverlap: only one ADT must be defined per hook.");
    RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(centerArray1, fc1);
    RELEASESHAREDS(coordArray2, f2); RELEASESHAREDS(centerArray2, fc2);
    return NULL;
  }

  K_INTERP::InterpAdt* interpData1 = (K_INTERP::InterpAdt*)(packet1[1]);
  K_INTERP::InterpAdt* interpData2 = (K_INTERP::InterpAdt*)(packet2[1]);

  E_Int api = f1->getApi();
  E_Float* xc1 = fc1->begin(posxc1);
  E_Float* yc1 = fc1->begin(posyc1);
  E_Float* zc1 = fc1->begin(poszc1);
  E_Float* vol1 = fc1->begin(posvol1);
  E_Float* cellN1 = fc1->begin(poscellN1);
  E_Float* xc2 = fc2->begin(posxc2);
  E_Float* yc2 = fc2->begin(posyc2);
  E_Float* zc2 = fc2->begin(poszc2);
  E_Float* vol2 = fc2->begin(posvol2);
  E_Float* cellN2 = fc2->begin(poscellN2);

  //modification du celln pour les zones 1 et 2
  if (prio1 == prio2) // meme priorite pour les 2 zones: compar. volumes
    modifyCellNWithVolCriterion(imc1, jmc1, kmc1, vol1, im1, jm1, km1, f1, interpData1, xc1,yc1,zc1,cellN1,
                                imc2, jmc2, kmc2, vol2, im2, jm2, km2, f2, interpData2, xc2,yc2,zc2,cellN2,
                                isDW);
  else // priorite
    modifyCellNWithPriority(prio1, imc1, jmc1, kmc1, vol1, im1, jm1, km1, f1,interpData1,xc1,yc1,zc1,cellN1,
                            prio2, imc2, jmc2, kmc2, vol2, im2, jm2, km2, f2,interpData2,xc2,yc2,zc2,cellN2,
                            isDW);

  RELEASESHAREDS(coordArray1, f1); RELEASESHAREDS(coordArray2, f2);

  //build Arrays
  PyObject* l = PyList_New(0);
  PyObject* tpl = K_ARRAY::buildArray3(*fc1, varStringc1, imc1, jmc1, kmc1, api);
  RELEASESHAREDS(centerArray1, fc1);
  PyList_Append(l, tpl); Py_DECREF(tpl);
  tpl = K_ARRAY::buildArray3(*fc2, varStringc2, imc2, jmc2, kmc2, api);
  RELEASESHAREDS(centerArray2, fc2);
  PyList_Append(l, tpl); Py_DECREF(tpl);

  return l;
}
//=============================================================================
/* Modification du cellN des deux zones base sur le critere de priorite */
//=============================================================================
void K_CONNECTOR::modifyCellNWithPriority(
  E_Int prio1, E_Int ni1, E_Int nj1, E_Int nk1, E_Float* vol1,
  E_Int nie1, E_Int nje1, E_Int nke1, FldArrayF* extCenters1,
  K_INTERP::InterpAdt* interpData1,
  E_Float* xc1, E_Float* yc1, E_Float* zc1, E_Float* celln1,
  E_Int prio2, E_Int ni2, E_Int nj2, E_Int nk2,E_Float* vol2,
  E_Int nie2, E_Int nje2, E_Int nke2, FldArrayF* extCenters2,
  K_INTERP::InterpAdt* interpData2,
  E_Float* xc2, E_Float* yc2, E_Float* zc2, E_Float* celln2,
  E_Int isDW)
{
  E_Int nature = 0;
  E_Int penalty = 0;

  E_Int ni1nj1 = ni1*nj1; E_Int ni2nj2 = ni2*nj2;
  // Interpolation type
  K_INTERP::InterpAdt::InterpolationType interpType = K_INTERP::InterpAdt::O2CF;
  E_Int nindi = 1; E_Int ncf = 8;

  E_Int ncells1 = ni1*nj1*nk1; E_Int ncells2 = ni2*nj2*nk2;
  FldArrayI interpCells1(ncells1,8); FldArrayI interpCells2(ncells2,8);
  interpCells1.setAllValuesAt(-1); interpCells2.setAllValuesAt(-1);
  E_Int* interpCells11 = interpCells1.begin(1);
  E_Int* interpCells12 = interpCells1.begin(2);
  E_Int* interpCells13 = interpCells1.begin(3);
  E_Int* interpCells14 = interpCells1.begin(4);
  E_Int* interpCells15 = interpCells1.begin(5);
  E_Int* interpCells16 = interpCells1.begin(6);
  E_Int* interpCells17 = interpCells1.begin(7);
  E_Int* interpCells18 = interpCells1.begin(8);
  E_Int* interpCells21 = interpCells2.begin(1);
  E_Int* interpCells22 = interpCells2.begin(2);
  E_Int* interpCells23 = interpCells2.begin(3);
  E_Int* interpCells24 = interpCells2.begin(4);
  E_Int* interpCells25 = interpCells2.begin(5);
  E_Int* interpCells26 = interpCells2.begin(6);
  E_Int* interpCells27 = interpCells2.begin(7);
  E_Int* interpCells28 = interpCells2.begin(8);
  FldArrayI tag1(ncells1); tag1.setAllValuesAtNull();
  FldArrayI tag2(ncells2); tag2.setAllValuesAtNull();
  E_Float constraint = 10.;//contrainte sur la somme des coeffs d extrapolation

  /* 1- Recherche des points interpolables du bloc non prioritaire */
#pragma omp parallel default(shared)
  {
    FldArrayI indi(nindi);// indice de la cellule d'interp
    FldArrayF cf(ncf);// coefs d'interp
    FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);
    
    short found;
    E_Float x, y, z;
    E_Float voli; E_Int type; E_Int noblk;
    FldArrayI indic(8);
    
#pragma omp for schedule(dynamic)
    for (E_Int ind1 = 0; ind1 < ncells1; ind1++)
    {
      if (celln1[ind1] != 0.)
      {
        // Recherche de la cellule d'interpolation
        x = xc1[ind1]; y = yc1[ind1]; z = zc1[ind1];
        voli = 0.; type = 0; noblk = 0;
        found = K_INTERP::getInterpolationCell(x, y, z, interpData2, extCenters2, &nie2, &nje2, &nke2, NULL,
                                               1, 2, 3, 0, voli, indi, cf, tmpIndi, tmpCf, type, noblk, interpType, nature, penalty);
        if (found <= 0 && isDW == 1)
        {
          found = K_INTERP::getExtrapolationCell(x, y, z, interpData2, extCenters2, &nie2, &nje2, &nke2, NULL,
                                                 1, 2, 3, 0, voli, indi, cf, type, noblk, interpType, nature, penalty, constraint);
        }
        if (found > 0) //ind1 interpolable de blk2
        {
          E_Int extrapB = 0;
          FldArrayI indTab;
          K_LOC::fromExtCenters2StdCenters(nie2, nje2, nke2, indi[0], type, indTab, extrapB);
          if ( type == 1 ) 
          {
            E_Int indloc = indTab[0] + indTab[1]*ni2+ indTab[2]*ni2nj2; 
            if (celln2[indloc] != 0. && celln2[indloc] != 2.)
            {
              interpCells11[ind1] = indloc;
              if (prio1 > prio2) tag1[ind1] = 1;     
            }
          }
          else 
          {
            E_Int isvalid = 1;
            E_Int count = 0;
            for (E_Int kc = 2*type; kc < 3*type; kc++)
              for (E_Int jc = type; jc < 2*type; jc++)
                for (E_Int ic = 0; ic < type; ic++)
                {
                  E_Int indloc = indTab[ic] + indTab[jc]* ni2 + indTab[kc]*ni2nj2;             
                  if (celln2[indloc] == 0. || celln2[indloc] == 2.) {isvalid = 0; break;}
                  indic[count] = indloc;
                  count++;
                }
            if (isvalid == 1)
            {
              interpCells11[ind1] = indic[0]; interpCells12[ind1] = indic[1];
              interpCells13[ind1] = indic[2]; interpCells14[ind1] = indic[3];
              interpCells15[ind1] = indic[4]; interpCells16[ind1] = indic[5];
              interpCells17[ind1] = indic[6]; interpCells18[ind1] = indic[7]; 
              if (prio1 > prio2) tag1[ind1] = 1;          
            }
          }
        }
      }
    }
  }


  // rechercher les pts a tagger ds blk2
#pragma omp parallel default(shared)
  {
    FldArrayI indi(nindi);// indice de la cellule d'interp
    FldArrayF cf(ncf);// coefs d'interp
    FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);
    
    short found;
    E_Float x, y, z;
    E_Float voli; E_Int type; E_Int noblk;
    FldArrayI indic(8);

#pragma omp for schedule(dynamic)
    for (E_Int ind2 = 0; ind2 < ncells2; ind2++)
    {
      if (celln2[ind2] != 0.)
      {
        // Recherche de la cellule d'interpolation
        x = xc2[ind2]; y = yc2[ind2]; z = zc2[ind2];
        voli = 0.; type = 0; noblk = 0;
        found = K_INTERP::getInterpolationCell(x, y, z, interpData1, extCenters1, &nie1, &nje1, &nke1, NULL,
                                               1, 2, 3, 0, voli, indi, cf, tmpIndi, tmpCf, type, noblk, interpType, nature, penalty);
        if (found<=0 && isDW == 1)
        {
          found = K_INTERP::getExtrapolationCell(x, y, z, interpData1, extCenters1, &nie1, &nje1, &nke1, NULL,
                                                 1, 2, 3, 0, voli, indi, cf, type, noblk, interpType, nature, penalty, constraint);
        }
        if (found > 0) //ind2 interpolable de blk1 
        {  
          E_Int extrapB = 0;
          FldArrayI indTab;
          K_LOC::fromExtCenters2StdCenters(nie1, nje1, nke1, indi[0], type, indTab, extrapB);    
          if ( type == 1 ) 
          {
            E_Int indloc = indTab[0] + indTab[1]*ni1+ indTab[2]*ni1nj1; 
            if (celln1[indloc] != 0. && celln1[indloc] != 2.)
            {
              interpCells21[ind2] = indloc;
              if (prio2 > prio1) tag2[ind2] = 1;
            }  
          }
          else 
          {     
            E_Int isvalid = 1;
            E_Int count = 0;
            for (E_Int kc = 2*type; kc < 3*type; kc++)
              for (E_Int jc = type; jc < 2*type; jc++)
                for (E_Int ic = 0; ic < type; ic++)
                {
                  E_Int indloc = indTab[ic] + indTab[jc]* ni1 + indTab[kc]*ni1nj1;
                  if ( celln1[indloc] == 0. || celln1[indloc] == 2 ) {isvalid = 0; break;}
                  indic[count] = indloc;
                  count++;
                }
            if (isvalid == 1)
            {
              interpCells21[ind2] = indic[0]; interpCells22[ind2] = indic[1];
              interpCells23[ind2] = indic[2]; interpCells24[ind2] = indic[3];
              interpCells25[ind2] = indic[4]; interpCells26[ind2] = indic[5];
              interpCells27[ind2] = indic[6]; interpCells28[ind2] = indic[7]; 
              if (prio2 > prio1) tag2[ind2] = 1;
            }
          }
        }
      }
    }
  }


  /* 2-interpolations croisees: celln = 1 si une cellule interpolee dans l'etape 1 est donneuse */
  modifyCellnForDonorCells(tag1, celln1, tag2, celln2, interpCells1);
  modifyCellnForDonorCells(tag2, celln2, tag1, celln1, interpCells2);
#pragma omp parallel default(shared)
  {
#pragma omp for schedule(dynamic)
    for (E_Int ind1 = 0; ind1 < ncells1; ind1++)
    {
      if (tag1[ind1] == 1) celln1[ind1] = 3.;
    }
#pragma omp for schedule(dynamic)
    for (E_Int ind2 = 0; ind2 < ncells2; ind2++)
    {
      if (tag2[ind2] == 1) celln2[ind2] = 3.;
    }
  }
}

//=============================================================================
/* Modification du cellN des deux zones base sur le critere de masquage de la
   cellule de plus grand volume */
//=============================================================================
void K_CONNECTOR::modifyCellNWithVolCriterion(
  E_Int ni1, E_Int nj1, E_Int nk1, E_Float* vol1,
  E_Int nie1, E_Int nje1, E_Int nke1, FldArrayF* extCenters1,
  K_INTERP::InterpAdt* interpData1,
  E_Float* xc1, E_Float* yc1, E_Float* zc1, E_Float* celln1,
  E_Int ni2, E_Int nj2, E_Int nk2, E_Float* vol2,
  E_Int nie2, E_Int nje2, E_Int nke2, FldArrayF* extCenters2,
  K_INTERP::InterpAdt* interpData2,
  E_Float* xc2, E_Float* yc2, E_Float* zc2, E_Float* celln2, E_Int isDW)
{
  FldArrayI interpCells1;
  FldArrayI interpCells2;
  E_Int ncells1 = ni1*nj1*nk1; E_Int ncells2 = ni2*nj2*nk2;
  FldArrayI tag1(ncells1); tag1.setAllValuesAtNull();
  FldArrayI tag2(ncells2); tag2.setAllValuesAtNull();

  /* 1-comparaison des volumes des cellules interpolees/ d interpolation
     cellN=2 pour la cellule la plus grossiere */
  compareInterpCells(ni1, nj1, nk1, nie1, nje1, nke1, extCenters1, interpData1,
                     xc1, yc1, zc1, celln1, vol1, interpCells1, tag1,
                     ni2, nj2, nk2, nie2, nje2, nke2, extCenters2, interpData2,
                     xc2, yc2, zc2, celln2, vol2, interpCells2, tag2, isDW);

  /* 2-interpolations croisees: celln = 1 si une cellule interpolee dans l'etape 1 est donneuse*/
  modifyCellnForDonorCells(tag1, celln1, tag2, celln2, interpCells1);
  modifyCellnForDonorCells(tag2, celln2, tag1, celln1, interpCells2);

#pragma omp parallel default(shared)
  {
#pragma omp for schedule(dynamic)
    for (E_Int ind1 = 0; ind1 < ncells1; ind1++)
    {
      if (tag1[ind1] == 1) celln1[ind1] = 3.;
    }
#pragma omp for schedule(dynamic)
    for (E_Int ind2 = 0; ind2 < ncells2; ind2++)
    {
      if (tag2[ind2] == 1) celln2[ind2] = 3.;
    }
  }
}
//=============================================================================
/* comparaison de la taille des cellules interpolables/d interpolation
    des blk1 et blk2. Interpolation en centres etendus effectuee.
   Retourne les celln modifies
   ni1, nj1, nk1 : en centres */
//=============================================================================
void K_CONNECTOR::compareInterpCells(
  E_Int ni1,E_Int nj1, E_Int nk1,
  E_Int nie1, E_Int nje1, E_Int nke1, FldArrayF* extCenters1,
  K_INTERP::InterpAdt* interpData1,
  E_Float* xc1, E_Float* yc1, E_Float* zc1, E_Float* celln1,
  E_Float* vol1, FldArrayI& interpCells1, FldArrayI& tag1,
  E_Int ni2, E_Int nj2, E_Int nk2,
  E_Int nie2, E_Int nje2, E_Int nke2, FldArrayF* extCenters2,
  K_INTERP::InterpAdt* interpData2,
  E_Float* xc2, E_Float* yc2, E_Float* zc2, E_Float* celln2,
  E_Float* vol2, FldArrayI& interpCells2, FldArrayI& tag2, E_Int isDW)
{
  // Interpolation type
  K_INTERP::InterpAdt::InterpolationType interpType = K_INTERP::InterpAdt::O2CF;
  E_Int nindi = 1; E_Int ncf = 8;

  //calcul des volumes des cellules de blk1 et blk2
  E_Int ncells1 = ni1*nj1*nk1; E_Int ncells2 = ni2*nj2*nk2;
  interpCells1.malloc(ncells1,8); interpCells1.setAllValuesAt(-1);
  interpCells2.malloc(ncells2,8); interpCells2.setAllValuesAt(-1);

  E_Int nature = 0;
  E_Int penalty = 0;
  //E_Int indi1, indi2, indi3, indi4, indi5, indi6, indi7, indi8;
  E_Int ni1nj1 = ni1*nj1;
  E_Int ni2nj2 = ni2*nj2;
  E_Int* interpCells11 = interpCells1.begin(1);
  E_Int* interpCells12 = interpCells1.begin(2);
  E_Int* interpCells13 = interpCells1.begin(3);
  E_Int* interpCells14 = interpCells1.begin(4);
  E_Int* interpCells15 = interpCells1.begin(5);
  E_Int* interpCells16 = interpCells1.begin(6);
  E_Int* interpCells17 = interpCells1.begin(7);
  E_Int* interpCells18 = interpCells1.begin(8);
  E_Int* interpCells21 = interpCells2.begin(1);
  E_Int* interpCells22 = interpCells2.begin(2);
  E_Int* interpCells23 = interpCells2.begin(3);
  E_Int* interpCells24 = interpCells2.begin(4);
  E_Int* interpCells25 = interpCells2.begin(5);
  E_Int* interpCells26 = interpCells2.begin(6);
  E_Int* interpCells27 = interpCells2.begin(7);
  E_Int* interpCells28 = interpCells2.begin(8);
  E_Float constraint = 10.;//contrainte sur la somme des coeffs d extrapolation

  // recherche de cellule d interpolation de ind1 sur blk2
#pragma omp parallel default(shared)
  {
    E_Float x, y, z, volc1, volc2;
    short found;
    FldArrayI indi(nindi);// indice de la cellule d interp
    FldArrayF cf(ncf);// coefs d interp
    FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);

#pragma omp for schedule(dynamic)
    for (E_Int ind1 = 0; ind1 < ncells1; ind1++)
    {
      if (celln1[ind1] != 0.)
      {
        volc1 = K_FUNC::E_abs(vol1[ind1]);

        // Recherche de la cellule d'interpolation
        x = xc1[ind1]; y = yc1[ind1]; z = zc1[ind1];
        E_Float voli = 0.; E_Int type = 0; E_Int noblk = 0;
        found = K_INTERP::getInterpolationCell(x, y, z, interpData2, extCenters2, &nie2, &nje2, &nke2, NULL,
                                               1, 2, 3, 0, voli, indi, cf, tmpIndi, tmpCf, type, noblk);
        if (found <= 0 && isDW == 1)
        {
          found = K_INTERP::getExtrapolationCell(x, y, z, interpData2, extCenters2, &nie2, &nje2, &nke2, NULL,
                                                 1, 2, 3, 0, voli, indi, cf, type, noblk, interpType, nature, penalty, constraint);
        }
        if (found > 0) //ind1 interpolable de blk2
        {
          E_Int extrapB = 0;
          FldArrayI indTab;
          FldArrayI indic(8); indic.setAllValuesAt(-1);
          K_LOC::fromExtCenters2StdCenters(nie2, nje2, nke2, indi[0], type, indTab, extrapB);  
          if ( type == 1 ) 
          {
            E_Int indloc = indTab[0] + indTab[1]*ni2+ indTab[2]*ni2nj2; 
            if (celln2[indloc] != 0. && celln2[indloc] != 2.)
            {
              interpCells11[ind1] = indloc;
              volc2 = vol2[indloc];
              if (volc1 > volc2) tag1[ind1] = 1;     
            } 
          }
          else 
          {
            E_Int isvalid = 1;
            E_Int count = 0;
            for (E_Int kc = 2*type; kc < 3*type; kc++)
              for (E_Int jc = type; jc < 2*type; jc++)
                for (E_Int ic = 0; ic < type; ic++)
                {
                  E_Int indloc = indTab[ic] + indTab[jc]* ni2 + indTab[kc]*ni2nj2;
                  if (celln2[indloc] == 0. || celln2[indloc] == 2) {isvalid = 0; break;}
                  indic[count] = indloc;
                  count++;
                }
          
            if (isvalid == 1)
            {
              volc2 = 0.;
              for (E_Int nocf = 0; nocf < 8; nocf++)
              {
                E_Int indcell = indic[nocf];
                volc2 += cf[nocf]*vol2[indcell];
              }
              volc2 = K_FUNC::E_abs(volc2);
              interpCells11[ind1] = indic[0]; interpCells12[ind1] = indic[1];
              interpCells13[ind1] = indic[2]; interpCells14[ind1] = indic[3];
              interpCells15[ind1] = indic[4]; interpCells16[ind1] = indic[5];
              interpCells17[ind1] = indic[6]; interpCells18[ind1] = indic[7]; 
              if (volc1 > volc2) tag1[ind1] = 1;
            }
          }
        }
      }
    }
  }

  // recherche de cellule d interpolation de ind2 sur blk1
#pragma omp parallel default(shared)
  {
    E_Float x, y, z, volc1, volc2;
    short found;
    FldArrayI indi(nindi);// indice de la cellule d interp
    FldArrayF cf(ncf);// coefs d interp
    FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);

#pragma omp for schedule(dynamic)
    for (E_Int ind2 = 0; ind2 < ncells2; ind2++)
    {
      if (celln2[ind2] != 0.)
      {
        volc2 = K_FUNC::E_abs(vol2[ind2]);

        // Recherche de la cellule d'interpolation
        x = xc2[ind2]; y = yc2[ind2]; z = zc2[ind2];
        E_Float voli = 0.; E_Int type = 0; E_Int noblk = 0;
        found = K_INTERP::getInterpolationCell(x, y, z, interpData1, extCenters1, &nie1, &nje1, &nke1, NULL,
                                               1, 2, 3, 0, voli, indi, cf, tmpIndi, tmpCf, type, noblk);
        if ( found <= 0 && isDW == 1)
        {
          found = K_INTERP::getExtrapolationCell(x, y, z, interpData1, extCenters1, &nie1, &nje1, &nke1, NULL,
                                                 1, 2, 3, 0, voli, indi, cf, type, noblk, interpType, nature, penalty, constraint);
        }
        if ( found > 0 )//ind2 interpolable de blk1
        {
          E_Int extrapB = 0;
          FldArrayI indTab;
          FldArrayI indic(8);
          K_LOC::fromExtCenters2StdCenters(nie1, nje1, nke1, indi[0], type, indTab, extrapB);
          if ( type == 1 ) 
          {
            E_Int indloc = indTab[0] + indTab[1]*ni1+ indTab[2]*ni1nj1; 
            if (celln1[indloc] != 0. && celln1[indloc] != 2.)
            {
              interpCells21[ind2] = indloc;
              volc1 = vol1[indloc];
              if (volc2 > volc1) tag2[ind2] = 1;     
            }
          }
          else 
          {
            E_Int isvalid = 1;
            E_Int count = 0;
            for (E_Int kc = 2*type; kc < 3*type; kc++)
              for (E_Int jc = type; jc < 2*type; jc++)
                for (E_Int ic = 0; ic < type; ic++)
                {
                  E_Int indloc = indTab[ic] + indTab[jc]* ni1 + indTab[kc]*ni1nj1;
                  if (celln1[indloc] == 0. || celln1[indloc] == 2) {isvalid = 0; break;}
                  indic[count] = indloc;
                  count++;
                }
            if (isvalid == 1)
            {
              volc1 = 0.;
              for (E_Int nocf = 0; nocf < 8; nocf++)
              {
                E_Int indcell = indic[nocf];
                volc1 += cf[nocf]*vol1[indcell];
              }
              volc1 = K_FUNC::E_abs(volc1);
              interpCells21[ind2] = indic[0]; interpCells22[ind2] = indic[1];
              interpCells23[ind2] = indic[2]; interpCells24[ind2] = indic[3];
              interpCells25[ind2] = indic[4]; interpCells26[ind2] = indic[5];
              interpCells27[ind2] = indic[6]; interpCells28[ind2] = indic[7];  
              if (volc2 > volc1) tag2[ind2] = 1;
            }
          }
        }
      }
    }
  }
}
//=============================================================================
/* Interpolations croisees: celln = 3 si une cellule interpolee du blk1 dans
   l'etape 1 est donneuse du blk2. */
//=============================================================================
void K_CONNECTOR::modifyCellnForDonorCells(FldArrayI& tag1, E_Float* celln1,
                                           FldArrayI& tag2, E_Float* celln2,
                                           FldArrayI& interpCells1)
{
  E_Int* interpCells1p = interpCells1.begin(1);
  E_Int* tag1p = tag1.begin();
  E_Int* tag2p = tag2.begin();

#pragma omp parallel default(shared)
  {
#pragma omp for schedule(dynamic)
    for (E_Int ind1 = 0; ind1 < tag1.getSize(); ind1++)
    {
      if (tag1p[ind1] == 1 || celln1[ind1] == 2.)
      {
        if (interpCells1p[ind1] != -1)
        {
          for (E_Int no = 1; no <= 8; no++)
          {
            E_Int indi2 = interpCells1(ind1,no);
            if ( indi2 == -1) goto nextind;
            tag2p[indi2] = 0;
          }
        }
      }
      nextind:;
    }
  }
}
