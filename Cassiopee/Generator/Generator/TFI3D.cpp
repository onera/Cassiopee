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

// TFI generator - 3D struct

# include "generator.h"

using namespace std;
using namespace K_FLD;
//===========================================================================
/* TFI 3D: structure 3D */
//===========================================================================
PyObject* K_GENERATOR::TFI3D(PyObject* arrays) 
{
  // Extract infos from arrays
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> fields; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true; 
  E_Boolean skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    fields, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  
  E_Int nzones = fields.size();
  E_Int nfld = 0;
  if (nzones != 0) nfld = fields[0]->getNfld();
  if (isOk == -1 || nfld < 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "TFI: invalid list of arrays.");
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objs[nos], fields[nos]);
    return NULL;
  }
  // verification que les arrays sont bien 2D
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int ni = nit[v]; E_Int nj = njt[v]; E_Int nk = nkt[v];
    if (ni < 2 || nj <2  || nk != 1)
    {
      for (E_Int nos = 0; nos < nzones; nos++)
        RELEASESHAREDS(objs[nos], fields[nos]);
      PyErr_SetString(PyExc_TypeError,
                      "TFI: one array is not valid: must be i and j varying only.");
      return NULL;
    }
  }
  // coordonnees deja verifiees dans getFromArrays
  char* varString = structVarString[0];
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); 
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); 
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;

  // reordonner les arrays le cas echeant
  FldArrayIS newOrder(6);
  E_Int ok = reorderTFI3D(posx, posy, posz, nit, njt, fields, newOrder);
  if (ok == 0) 
  {
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objs[nos], fields[nos]);
    PyErr_SetString(PyExc_TypeError,
                    "TFI: input arrays must be C0.");
    return NULL;
  }
  if ( ok == -1)
  {
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objs[nos], fields[nos]);
    PyErr_SetString(PyExc_TypeError,
                    "TFI: error: two arrays might have the same coordinates.");
    return NULL;
  }
  E_Int imin = newOrder[0]; E_Int imax = newOrder[5];
  E_Int jmin = newOrder[1]; E_Int jmax = newOrder[2];
  E_Int kmin = newOrder[3]; E_Int kmax = newOrder[4];
  E_Int ni = nit[jmin]; E_Int nj = nit[imin]; E_Int nk = njt[jmin]; 
  E_Int npts = ni*nj*nk;
  PyObject* tpl = K_ARRAY::buildArray(nfld, varString, ni, nj, nk);
  E_Float* coordp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF coord(npts, nfld, coordp, true);

  //short isok = TFIstruct3D(ni, nj, nk, nfld, imin, imax, jmin, jmax, kmin, kmax,
  //                         fields, coord);
  short isok = TFIstruct3D2(ni, nj, nk, nfld, posx, posy, posz, 
                            imin, imax, jmin, jmax, kmin, kmax,
                            fields, coord);
  
  if (isok == 0) //echec
  {
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objs[nos], fields[nos]);
    if ( isok == -1)
      PyErr_SetString(PyExc_TypeError,
                      "TFI: imin and imax borders are not of same size ni.");
    else if ( isok == -2)
      PyErr_SetString(PyExc_TypeError,
                      "TFI: jmin and jmax borders are not of same size nj.");
    else if ( isok == -3)
      PyErr_SetString(PyExc_TypeError,
                      "TFI: kmin and kmax borders are not of same size nk.");    
    else    
      PyErr_SetString(PyExc_TypeError,
                      "TFI: input arrays are not valid.");
    return NULL;
  }
  for (E_Int nos = 0; nos < nzones; nos++)
    RELEASESHAREDS(objs[nos], fields[nos]);
  return tpl;
}

//============================================================================
/* TFI 3D structure.*/
//============================================================================
short K_GENERATOR::TFIstruct3D(
  E_Int ni, E_Int nj, E_Int nk, E_Int nfld, 
  E_Int imin, E_Int imax, E_Int jmin, E_Int jmax, E_Int kmin, E_Int kmax, 
  std::vector<FldArrayF*>& fields, FldArrayF& coords)
{ 
  E_Int ni1 = ni-1;
  E_Int nj1 = nj-1;
  E_Int nk1 = nk-1;
  E_Float invni1 = 1./ni1;
  E_Float invnj1 = 1./nj1;
  E_Float invnk1 = 1./nk1;
  E_Int ind, indij, indjk, indik;
  E_Float ip, jp, kp, ip1, jp1, kp1, t1x, t2x, t3x;
  E_Int ninj = ni*nj;
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* ximin = fields[imin]->begin(eq);
    E_Float* ximax = fields[imax]->begin(eq);
    E_Float* xjmin = fields[jmin]->begin(eq);
    E_Float* xjmax = fields[jmax]->begin(eq);
    E_Float* xkmin = fields[kmin]->begin(eq);
    E_Float* xkmax = fields[kmax]->begin(eq);
    E_Float* xt = coords.begin(eq);
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j * ni + k * ninj;
          ip = i*invni1;ip1 = 1.-ip;
          jp = j*invnj1;jp1 = 1.-jp;
          kp = k*invnk1; kp1 = 1.-kp;
          indij = i+j*ni;
          indjk = j+k*nj;
          indik = i+k*ni;
          t1x = 
            ip1 * ximin[indjk] + ip * ximax[indjk] +
            jp1 * xjmin[indik] + jp * xjmax[indik] +
            kp1 * xkmin[indij] + kp * xkmax[indij];

         t2x = 
           ip1 * (jp1*ximin[k*nj] + jp*ximin[nj1+k*nj]) + ip * (jp1*ximax[k*nj] + jp*ximax[nj1+k*nj]) +
           jp1 * (kp1*xjmin[i]    + kp*xjmin[i+nk1*ni]) + jp * (kp1*xjmax[i]    + kp*xjmax[i+nk1*ni]) +
           kp1 * (ip1*xkmin[j*ni] + ip*xkmin[ni1+j*ni]) + kp * (ip1*xkmax[j*ni] + ip*xkmax[ni1+j*ni]);


         t3x = 
           ip1*jp1* ( kp1 * xkmin[0]     + kp * xkmax[0] ) + 
           ip1*jp * ( kp1 * xkmin[nj1*ni]+ kp * xkmax[nj1*ni] ) +
           ip*jp1 * ( kp1 * xkmin[ni1]   + kp * xkmax[ni1]) + 
           ip*jp  * ( kp1 * xkmin[ni1+nj1*ni] + kp * xkmax[ni1+nj1*ni]);
        
         xt[ind] = t1x - t2x + t3x;
        }
  }
  return 1;
}
//============================================================================
/* TFI 3D structure.*/
//============================================================================
short K_GENERATOR::TFIstruct3D2(
  E_Int ni, E_Int nj, E_Int nk, E_Int nfld,
  E_Int posx, E_Int posy, E_Int posz,
  E_Int imin, E_Int imax, E_Int jmin, E_Int jmax, E_Int kmin, E_Int kmax,
  std::vector<FldArrayF*>& fields, FldArrayF& coords)
{
  E_Int ni1 = ni-1;
  E_Int nj1 = nj-1;
  E_Int nk1 = nk-1;
  E_Int ind, indij, indjk, indik;
  E_Float t1x, t3x, t2x12, t2x13, t2x23;
  E_Int ninj = ni*nj;
  if ( fields[imin]->getSize() != fields[imax]->getSize()) return -1;
  if ( fields[jmin]->getSize() != fields[jmax]->getSize()) return -2;
  if ( fields[kmin]->getSize() != fields[kmax]->getSize()) return -3;

  //  x,y,z pour 6 faces
  E_Float* pondximin = fields[imin]->begin(posx);
  E_Float* pondximax = fields[imax]->begin(posx);
  E_Float* pondxjmin = fields[jmin]->begin(posx);
  E_Float* pondxjmax = fields[jmax]->begin(posx);
  E_Float* pondxkmin = fields[kmin]->begin(posx);
  E_Float* pondxkmax = fields[kmax]->begin(posx);

  E_Float* pondyimin = fields[imin]->begin(posy);
  E_Float* pondyimax = fields[imax]->begin(posy);
  E_Float* pondyjmin = fields[jmin]->begin(posy);
  E_Float* pondyjmax = fields[jmax]->begin(posy);
  E_Float* pondykmin = fields[kmin]->begin(posy);
  E_Float* pondykmax = fields[kmax]->begin(posy);

  E_Float* pondzimin = fields[imin]->begin(posz);
  E_Float* pondzimax = fields[imax]->begin(posz);
  E_Float* pondzjmin = fields[jmin]->begin(posz);
  E_Float* pondzjmax = fields[jmax]->begin(posz);
  E_Float* pondzkmin = fields[kmin]->begin(posz);
  E_Float* pondzkmax = fields[kmax]->begin(posz);

  //  permutations min/max pour 3 couple d'indices
  //E_Int indiceiminjmin = 0+0*ni;
  E_Int indicejminkmin = 0+0*nj;
  //E_Int indiceiminkmin = 0+0*ni;

  //E_Int indiceimaxjmin = ni1+0*ni;
  E_Int indicejmaxkmin = nj1+0*nj;
  //E_Int indiceimaxkmin = ni1+0*ni;

  //E_Int indiceiminjmax = 0+nj1*ni;
  E_Int indicejminkmax = 0+nk1*nj;
  //E_Int indiceiminkmax = 0+nk1*ni;

  //E_Int indiceimaxjmax = ni1+nj1*ni;
  E_Int indicejmaxkmax = nj1+nk1*nj;
  //E_Int indiceimaxkmax = ni1+nk1*ni;

  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* ximin = fields[imin]->begin(eq);
    E_Float* ximax = fields[imax]->begin(eq);
    E_Float* xjmin = fields[jmin]->begin(eq);
    E_Float* xjmax = fields[jmax]->begin(eq);
    E_Float* xkmin = fields[kmin]->begin(eq);
    E_Float* xkmax = fields[kmax]->begin(eq);
    E_Float* xt = coords.begin(eq);
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j * ni + k * ninj;
          indij = i+j*ni;
          indjk = j+k*nj;
          indik = i+k*ni;

          //  permutations courant/fixe min/max pour 3 couple d'indices
          E_Int indicei_jmin = i+0*ni;
          E_Int indicej_kmin = j+0*nj;
          E_Int indicei_kmin = i+0*ni;

          E_Int indiceimin_j = 0+j*ni;
          E_Int indicejmin_k = 0+k*nj;
          E_Int indiceimin_k = 0+k*ni;

          E_Int indiceimax_j = ni1+j*ni;
          E_Int indicejmax_k = nj1+k*nj;
          E_Int indiceimax_k = ni1+k*ni;

          E_Int indicei_jmax = i+nj1*ni;
          E_Int indicej_kmax = j+nk1*nj;
          E_Int indicei_kmax = i+nk1*ni;

          E_Int indice_i_j = i+j*ni;
          E_Int indice_j_k = j+k*nj;
          E_Int indice_i_k = i+k*ni;

          //  2 ponderations par faces donc 12 ponderations
          E_Float pondximindisttotalj = pow( pow((pondximin[indicejmax_k] - pondximin[indicejmin_k]) , 2.0) + pow((pondyimin[indicejmax_k] - pondyimin[indicejmin_k]) , 2.0) + pow((pondzimin[indicejmax_k] - pondzimin[indicejmin_k]) , 2.0) , 0.5);
          E_Float pondximindisttotalk = pow( pow((pondximin[indicej_kmax] - pondximin[indicej_kmin]) , 2.0) + pow((pondyimin[indicej_kmax] - pondyimin[indicej_kmin]) , 2.0) + pow((pondzimin[indicej_kmax] - pondzimin[indicej_kmin]) , 2.0) , 0.5);
          E_Float pondximaxdisttotalj = pow( pow((pondximax[indicejmax_k] - pondximax[indicejmin_k]) , 2.0) + pow((pondyimax[indicejmax_k] - pondyimax[indicejmin_k]) , 2.0) + pow((pondzimax[indicejmax_k] - pondzimax[indicejmin_k]) , 2.0) , 0.5);
          E_Float pondximaxdisttotalk = pow( pow((pondximax[indicej_kmax] - pondximax[indicej_kmin]) , 2.0) + pow((pondyimax[indicej_kmax] - pondyimax[indicej_kmin]) , 2.0) + pow((pondzimax[indicej_kmax] - pondzimax[indicej_kmin]) , 2.0) , 0.5);

          E_Float pondxjmindisttotali = pow( pow((pondxjmin[indiceimax_k] - pondxjmin[indiceimin_k]) , 2.0) + pow((pondyjmin[indiceimax_k] - pondyjmin[indiceimin_k]) , 2.0) + pow((pondzjmin[indiceimax_k] - pondzjmin[indiceimin_k]) , 2.0) , 0.5);
          E_Float pondxjmindisttotalk = pow( pow((pondxjmin[indicei_kmax] - pondxjmin[indicei_kmin]) , 2.0) + pow((pondyjmin[indicei_kmax] - pondyjmin[indicei_kmin]) , 2.0) + pow((pondzjmin[indicei_kmax] - pondzjmin[indicei_kmin]) , 2.0) , 0.5);
          E_Float pondxjmaxdisttotali = pow( pow((pondxjmax[indiceimax_k] - pondxjmax[indiceimin_k]) , 2.0) + pow((pondyjmax[indiceimax_k] - pondyjmax[indiceimin_k]) , 2.0) + pow((pondzjmax[indiceimax_k] - pondzjmax[indiceimin_k]) , 2.0) , 0.5);
          E_Float pondxjmaxdisttotalk = pow( pow((pondxjmax[indicei_kmax] - pondxjmax[indicei_kmin]) , 2.0) + pow((pondyjmax[indicei_kmax] - pondyjmax[indicei_kmin]) , 2.0) + pow((pondzjmax[indicei_kmax] - pondzjmax[indicei_kmin]) , 2.0) , 0.5);

          E_Float pondxkmindisttotali = pow( pow((pondxkmin[indiceimax_j] - pondxkmin[indiceimin_j]) , 2.0) + pow((pondykmin[indiceimax_j] - pondykmin[indiceimin_j]) , 2.0) + pow((pondzkmin[indiceimax_j] - pondzkmin[indiceimin_j]) , 2.0) , 0.5);
          E_Float pondxkmindisttotalj = pow( pow((pondxkmin[indicei_jmax] - pondxkmin[indicei_jmin]) , 2.0) + pow((pondykmin[indicei_jmax] - pondykmin[indicei_jmin]) , 2.0) + pow((pondzkmin[indicei_jmax] - pondzkmin[indicei_jmin]) , 2.0) , 0.5);
          E_Float pondxkmaxdisttotali = pow( pow((pondxkmax[indiceimax_j] - pondxkmax[indiceimin_j]) , 2.0) + pow((pondykmax[indiceimax_j] - pondykmax[indiceimin_j]) , 2.0) + pow((pondzkmax[indiceimax_j] - pondzkmax[indiceimin_j]) , 2.0) , 0.5);
          E_Float pondxkmaxdisttotalj = pow( pow((pondxkmax[indicei_jmax] - pondxkmax[indicei_jmin]) , 2.0) + pow((pondykmax[indicei_jmax] - pondykmax[indicei_jmin]) , 2.0) + pow((pondzkmax[indicei_jmax] - pondzkmax[indicei_jmin]) , 2.0) , 0.5);


          E_Float inv_pondximindisttotalj = 1.0 / pondximindisttotalj ;
          E_Float inv_pondximindisttotalk = 1.0 / pondximindisttotalk ;
          E_Float inv_pondximaxdisttotalj = 1.0 / pondximaxdisttotalj ;
          E_Float inv_pondximaxdisttotalk = 1.0 / pondximaxdisttotalk ;
          E_Float inv_pondxjmindisttotali = 1.0 / pondxjmindisttotali ;
          E_Float inv_pondxjmindisttotalk = 1.0 / pondxjmindisttotalk ;
          E_Float inv_pondxjmaxdisttotali = 1.0 / pondxjmaxdisttotali ;
          E_Float inv_pondxjmaxdisttotalk = 1.0 / pondxjmaxdisttotalk ;
          E_Float inv_pondxkmindisttotali = 1.0 / pondxkmindisttotali ;
          E_Float inv_pondxkmindisttotalj = 1.0 / pondxkmindisttotalj ;
          E_Float inv_pondxkmaxdisttotali = 1.0 / pondxkmaxdisttotali ;
          E_Float inv_pondxkmaxdisttotalj = 1.0 / pondxkmaxdisttotalj ;


          //  2 ponderations par faces donc 12 ponderations
          E_Float pondximindistpartialj = pow( pow((pondximin[indice_j_k] - pondximin[indicejmin_k]) , 2.0) + pow((pondyimin[indice_j_k] - pondyimin[indicejmin_k]) , 2.0) + pow((pondzimin[indice_j_k] - pondzimin[indicejmin_k]) , 2.0) , 0.5);
          E_Float pondximindistpartialk = pow( pow((pondximin[indice_j_k] - pondximin[indicej_kmin]) , 2.0) + pow((pondyimin[indice_j_k] - pondyimin[indicej_kmin]) , 2.0) + pow((pondzimin[indice_j_k] - pondzimin[indicej_kmin]) , 2.0) , 0.5);
          E_Float pondximaxdistpartialj = pow( pow((pondximax[indice_j_k] - pondximax[indicejmin_k]) , 2.0) + pow((pondyimax[indice_j_k] - pondyimax[indicejmin_k]) , 2.0) + pow((pondzimax[indice_j_k] - pondzimax[indicejmin_k]) , 2.0) , 0.5);
          E_Float pondximaxdistpartialk = pow( pow((pondximax[indice_j_k] - pondximax[indicej_kmin]) , 2.0) + pow((pondyimax[indice_j_k] - pondyimax[indicej_kmin]) , 2.0) + pow((pondzimax[indice_j_k] - pondzimax[indicej_kmin]) , 2.0) , 0.5);

          E_Float pondxjmindistpartiali = pow( pow((pondxjmin[indice_i_k] - pondxjmin[indiceimin_k]) , 2.0) + pow((pondyjmin[indice_i_k] - pondyjmin[indiceimin_k]) , 2.0) + pow((pondzjmin[indice_i_k] - pondzjmin[indiceimin_k]) , 2.0) , 0.5);
          E_Float pondxjmindistpartialk = pow( pow((pondxjmin[indice_i_k] - pondxjmin[indicei_kmin]) , 2.0) + pow((pondyjmin[indice_i_k] - pondyjmin[indicei_kmin]) , 2.0) + pow((pondzjmin[indice_i_k] - pondzjmin[indicei_kmin]) , 2.0) , 0.5);
          E_Float pondxjmaxdistpartiali = pow( pow((pondxjmax[indice_i_k] - pondxjmax[indiceimin_k]) , 2.0) + pow((pondyjmax[indice_i_k] - pondyjmax[indiceimin_k]) , 2.0) + pow((pondzjmax[indice_i_k] - pondzjmax[indiceimin_k]) , 2.0) , 0.5);
          E_Float pondxjmaxdistpartialk = pow( pow((pondxjmax[indice_i_k] - pondxjmax[indicei_kmin]) , 2.0) + pow((pondyjmax[indice_i_k] - pondyjmax[indicei_kmin]) , 2.0) + pow((pondzjmax[indice_i_k] - pondzjmax[indicei_kmin]) , 2.0) , 0.5);

          E_Float pondxkmindistpartiali = pow( pow((pondxkmin[indice_i_j] - pondxkmin[indiceimin_j]) , 2.0) + pow((pondykmin[indice_i_j] - pondykmin[indiceimin_j]) , 2.0) + pow((pondzkmin[indice_i_j] - pondzkmin[indiceimin_j]) , 2.0) , 0.5);
          E_Float pondxkmindistpartialj = pow( pow((pondxkmin[indice_i_j] - pondxkmin[indicei_jmin]) , 2.0) + pow((pondykmin[indice_i_j] - pondykmin[indicei_jmin]) , 2.0) + pow((pondzkmin[indice_i_j] - pondzkmin[indicei_jmin]) , 2.0) , 0.5);
          E_Float pondxkmaxdistpartiali = pow( pow((pondxkmax[indice_i_j] - pondxkmax[indiceimin_j]) , 2.0) + pow((pondykmax[indice_i_j] - pondykmax[indiceimin_j]) , 2.0) + pow((pondzkmax[indice_i_j] - pondzkmax[indiceimin_j]) , 2.0) , 0.5);
          E_Float pondxkmaxdistpartialj = pow( pow((pondxkmax[indice_i_j] - pondxkmax[indicei_jmin]) , 2.0) + pow((pondykmax[indice_i_j] - pondykmax[indicei_jmin]) , 2.0) + pow((pondzkmax[indice_i_j] - pondzkmax[indicei_jmin]) , 2.0) , 0.5);

          //  2 ratio par couple de face donc 12 ratio
          E_Float pondximindistratioj =  pondximindistpartialj * inv_pondximindisttotalj ;
          E_Float pondximindistratiok =  pondximindistpartialk * inv_pondximindisttotalk ;
          E_Float pondximaxdistratioj =  pondximaxdistpartialj * inv_pondximaxdisttotalj ;
          E_Float pondximaxdistratiok =  pondximaxdistpartialk * inv_pondximaxdisttotalk ;
          E_Float pondxjmindistratioi =  pondxjmindistpartiali * inv_pondxjmindisttotali ;
          E_Float pondxjmindistratiok =  pondxjmindistpartialk * inv_pondxjmindisttotalk ;
          E_Float pondxjmaxdistratioi =  pondxjmaxdistpartiali * inv_pondxjmaxdisttotali ;
          E_Float pondxjmaxdistratiok =  pondxjmaxdistpartialk * inv_pondxjmaxdisttotalk ;
          E_Float pondxkmindistratioi =  pondxkmindistpartiali * inv_pondxkmindisttotali ;
          E_Float pondxkmindistratioj =  pondxkmindistpartialj * inv_pondxkmindisttotalj ;
          E_Float pondxkmaxdistratioi =  pondxkmaxdistpartiali * inv_pondxkmaxdisttotali ;
          E_Float pondxkmaxdistratioj =  pondxkmaxdistpartialj * inv_pondxkmaxdisttotalj ;

          //  1 ratio moyen par couple de face donc 3 ratio moyen

          E_Float pondxdir_i_ratio = (pondxjmindistratioi + pondxjmaxdistratioi +  pondxkmindistratioi +  pondxkmaxdistratioi) / 4.0 ;
          E_Float pondxdir_j_ratio = (pondximindistratioj + pondximaxdistratioj +  pondxkmindistratioj +  pondxkmaxdistratioj) / 4.0 ;
          E_Float pondxdir_k_ratio = (pondximindistratiok + pondximaxdistratiok +  pondxjmindistratiok +  pondxjmaxdistratiok) / 4.0 ;

          E_Float pondxdir_i_ratio1 = 1.0 - pondxdir_i_ratio ;
          E_Float pondxdir_j_ratio1 = 1.0 - pondxdir_j_ratio ;
          E_Float pondxdir_k_ratio1 = 1.0 - pondxdir_k_ratio ;

          // E_Float pondxdir_i_ratio1 = (pondxjmindistratioi + pondxjmaxdistratioi +  pondxkmindistratioi +  pondxkmaxdistratioi) / 4.0 ;
          // E_Float pondxdir_j_ratio1 = (pondximindistratioj + pondximaxdistratioj +  pondxkmindistratioj +  pondxkmaxdistratioj) / 4.0 ;
          // E_Float pondxdir_k_ratio1 = (pondximindistratiok + pondximaxdistratiok +  pondxjmindistratiok +  pondxjmaxdistratiok) / 4.0 ;

          // E_Float pondxdir_i_ratio = 1.0 - pondxdir_i_ratio1 ;
          // E_Float pondxdir_j_ratio = 1.0 - pondxdir_j_ratio1 ;
          // E_Float pondxdir_k_ratio = 1.0 - pondxdir_k_ratio1 ;

          // U + V + W
          t1x =
            pondxdir_i_ratio1 * ximin[indjk] + pondxdir_i_ratio * ximax[indjk] +
            pondxdir_j_ratio1 * xjmin[indik] + pondxdir_j_ratio * xjmax[indik] +
            pondxdir_k_ratio1 * xkmin[indij] + pondxdir_k_ratio * xkmax[indij];

          // UV
          t2x12 =
            pondxdir_i_ratio1 * pondxdir_j_ratio1 * ximin[indicejmin_k] +
            pondxdir_i_ratio1 * pondxdir_j_ratio  * ximin[indicejmax_k] +
            pondxdir_i_ratio  * pondxdir_j_ratio1 * ximax[indicejmin_k] +
            pondxdir_i_ratio  * pondxdir_j_ratio  * ximax[indicejmax_k] ;


          // UW
          t2x13 =
            pondxdir_i_ratio1 * pondxdir_k_ratio1 * ximin[indicej_kmin] +
            pondxdir_i_ratio1 * pondxdir_k_ratio  * ximin[indicej_kmax] +
            pondxdir_i_ratio  * pondxdir_k_ratio1 * ximax[indicej_kmin] +
            pondxdir_i_ratio  * pondxdir_k_ratio  * ximax[indicej_kmax] ;


          // UW
          t2x23 =
            pondxdir_j_ratio1 * pondxdir_k_ratio1 * xjmin[indicei_kmin] +
            pondxdir_j_ratio1 * pondxdir_k_ratio  * xjmin[indicei_kmax] +
            pondxdir_j_ratio  * pondxdir_k_ratio1 * xjmax[indicei_kmin] +
            pondxdir_j_ratio  * pondxdir_k_ratio  * xjmax[indicei_kmax] ;


          // UVW
          t3x =
            pondxdir_i_ratio1 * pondxdir_j_ratio1 * pondxdir_k_ratio1 * ximin[indicejminkmin] + pondxdir_i_ratio1 * pondxdir_j_ratio1 * pondxdir_k_ratio * ximin[indicejminkmax] +
            pondxdir_i_ratio1 * pondxdir_j_ratio  * pondxdir_k_ratio1 * ximin[indicejmaxkmin] + pondxdir_i_ratio1 * pondxdir_j_ratio  * pondxdir_k_ratio * ximin[indicejmaxkmax] +
            pondxdir_i_ratio  * pondxdir_j_ratio1 * pondxdir_k_ratio1 * ximax[indicejminkmin] + pondxdir_i_ratio  * pondxdir_j_ratio1 * pondxdir_k_ratio * ximax[indicejminkmax] +
            pondxdir_i_ratio  * pondxdir_j_ratio  * pondxdir_k_ratio1 * ximax[indicejmaxkmin] + pondxdir_i_ratio  * pondxdir_j_ratio  * pondxdir_k_ratio * ximax[indicejmaxkmax] ;


          // U + V + W - UV -UW - VW + UVW
          xt[ind] = t1x - t2x12 - t2x13 - t2x23 + t3x;

        }
  }
  return 1;
}


//=========================================================================
/* Reordonne les arrays d entree - cas 3D structure 
   newOrder[0] : imin, newOrder[5] : imax
   newOrder[1] : jmin, newOrder[2] : jmax
   newOrder[3] : kmin, newOrder[4] : kmax 
*/
//=========================================================================
E_Int K_GENERATOR::reorderTFI3D(E_Int posx, E_Int posy, E_Int posz,
                                vector<E_Int>& nit, vector<E_Int>& njt,
                                vector<FldArrayF*>& fields, 
                                FldArrayIS& newOrder, E_Float eps)
{
  E_Int nk0 = 1;
  
  newOrder.setAllValuesAt(-1);//fournit le numero de la face
  newOrder[0] = 0;
  E_Int ni0 = nit[0]; E_Int nj0 = njt[0];
  E_Int c = 1;
  E_Int nof1 = 0; E_Int nof2 = 0;
  FldArrayIS nof1t(4); nof1t.setAllValuesAt(1);
  FldArrayIS nof2t(4); nof2t.setAllValuesAt(1);
  for (E_Int v = 1; v < 6; v++)
  {
    E_Int& ni = nit[v]; E_Int& nj = njt[v];
    E_Int ok = 
      K_CONNECT::detectMatchInterface(ni0, nj0, nk0, posx, posy, posz,
                                      ni, nj, nk0, posx, posy, posz,
                                      *fields[0], *fields[v],
                                      nof1, nof2, eps);
    if (ok == 1)
    {nof1t[nof1-1] = v; nof2t[nof1-1] = nof2; c++;}
    else newOrder[5] = v;
  }
  if (c < 4) return 0;

  c = 1;
  for (E_Int nof1 = 1; nof1 < 5; nof1++)
  {
    E_Int v = nof1t[nof1-1];
    newOrder[c] = v;
    c++;
  }
  if ( c < 5 ) return 0;

  E_Float dx, dy, dz;
  
  //-----------------------//
  // imin : nj x nk points //
  //-----------------------//
  E_Int fimin = newOrder[0]; E_Int fimax = newOrder[5];
  E_Int fjmin = newOrder[1]; E_Int fjmax = newOrder[2];
  E_Int fkmin = newOrder[3]; E_Int fkmax = newOrder[4];
  if ( fimin==-1 || fimax==-1 || fjmin==-1 ||fjmax==-1 ||fkmin==-1 || fkmax==-1) 
    return -1;
    
  /*---------------------------------------------------------*/
  /* verification imin(i,1)  = kmin(1,i), reorder kmin sinon */
  /*---------------------------------------------------------*/
  K_CONNECT::detectMatchInterface(
    nit[fimin], njt[fimin], nk0, posx, posy, posz,
    nit[fkmin], njt[fkmin], nk0, posx, posy, posz,              
    *fields[fimin], *fields[fkmin], nof1, nof2, eps);
  E_Int ind1 = 0;
  if (nof1 != 3) 
  {
    printf("Error: TFI3D: imin connected to kmin must be j=1 for imin face.\n"); return 0;
  }
  switch (nof2)
  {
    case 1: 
      break;
    case 2: 
      K_CONNECT::reorderStructField(nit[fkmin], njt[fkmin],nk0,*fields[fkmin],-1,2,3);
      break;
    case 3:
      K_CONNECT::reorderStructField(nit[fkmin], njt[fkmin],nk0,*fields[fkmin], 2,1,3);      
      break;
    case 4:
      K_CONNECT::reorderStructField(nit[fkmin], njt[fkmin],nk0,*fields[fkmin],-2,1,3);
      break;
  }
  //verif ds le bon sens
  K_CONNECT::detectMatchInterface(nit[fimin], njt[fimin], nk0, posx, posy, posz,
				  nit[fkmin], njt[fkmin], nk0, posx, posy, posz,              
				  *fields[fimin], *fields[fkmin], nof1, nof2, eps);
  if ( nof2 == 2 ) K_CONNECT::reorderStructField(nit[fkmin], njt[fkmin],nk0,*fields[fkmin],-1,2,3);

  dx = (*fields[fimin])(0,posx)-(*fields[fkmin])(0,posx);
  dy = (*fields[fimin])(0,posy)-(*fields[fkmin])(0,posy);
  dz = (*fields[fimin])(0,posz)-(*fields[fkmin])(0,posz);
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps )
    K_CONNECT::reorderStructField(nit[fkmin], njt[fkmin],nk0,*fields[fkmin],1,-2,3);

  /*---------------------------------------------------------*/
  /* verification imin(i,nk) = kmax(1,i), reorder kmax sinon */
  /*---------------------------------------------------------*/
  K_CONNECT::detectMatchInterface(nit[fimin], njt[fimin], nk0, posx, posy, posz,
				  nit[fkmax], njt[fkmax], nk0, posx, posy, posz,              
				  *fields[fimin], *fields[fkmax], nof1, nof2, eps);
  if ( nof1 != 4) 
  {
    printf("Error: TFI3D: imin connected to kmax must be j=jmax for imin face.\n"); return 0;
  }
  switch (nof2)
  {
    case 1: 
      break;
    case 2: 
      K_CONNECT::reorderStructField(nit[fkmax], njt[fkmax],nk0,*fields[fkmax],-1,2,3);
      break;
    case 3:
      K_CONNECT::reorderStructField(nit[fkmax], njt[fkmax],nk0,*fields[fkmax], 2,1,3);      
      break;
    case 4:
      K_CONNECT::reorderStructField(nit[fkmax], njt[fkmax],nk0,*fields[fkmax],-2,1,3);
      break;
  }
  //verif ds le bon sens
  K_CONNECT::detectMatchInterface(nit[fimin], njt[fimin], nk0, posx, posy, posz,
				  nit[fkmax], njt[fkmax], nk0, posx, posy, posz,              
				  *fields[fimin], *fields[fkmax], nof1, nof2, eps);
  if ( nof2 == 2 ) K_CONNECT::reorderStructField(nit[fkmax], njt[fkmax],nk0,*fields[fkmax],-1,2,3);
  ind1 = (njt[fimin]-1)*nit[fimin];
  dx = (*fields[fimin])(ind1,posx)-(*fields[fkmax])(0,posx);
  dy = (*fields[fimin])(ind1,posy)-(*fields[fkmax])(0,posy);
  dz = (*fields[fimin])(ind1,posz)-(*fields[fkmax])(0,posz);
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps )
    K_CONNECT::reorderStructField(nit[fkmax], njt[fkmax],nk0,*fields[fkmax],1,-2,3);

  /*---------------------------------------------------------*/
  /* verification imin(1,i)  = jmin(1,i), reorder jmin sinon */
  /*---------------------------------------------------------*/
  K_CONNECT::detectMatchInterface(nit[fimin], njt[fimin], nk0, posx, posy, posz,
				  nit[fjmin], njt[fjmin], nk0, posx, posy, posz,              
				  *fields[fimin], *fields[fjmin], nof1, nof2, eps);
  if ( nof1 != 1) 
  {
    printf("Error: TFI3D: imin connected to jmin must be i=1 for imin face.\n"); return 0;
  }
  switch (nof2)
  {
    case 1: 
      break;
    case 2: 
      K_CONNECT::reorderStructField(nit[fjmin], njt[fjmin],nk0,*fields[fjmin],-1,2,3);
      break;
    case 3:
      K_CONNECT::reorderStructField(nit[fjmin], njt[fjmin],nk0,*fields[fjmin], 2,1,3);      
      break;
    case 4:
      K_CONNECT::reorderStructField(nit[fjmin], njt[fjmin],nk0,*fields[fjmin],-2,1,3);
      break;
  }
  //verif ds le bon sens
  K_CONNECT::detectMatchInterface(nit[fimin], njt[fimin], nk0, posx, posy, posz,
				  nit[fjmin], njt[fjmin], nk0, posx, posy, posz,              
				  *fields[fimin], *fields[fjmin], nof1, nof2, eps);
  if ( nof2 == 2 ) K_CONNECT::reorderStructField(nit[fjmin], njt[fjmin],nk0,*fields[fjmin],-1,2,3);
  dx = (*fields[fimin])(0,posx)-(*fields[fjmin])(0,posx);
  dy = (*fields[fimin])(0,posy)-(*fields[fjmin])(0,posy);
  dz = (*fields[fimin])(0,posz)-(*fields[fjmin])(0,posz);
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps )
    K_CONNECT::reorderStructField(nit[fjmin], njt[fjmin],nk0,*fields[fjmin],1,-2,3);

  /*---------------------------------------------------------*/
  /* verification imin(nj,i) = jmax(1,i), reorder jmax sinon */
  /*---------------------------------------------------------*/
  K_CONNECT::detectMatchInterface(nit[fimin], njt[fimin], nk0, posx, posy, posz,
				  nit[fjmax], njt[fjmax], nk0, posx, posy, posz,              
				  *fields[fimin], *fields[fjmax], nof1, nof2, eps);
  if ( nof1 != 2) 
  {
    printf("Error: TFI3D: imin connected to jmax must be i=imax for imin face.\n"); return 0;
  }
  switch (nof2)
  {
    case 1: 
      break;
    case 2: 
      K_CONNECT::reorderStructField(nit[fjmax], njt[fjmax],nk0,*fields[fjmax],-1,2,3);
      break;
    case 3:
      K_CONNECT::reorderStructField(nit[fjmax], njt[fjmax],nk0,*fields[fjmax], 2,1,3);      
      break;
    case 4:
      K_CONNECT::reorderStructField(nit[fjmax], njt[fjmax],nk0,*fields[fjmax],-2,1,3);
      break;
  }
  //verif ds le bon sens
  K_CONNECT::detectMatchInterface(nit[fimin], njt[fimin], nk0, posx, posy, posz,
				  nit[fjmax], njt[fjmax], nk0, posx, posy, posz,              
				  *fields[fimin], *fields[fjmax], nof1, nof2, eps);
  if ( nof2 == 2 ) K_CONNECT::reorderStructField(nit[fjmax], njt[fjmax],nk0,*fields[fjmax],-1,2,3);
  ind1 = nit[fimin]-1;
  dx = (*fields[fimin])(ind1,posx)-(*fields[fjmax])(0,posx);
  dy = (*fields[fimin])(ind1,posy)-(*fields[fjmax])(0,posy);
  dz = (*fields[fimin])(ind1,posz)-(*fields[fjmax])(0,posz);
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps )
    K_CONNECT::reorderStructField(nit[fjmax], njt[fjmax],nk0,*fields[fjmax],1,-2,3);

  /*---------------------------------------------------------*/
  /* imax */
  /*---------------------------------------------------------*/
  K_CONNECT::detectMatchInterface(nit[fkmin], njt[fkmin], nk0, posx, posy, posz,
				  nit[fimax], njt[fimax], nk0, posx, posy, posz,              
				  *fields[fkmin], *fields[fimax], nof1, nof2, eps);  
  switch (nof2)
  {
    case 1:
      K_CONNECT::reorderStructField(nit[fimax], njt[fimax],nk0,*fields[fimax],2,1,3);
      break;
    case 2:
      K_CONNECT::reorderStructField(nit[fimax], njt[fimax],nk0,*fields[fimax],-2,1,3);
      break;
    case 3:
      break;
    case 4:
      K_CONNECT::reorderStructField(nit[fimax], njt[fimax],nk0,*fields[fimax],1,-2,3);
      break;
  }
  K_CONNECT::detectMatchInterface(nit[fkmax], njt[fkmax], nk0, posx, posy, posz,
				  nit[fimax], njt[fimax], nk0, posx, posy, posz,              
				  *fields[fkmax], *fields[fimax], nof1, nof2, eps);
  //verif
  K_CONNECT::detectMatchInterface(nit[fkmin], njt[fkmin], nk0, posx, posy, posz,
				  nit[fimax], njt[fimax], nk0, posx, posy, posz,              
				  *fields[fkmin], *fields[fimax], nof1, nof2, eps);
  if ( nof2 == 4 ) K_CONNECT::reorderStructField(nit[fimax], njt[fimax],nk0,*fields[fjmax],1,-2,3);
  // verif des pts
  ind1 = nit[fkmin]-1;
  dx = (*fields[fimax])(0,posx)-(*fields[fkmin])(ind1,posx);
  dy = (*fields[fimax])(0,posy)-(*fields[fkmin])(ind1,posy);
  dz = (*fields[fimax])(0,posz)-(*fields[fkmin])(ind1,posz);
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps )
    K_CONNECT::reorderStructField(nit[fimax], njt[fimax],nk0,*fields[fimax],1,-2,3);

  return 1;
}
