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

// TFI generator - TRI

# include "generator.h"

using namespace std;
using namespace K_FLD;

//===========================================================================
/* TFI TRI */
//===========================================================================
PyObject* K_GENERATOR::TFITRI(PyObject* arrays)
{
  // Extract infos from arrays
  vector<E_Int> res;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> fields;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; 
  vector<E_Int> njt; 
  vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true; 
  E_Boolean skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    fields, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured);

  E_Int nzones = fields.size();
  E_Int nfld = 0;
  if (nzones != 0) nfld = fields[0]->getNfld();
  if (isOk == -1 || nfld < 3 || nzones != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "TFI: invalid list of arrays.");
    for (E_Int v = 0; v < nzones; v++) delete fields[v];
    return NULL;
  }
  // verification que les arrays sont bien 1D
  E_Int ni0 = nit[0];
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int ni = nit[v]; E_Int nj = njt[v]; E_Int nk = nkt[v];
    if ( ni < 2 || nj != 1 || nk != 1 )
    {
       for (E_Int v2 = 0; v2 < nzones; v2++)
         delete fields[v2];
      PyErr_SetString(PyExc_TypeError,
                      "TFI: one array is invalid: must be i-varying only.");
      return NULL;
    }
    if (ni != ni0) 
    {
      for (E_Int v2 = 0; v2 < nzones; v2++) delete fields[v2];
      PyErr_SetString(PyExc_TypeError,
                      "TFI: all arrays must have the same size.");
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
  FldArrayIS newOrder(3);
  E_Int ok = reorderTFITRI(posx, posy, posz, nit, fields, newOrder);
  if (ok == 0) return 0;
  E_Int imin = newOrder[0];
  E_Int jmin = newOrder[2];
  E_Int diag = newOrder[1];
  FldArrayF* Fimin = fields[imin];
  FldArrayF* Fjmin = fields[jmin];
  FldArrayF* Fdiag = fields[diag];
  
  E_Float* ximin = Fimin->begin(posx);
  E_Float* yimin = Fimin->begin(posy);
  E_Float* zimin = Fimin->begin(posz);

  E_Float* xjmin = Fjmin->begin(posx);
  E_Float* yjmin = Fjmin->begin(posy);
  E_Float* zjmin = Fjmin->begin(posz);

  E_Float* xdiag = Fdiag->begin(posx);
  E_Float* ydiag = Fdiag->begin(posy);
  E_Float* zdiag = Fdiag->begin(posz);
  E_Int ni = Fimin->getSize(); E_Int ni1 = ni-1;

  E_Float ip, jp, ip1; E_Float invni1 = 1./ni1;
  E_Int npts = (ni+1)*ni/2; E_Int nelts = ni1*ni1;

  PyObject* tpl = K_ARRAY::buildArray(3, varString, npts, nelts, -1, "TRI");
  E_Float* coordp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF coord(npts, 3, coordp, true);
  E_Int* cnp = K_ARRAY::getConnectPtr(tpl);
  FldArrayI cn(nelts,3, cnp, true);

  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);
  
  E_Int ind = 0;
  for (E_Int j = 0; j < ni; j++)
    for (E_Int i = 0; i < ni-j; i++)
    {
      ip = i*invni1;
      jp = j*invni1;
      ip1 = 1.-ip-jp;
      
      xt[ind] = ip1 * ( xjmin[i] + ximin[j]  - xjmin[0]) + 
        ip  * ( xdiag[j] + xjmin[i+j]  - xdiag[0] )+ 
        jp  * ( ximin[i+j] + xdiag[ni1-i] - ximin[ni1] );
      yt[ind] = ip1 * ( yjmin[i] + yimin[j]  - yjmin[0]) + 
        ip  * ( ydiag[j] + yjmin[i+j]  - ydiag[0] )+ 
        jp  * ( yimin[i+j] + ydiag[ni1-i] - yimin[ni1] );
      zt[ind] = ip1 * ( zjmin[i] + zimin[j]  - zjmin[0]) + 
        ip  * ( zdiag[j] + zjmin[i+j]  - zdiag[0] )+ 
        jp  * ( zimin[i+j] + zdiag[ni1-i] - zimin[ni1] );
      ind++;
    } 
  

  E_Int et = 0;
  E_Int off1 = 0;
  E_Int off2;
  E_Int* cn1 = cn.begin(1);
  E_Int* cn2 = cn.begin(2);
  E_Int* cn3 = cn.begin(3);

  for (E_Int j = ni1-1; j >= 0; j--)
  {
    off2 = off1 + j + 2; 
    for (E_Int i = 0; i < j; i++)
    {
      cn1[et] = off1 + i + 1;
      cn2[et] = off1 + i + 2;
      cn3[et] = off2 + i + 1;
      et++;
      cn1[et] = off1 + i + 2;
      cn2[et] = off2 + i + 1;
      cn3[et] = off2 + i + 2;
      et++;
    }
    cn1[et] = off1 + j + 1;
    cn2[et] = off1 + j + 2;
    cn3[et] = off2 + j + 1;
    et++;
    off1 = off2;
  }
  delete Fimin; delete Fjmin; delete Fdiag;
  return tpl;
}

//=========================================================================
/* Reordonne les arrays d'entree - cas TRI */
//=========================================================================
E_Int K_GENERATOR::reorderTFITRI(E_Int posx, E_Int posy, E_Int posz,
                                 vector<E_Int>& nit, vector<FldArrayF*>& fields, 
                                 FldArrayIS& newOrder)
{
  E_Float eps = 1.e-6;
  E_Int jm = 1; E_Int km = 1;
  FldArrayIS dejaVu(3); dejaVu.setAllValuesAtNull();

  E_Int ok = 0; E_Int i2 = -1; E_Int j2 = -1; E_Int k2 = -1;
  //1- depart : f1 -> detection raccord en imax
  E_Int vp = 0; E_Int vpp = 0;
  E_Int c = 1;
  while ( c < 4 ) 
  {
    i2 = -1; j2 = -1; k2 = -1;
    for (E_Int vn = 0; vn < 3; vn++)
    {
      if ( dejaVu[vn] == 0 && vp != vn && vpp != vn) 
      {
        i2 = -1; j2 = -1; k2 = -1;        
        E_Int ind1 = 1; 
        if ( c == 1 ) ind1 = nit[vp];
        ok = K_CONNECT::getCoincidentVertexIndices(
          ind1, jm, km, nit[vp], jm, km, nit[vn], jm, km,
          posx, posy, posz, posx, posy, posz,
          *fields[vp], *fields[vn], i2, j2, k2, eps);
        if ( ok > 0 )
        { 
          newOrder[c-1] = vp;
          dejaVu[vn] = 1; vpp = vp; vp = vn;
          goto next;
        }

        if ( c != 1 ) 
        {
          ok = K_CONNECT::getCoincidentVertexIndices(
            1, jm, km, nit[vp], jm, km, nit[vn], jm, km,
            posx, posy, posz, posx, posy, posz,
            *fields[vp], *fields[vn], i2, j2, k2, eps);

          if ( ok > 0 )
          { 
            newOrder[c-1] = vp;
            dejaVu[vn] = 1; vpp = vp; vp = vn; 
            goto next;
          }
        }
      }
    }
    printf("Error: TFI: some arrays are not continuous."); return 0;
    next:;
    if ( ( c == 1 && i2 != nit[vp]) || 
         ( c == 2 && i2 != nit[vp]) )
    {
      K_CONNECT::reorderStructField(nit[vp], jm, km, *fields[vp],-1,2,3);
    }
    c++;
  }
  return 1;
}

