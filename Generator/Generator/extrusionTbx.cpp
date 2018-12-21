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

// ============================================================================
/* calcul du facteur d'amplification pour la hauteur h a extruder 
   normales contient les coordonnees des centres + les normales aux centres
*/
// ============================================================================
PyObject* K_GENERATOR::getLocalStepFactor(PyObject* self, PyObject* args)
{
  E_Float tolps = 0.1;
  PyObject *array, *normales; 
  if (!PyArg_ParseTuple(args, "OO", &array, &normales)) return NULL;

  // Check arrays : surface + normales sx,sy,sz
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, 
                                    im, jm, km, cn, eltType, true);
  if (res != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLocalStepFactor: array must unstructured.");
    if (res == 1) RELEASESHAREDS(array,f); return NULL;
  }
  E_Int posx0 = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy0 = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz0 = K_ARRAY::isCoordinateZPresent(varString);
  if (posx0 == -1 || posy0 == -1 || posz0 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLocalStepFactor: array must contain coordinates.");
    RELEASESHAREDU(array,f,cn); return NULL;
  }
  posx0++; posy0++; posz0++;
  E_Float* xn = f->begin(posx0);
  E_Float* yn = f->begin(posy0);
  E_Float* zn = f->begin(posz0);

  FldArrayF* fn;
  char* varStringn;
  char* eltTypen;
  FldArrayI* cnn;
  res = K_ARRAY::getFromArray(normales, varStringn, fn, 
                              im, jm, km, cnn, eltTypen, true);
  if (res != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLocalStepFactor: array must unstructured.");
    RELEASESHAREDU(array, f, cn); 
    if (res == 1) RELEASESHAREDS(normales, fn);
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varStringn);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varStringn);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varStringn);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLocalStepFactor: array must contain coordinates.");
    RELEASESHAREDU(array,f,cn); RELEASESHAREDU(normales,fn,cnn); 
    return NULL;
  }
  posx++; posy++; posz++;
  E_Int possx = K_ARRAY::isNamePresent("sx",varStringn);
  E_Int possy = K_ARRAY::isNamePresent("sy",varStringn);
  E_Int possz = K_ARRAY::isNamePresent("sz",varStringn);
  if (possx == -1 || possy == -1 || possz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLocalStepFactor: array must contain sx,sy and sz variables.");
    RELEASESHAREDU(array,f,cn);  RELEASESHAREDU(normales,fn,cnn); 
    return NULL;
  }
  possx++; possy++; possz++;

  E_Int npts = f->getSize();
  FldArrayF fout(npts,1);
  vector< vector<E_Int> > cVE(npts);
  K_CONNECT::connectEV2VE(*cn, cVE);
  E_Float* xc = fn->begin(posx);
  E_Float* yc = fn->begin(posy);
  E_Float* zc = fn->begin(posz);
  E_Float* sx = fn->begin(possx);
  E_Float* sy = fn->begin(possy);
  E_Float* sz = fn->begin(possz);
  fout.setAllValuesAt(1.);
  E_Int nvert = cn->getNfld();
  
#pragma omp parallel default(shared)
  {
    E_Int indP1, indP2;
    E_Float snc0, snc1, snc2, xpm0, xpm1, xpm2;

#pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      vector<E_Int>& voisins = cVE[ind];
      E_Int nvoisins = voisins.size();
      E_Float psmax = tolps;
      for (E_Int noe1 = 0; noe1 < nvoisins; noe1++)
      {
        E_Int inde1 = voisins[noe1];
        E_Float sx1 = sx[inde1];
        E_Float sy1 = sy[inde1];
        E_Float sz1 = sz[inde1];
        E_Float xc1 = xc[inde1];// M1
        E_Float yc1 = yc[inde1];
        E_Float zc1 = zc[inde1];

        for (E_Int noe2 = noe1+1; noe2 < nvoisins; noe2++)        
        {
          E_Int inde2 = voisins[noe2];
          E_Float sx2 = sx[inde2];
          E_Float sy2 = sy[inde2];
          E_Float sz2 = sz[inde2];
          E_Float xc2 = xc[inde2];//M2
          E_Float yc2 = yc[inde2];
          E_Float zc2 = zc[inde2];
          E_Float ps = sx1*sx2+sy1*sy2+sz1*sz2;       
          if (ps < psmax) 
          {
            psmax = ps;
            //recherche de l'arete commune
            indP1 = -1; indP2 = -1;
            for (E_Int nov1 = 1; nov1 <= nvert; nov1++)
            {
              for (E_Int nov2 = 1; nov2 <= nvert; nov2++)
              {
                if ((*cn)(inde1,nov1) == (*cn)(inde2,nov2)) 
                {
                  if (indP1 == -1) {indP1 =(*cn)(inde1,nov1)-1; goto nextv1;}
                  else if (indP2 == -1) 
                  {
                    indP2 = (*cn)(inde1,nov1)-1; goto fin;
                  }
                }
              }
              nextv1:;
            }
            fin:;
            if (indP1 != -1 && indP2 != -1)
            {
              E_Float xp = 0.5*(xn[indP1]+xn[indP2]); E_Float xm = 0.5*(xc2+xc1);
              E_Float yp = 0.5*(yn[indP1]+yn[indP2]); E_Float ym = 0.5*(yc2+yc1); 
              E_Float zp = 0.5*(zn[indP1]+zn[indP2]); E_Float zm = 0.5*(zc2+zc1);
              snc0 = 0.5*(sx1+sx2); snc1 = 0.5*(sy1+sy2); snc2 = 0.5*(sz1+sz2);
              xpm0 = xm-xp; xpm1 = ym-yp; xpm2 = zm-zp;
              E_Float ps2 = snc0*xpm0+snc1*xpm1+snc2*xpm2;
            
              if (ps2 > -tolps) // concavite
              {
                E_Float alpha = acos(ps)/2.;
                E_Float cas2 = K_FUNC::E_abs(cos(alpha));
                E_Float sas2 = K_FUNC::E_abs(sin(alpha));
                fout[ind] = K_FUNC::E_min(1.5,1./K_FUNC::E_max(cas2,sas2));
              }
            }
          }
        }
      }
    }
  }
  PyObject* tpl = K_ARRAY::buildArray(fout, "ht", *cn, -1, eltType);
  RELEASESHAREDU(array, f, cn);
  RELEASESHAREDU(normales, fn, cnn);
  return tpl;
}
