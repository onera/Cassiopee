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
                    "getLocalStepFactor: array must be unstructured.");
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

  FldArrayF* fn; FldArrayI* cnn;
  char* varStringn; char* eltTypen;
  res = K_ARRAY::getFromArray(normales, varStringn, fn, 
                              im, jm, km, cnn, eltTypen, true);
  if (res != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLocalStepFactor: array must be unstructured.");
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
    E_Int nvoisins, inde1, inde2;
    E_Float psmax, sx1, sy1, sz1, xc1, yc1, zc1;
    E_Float sx2, sy2, sz2, xc2, yc2, zc2, ps, ps2;
    E_Float xp, yp, zp, alpha, cas2, sas2;

    #pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      vector<E_Int>& voisins = cVE[ind];
      nvoisins = voisins.size();
      psmax = tolps;
      for (E_Int noe1 = 0; noe1 < nvoisins; noe1++)
      {
        inde1 = voisins[noe1];
        sx1 = sx[inde1];
        sy1 = sy[inde1];
        sz1 = sz[inde1];
        xc1 = xc[inde1];// M1
        yc1 = yc[inde1];
        zc1 = zc[inde1];

        for (E_Int noe2 = noe1+1; noe2 < nvoisins; noe2++)        
        {
          inde2 = voisins[noe2];
          sx2 = sx[inde2];
          sy2 = sy[inde2];
          sz2 = sz[inde2];
          xc2 = xc[inde2];//M2
          yc2 = yc[inde2];
          zc2 = zc[inde2];
          ps = sx1*sx2+sy1*sy2+sz1*sz2;       
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
              xp = 0.5*(xn[indP1]+xn[indP2]); E_Float xm = 0.5*(xc2+xc1);
              yp = 0.5*(yn[indP1]+yn[indP2]); E_Float ym = 0.5*(yc2+yc1); 
              zp = 0.5*(zn[indP1]+zn[indP2]); E_Float zm = 0.5*(zc2+zc1);
              snc0 = 0.5*(sx1+sx2); snc1 = 0.5*(sy1+sy2); snc2 = 0.5*(sz1+sz2);
              xpm0 = xm-xp; xpm1 = ym-yp; xpm2 = zm-zp;
              ps2 = snc0*xpm0+snc1*xpm1+snc2*xpm2;
            
              if (ps2 > -tolps) // concavite
              {
                alpha = acos(ps)/2.;
                cas2 = K_FUNC::E_abs(cos(alpha));
                sas2 = K_FUNC::E_abs(sin(alpha));
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

// ============================================================================
/* calcul du facteur d'amplification pour la hauteur h a extruder 
   normales contient les coordonnees des centres + les normales aux centres
   sortie: champ 1 (ht) : hauteur d'extrusion
           champ 2 (hl): champ utilise dans le lissage
*/
// ============================================================================
PyObject* K_GENERATOR::getLocalStepFactor2(PyObject* self, PyObject* args)
{
  //E_Float tolps = 0.5;
  PyObject *array, *normales;
  E_Int kappaType;
  E_Float kappaL, kappaP;
  if (!PyArg_ParseTuple(args, "OOldd", &array, &normales, &kappaType, &kappaL, &kappaP)) return NULL;

  // Check arrays : surface + normales sx,sy,sz
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, 
                                    im, jm, km, cn, eltType, true);
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLocalStepFactor: array must be unstructured.");
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

  FldArrayF* fn; FldArrayI* cnn;
  char* varStringn; char* eltTypen;
  res = K_ARRAY::getFromArray(normales, varStringn, fn, 
                              im, jm, km, cnn, eltTypen, true);
  if (res != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLocalStepFactor: array must be unstructured.");
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
    RELEASESHAREDU(array,f,cn); RELEASESHAREDU(normales,fn,cnn); 
    return NULL;
  }
  possx++; possy++; possz++;

  E_Int npts = f->getSize();
  FldArrayF fout(npts,2);
  vector< vector<E_Int> > cVE(npts);
  K_CONNECT::connectEV2VE(*cn, cVE);
  E_Float* xc = fn->begin(posx);
  E_Float* yc = fn->begin(posy);
  E_Float* zc = fn->begin(posz);
  E_Float* sx = fn->begin(possx);
  E_Float* sy = fn->begin(possy);
  E_Float* sz = fn->begin(possz);
  fout.setAllValuesAt(1.);
  E_Float* fouth = fout.begin(1);
  E_Float* foute = fout.begin(2);
  
  E_Int nvert = cn->getNfld();
  
//#pragma omp parallel default(shared)
  {
    E_Int indP1, indP2;
    //E_Float snc0, snc1, snc2, xpm0, xpm1, xpm2;

//#pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      vector<E_Int>& voisins = cVE[ind];
      E_Int nvoisins = voisins.size();
      //E_Float psmax = tolps;
      for (E_Int noe1 = 0; noe1 < nvoisins; noe1++)
      {
        E_Int inde1 = voisins[noe1];
        //E_Float xc1 = xc[inde1];// M1
        //E_Float yc1 = yc[inde1];
        //E_Float zc1 = zc[inde1];

        for (E_Int noe2 = noe1+1; noe2 < nvoisins; noe2++)        
        {
          E_Int inde2 = voisins[noe2];
          //E_Float sx2 = sx[inde2];
          //E_Float sy2 = sy[inde2];
          //E_Float sz2 = sz[inde2];
          //E_Float xc2 = xc[inde2];//M2
          //E_Float yc2 = yc[inde2];
          //E_Float zc2 = zc[inde2];
          
          //E_Float ps = sx1*sx2+sy1*sy2+sz1*sz2;
          //printf("n1: %f %f %f\n",sx1,sy1,sz1);
          //printf("n2: %f %f %f\n",sx2,sy2,sz2);
          //printf("ps: %f\n",ps);

          //E_Float alpha = acos(ps); // entre 0 et pi
          //if (ps > 0) alpha += K_CONST::pi*0.5;
          
          indP1 = -1; indP2 = -1;
          for (E_Int nov1 = 1; nov1 <= nvert; nov1++)
          {
            for (E_Int nov2 = 1; nov2 <= nvert; nov2++)
            {
              if ((*cn)(inde1,nov1) == (*cn)(inde2,nov2)) 
              {
                if (indP1 == -1) { indP1 =(*cn)(inde1,nov1)-1; goto nextv1; }
                else if (indP2 == -1) { indP2 = (*cn)(inde1,nov1)-1; goto fin; }
              }
            }
            nextv1:;
          }
          fin:;
          if (indP1 != -1 && indP2 != -1)
          {
            //printf("ind=%d, edge=%d %d, elt1=%d, elt2=%d\n", ind, indP1, indP2, inde1, inde2);
            
            // Coordonnees de l'arete
            E_Float x1 = xn[indP1];
            E_Float y1 = yn[indP1];
            E_Float z1 = zn[indP1];
            E_Float x2 = xn[indP2];
            E_Float y2 = yn[indP2];
            E_Float z2 = zn[indP2];
            // point milieu de l'arete
            E_Float xm = 0.5*(x1+x2);
            E_Float ym = 0.5*(y1+y2);
            E_Float zm = 0.5*(z1+z2);
            // Centre des elements adajacents
            E_Float xc1 = xc[inde1];
            E_Float yc1 = yc[inde1];
            E_Float zc1 = zc[inde1];
            E_Float xc2 = xc[inde2];
            E_Float yc2 = yc[inde2];
            E_Float zc2 = zc[inde2];
            // normales aux elements adjacents    
            E_Float sx1 = sx[inde1];
            E_Float sy1 = sy[inde1];
            E_Float sz1 = sz[inde1];
            //E_Float sx2 = sx[inde2];
            //E_Float sy2 = sy[inde2];
            //E_Float sz2 = sz[inde2];
            // Calcul de l'angle alpha
            E_Float v1x = xc1-xm;
            E_Float v1y = yc1-ym;
            E_Float v1z = zc1-zm;
            E_Float v2x = xc2-xm;
            E_Float v2y = yc2-ym;
            E_Float v2z = zc2-zm;
            E_Float norm = 1./K_FUNC::E_max(sqrt(v1x*v1x+v1y*v1y+v1z*v1z), 1.e-10);
            v1x = v1x * norm;
            v1y = v1y * norm;
            v1z = v1z * norm;
            norm = 1./K_FUNC::E_max(sqrt(v2x*v2x+v2y*v2y+v2z*v2z), 1.e-10);
            v2x = v2x * norm;
            v2y = v2y * norm;
            v2z = v2z * norm;
                        
            //printf("ind=%d, xm=%f %f %f, xc1=%f %f %f, xc2=%f %f %f\n", ind, xm, ym, zm, xc1, yc1, zc1, xc2, yc2, zc2);
            //printf("ind=%d, norm1=%f %f %f\n", ind, sx1, sy1,sz1);
            
            E_Float ps = v1x*v2x+v1y*v2y+v1z*v2z;
            // si ps est positif, alors alpha est entre [0,90] ou [270,360]
            
            E_Float pvx = v1y*v2z-v1z*v2y;
            E_Float pvy = v1z*v2x-v1x*v2z;
            E_Float pvz = v1x*v2y-v1y*v2x;
            
            E_Float pv2x = v1y*sz1-v1z*sy1;
            E_Float pv2y = v1z*sx1-v1x*sz1;
            E_Float pv2z = v1x*sy1-v1y*sx1;
            
            E_Float ps2 = pvx*pv2x+pvy*pv2y+pvz*pv2z;
            // si ps2 est positif, alpha est entre [0-180]
            
            // acos retourne un angle entre [0,180]
            E_Float alpha = acos(ps);
            //printf("ind=%d, rough=%f, ps=%f, ps2=%f\n", ind, alpha*180./K_CONST::E_PI, ps, ps2);
            if (K_FUNC::E_abs(ps2) < 1.e-3 && K_FUNC::E_abs(ps+1.) < 1.e-3) alpha = K_CONST::E_PI;
            else if (K_FUNC::E_abs(ps2) < 1.e-3 && K_FUNC::E_abs(ps-1.) < 1.e-3) alpha = 0.;
            else if (ps >=0 && ps2 >=0)
            {
              // alpha must be in 0-90
              if (alpha > K_CONST::E_PI*0.5) printf("warning1\n");
            }
            else if (ps <=0 && ps2 >=0)
            {
              // alpha must be in 90-180
              if (alpha < K_CONST::E_PI*0.5) printf("warning2\n"); 
            }
            else if (ps >=0 && ps2 <=0)
            {
              // alpha must be in 180-270
              alpha = 2*K_CONST::E_PI-alpha;
              if (alpha < 3.*K_CONST::E_PI*0.5) { printf("warning3\n");  }
            }
            else
            {
              // alpha must be in 270-360
              alpha = 2.*K_CONST::E_PI-alpha;
              if (alpha > 3.*K_CONST::E_PI*0.5) { printf("warning4\n"); }
            }
            
            //printf("ind=%d: alpha=%f %f %f\n", ind, alpha*180./K_CONST::E_PI, ps, ps2);
            E_Float sas2 = K_FUNC::E_abs(sin(alpha*0.5));
            sas2 = std::max(sas2, 0.5);
            sas2 = std::min(sas2, 2.);
            fouth[ind] = 1./sas2;
            //printf("ind=%d: sas2=%f, fout=%f\n", ind,sas2,fout[ind]);
            
            // calcul de kappa
            //E_Float kappaL = 0.2;
            //E_Float kappaP = 1.6;
            E_Float kappa;
            if (alpha <= K_CONST::E_PI) kappa = ((1.-kappaP)/K_CONST::E_PI)*alpha+kappaP;
            else kappa = ((kappaL-1.)/K_CONST::E_PI)*alpha+2.-kappaL;
            //E_Float kappa = 1.3 + (0.7-1.3)/(2.*K_CONST::E_PI)*alpha;
            //E_Float kappa = 0.6 + (1.4-0.6)/(2.*K_CONST::E_PI)*alpha;
            //printf("ind=%d, kappa=%f, alpha=%f\n",ind, kappa, alpha*180./K_CONST::E_PI);
            if (kappaType == 1) fouth[ind] = kappa*fouth[ind];

            // Champ 2
            // Lissage en fonction de alpha (pas suffisant car trop ponctuel)
            //if (alpha > K_CONST::E_PI) foute[ind] = 0.;
            //else foute[ind] = 3.*(K_CONST::E_PI-alpha)/K_CONST::E_PI;
            // Pas de lissage
            //foute[ind] = 0.;
            // Lissage en fonction de alpha (pas suffisant car trop ponctuel)
            foute[ind] = 3.*K_FUNC::E_abs(K_CONST::E_PI-alpha)/K_CONST::E_PI;
          }
        }
      }
    }
  }
  PyObject* tpl = K_ARRAY::buildArray(fout, "ht,hl", *cn, -1, eltType);
  RELEASESHAREDU(array, f, cn);
  RELEASESHAREDU(normales, fn, cnn);
  return tpl;
}
