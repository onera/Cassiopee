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
// compute variables such as p, primitive variables

// Recupere les pointeurs si grandeurs conservatives
#define GETPTRS \
  rop = (E_Float*)f.begin(posro); \
  roup = (E_Float*)f.begin(posrou); \
  rovp = (E_Float*)f.begin(posrov); \
  rowp = (E_Float*)f.begin(posrow); \
  roEp = (E_Float*)f.begin(posroe);

// Recupere les pointeurs si grandeurs ro, u, T
#define GETPTRS2 \
  rop = (E_Float*)f.begin(posro); \
  up = (E_Float*)f.begin(posu); \
  vp = (E_Float*)f.begin(posv); \
  wp = (E_Float*)f.begin(posw); \
  tp = (E_Float*)f.begin(post);

#define GETCONS \
  ro = rop[i]; rou = roup[i]; rov = rovp[i]; row = rowp[i]; roE = roEp[i];

#define GETPRIM \
  ro = rop[i]; vx = up[i]; vy = vp[i]; vz = wp[i]; t = tp[i];

#define VELOCITYX vx = rou/ro;
#define VELOCITYY vy = rov/ro;
#define VELOCITYZ vz = row/ro;
#define VELOCITY vx = rou/ro; vy = rov/ro; vz = row/ro;
#define MAGNITUDE(vx,vy,vz) sqrt(vx*vx+vy*vy+vz*vz);
#define PRESSURE p = (gamma-1.)*(roE-0.5*ro*(vx*vx+vy*vy+vz*vz)); // fonction de vx,vy,vz
#define TEMPERATURE t = p / (ro*rgp);
#define PRESSURE2 p = t*ro*rgp;
#define GAM4 E_Float gam4 = 0.5*(gamma-1.);
#define GAM5 E_Float gam5 = rgp *gamma/(gamma-1);
#define GAM6 E_Float gam6 = gamma/(gamma-1.);
#define ENTROPY s = s0 + gam5*log(t)-rgp*log(p);
#define ENTHALPY h = gam6*p/ro;
#define MACH mach = sqrt((vx*vx+vy*vy+vz*vz)*ro/(gamma*p));
#define MU mu = betas*sqrt(t)/(1.+Cs/t);

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include "post.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
void computeVelocity__(E_Int npts, 
  const E_Float* ro, const E_Float* rou, const E_Float* rov, const E_Float* row,
  E_Float* vx, E_Float* vy, E_Float* vz)
{
  #pragma omp parallel
  {
    E_Float roi;
    #pragma omp for
    for (E_Int i = 0; i < npts; i++)
    {
      roi = 1. / ro[i];
      vx[i] = rou[i] * roi;
      vy[i] = rov[i] * roi;
      vz[i] = row[i] * roi;
    }
  }
}

//=============================================================================
void computePressure__(E_Int npts, E_Float gamma, 
                       const E_Float* vx, const E_Float* vy, const E_Float* vz, 
                       const E_Float* ro, const E_Float* roe,
                       E_Float* p)
{
  E_Float gam1 = gamma - 1.;
  E_Float ONEHALF = 0.5;

  #pragma omp parallel
  {
    E_Float vx2, vy2, vz2;
    #pragma omp for
    for (E_Int i = 0; i < npts; i++)
    {
      vx2 = vx[i] * vx[i];
      vy2 = vy[i] * vy[i];
      vz2 = vz[i] * vz[i];
      p[i] = gam1 * (roe[i]-ONEHALF*ro[i]*(vx2 + vy2 + vz2));
    }
  }
}

//=============================================================================
/* Compute variables */
//=============================================================================
PyObject* K_POST::computeVariables(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars0;
  E_Float gamma, rgp, s0, betas, Cs, mus, Ts;
  if (!PYPARSETUPLE_(args, OO_ RRRR_ RRR_, 
                    &array, &vars0, &gamma, &rgp, &s0, &betas, &Cs, &mus, &Ts))
  {
      return NULL;
  }

  // calcul de betas donne par mus,Ts
  if (K_FUNC::fEqualZero(mus) == false && 
      K_FUNC::fEqualZero(Ts) == false)
  {
    betas = mus*(Ts+Cs)/(Ts*sqrt(Ts));
    //printf("Info: computeVariables: betas = %12.15e\n", betas);
  }

  // Check array
  char* varString; char* eltType;
  char varStringOut[K_ARRAY::VARSTRINGLENGTH]; varStringOut[0] = '\0';
  FldArrayF* f; FldArrayI* c;
  E_Int res; 
  E_Int ni, nj, nk; // number of points of array
  res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, c, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeVariables: unknown type of array.");
    return NULL;
  } 
     
  // Extrait les variables a calculer de la chaine vars0. 
  // Insert dans vars uniquement celles qui seront effectivement calculees
  vector<char*> vars;
  if (PyList_Check(vars0) != 0)
  {
    for (int i = 0; i < PyList_Size(vars0); i++)
    {
      PyObject* tpl0 = PyList_GetItem(vars0, i);
      if (PyString_Check(tpl0))
      {
        char* str = PyString_AsString(tpl0);
        char tmpVarString[K_ARRAY::VARSTRINGLENGTH];
        short ok = checkAndExtractVariables(str, vars, tmpVarString);
        if (ok != 0) 
        {
          if (varStringOut[0] == '\0') strcpy(varStringOut, tmpVarString);
          else {strcat(varStringOut, ","); strcat(varStringOut, tmpVarString);}
        }
      }
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(tpl0)) 
      {
        char* str = (char*)PyUnicode_AsUTF8(tpl0);
        char tmpVarString[K_ARRAY::VARSTRINGLENGTH];
        short ok = checkAndExtractVariables(str, vars, tmpVarString);
        if (ok != 0) 
        {
          if (varStringOut[0] == '\0') strcpy(varStringOut, tmpVarString);
          else {strcat(varStringOut, ","); strcat(varStringOut, tmpVarString);}
        }
      }
#endif
      else  
      {
        printf("Warning: computeVariables: varname must be a string. Skipped...\n");
      }
    }
  }
  else 
  {
    if (PyString_Check(vars0)) 
    {  
      char* str = PyString_AsString(vars0);
      checkAndExtractVariables(str, vars, varStringOut);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(vars0)) 
    {
      char* str = (char*)PyUnicode_AsUTF8(vars0);
      checkAndExtractVariables(str, vars, varStringOut);
    }
#endif
    else
      printf("Warning: computeVariables: varname must be a string. Skipped...\n");
  }
  E_Int nvarout = vars.size();// variables a calculer
  if (nvarout == 0)
  {
    RELEASESHAREDB(res, array, f, c); 
    PyErr_SetString(PyExc_TypeError,
                    "computeVariables: no variable computed.");
    return NULL;
  }

  // Determine les variables a notre disposition
  E_Int typeVars = -1; // 0: conservatives, 1: ro, u, T
  E_Int posro, posrou, posrov, posrow, posroe;
  E_Int posu, posv, posw, post;
  posro = K_ARRAY::isDensityPresent(varString); 
  posrou = K_ARRAY::isMomentumXPresent(varString); 
  posrov = K_ARRAY::isMomentumYPresent(varString); 
  posrow = K_ARRAY::isMomentumZPresent(varString); 
  posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
  if (posro >= 0 && posrou >= 0 && posrov >= 0 && posrow >= 0 && posroe >= 0)
  {
    posro++; posrou++; posrov++; posrow++; posroe++;
    typeVars = 0;
  }
  else
  {
    // Essai ro, u, T
    posu = K_ARRAY::isVelocityXPresent(varString);
    posv = K_ARRAY::isVelocityYPresent(varString); 
    posw = K_ARRAY::isVelocityZPresent(varString); 
    post = K_ARRAY::isTemperaturePresent(varString);

    if (posro >= 0 && posu >= 0 && posv >= 0 && posw >= 0 && post >= 0)
    {
      posro++; posu++; posv++; posw++; post++;
      typeVars = 1;
    }
    else
    {
      printf("Warning: computeVariables: one conservative or primitive field was not found. Variables are not computed.\n");
      RELEASESHAREDB(res, array, f, c); return array;
    }
  }

  FldArrayF* fnew = new FldArrayF(f->getSize(), nvarout);

  // calcul des variables composees
  E_Int ok = 0;
  if (typeVars == 0) // cons
    ok = computeCompVariables(*f, posro, posrou, posrov, posrow, posroe, 
                              gamma, rgp, s0, betas, Cs, vars, *fnew); 
  else if (typeVars == 1) // ro,u,T
    ok = computeCompVariables2(*f, posro, posu, posv, posw, post, 
                               gamma, rgp, s0, betas, Cs, vars, *fnew); 
  // nettoyage...
  E_Int varsSize = vars.size();
  for (E_Int v = 0; v < varsSize; v++) delete [] vars[v];

  if (ok == 0) // erreur de developpt
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeVariables: invalid string.\n");
    delete fnew; RELEASESHAREDB(res, array, f, c); return NULL;
  }

  if (res == 1)
  {
    PyObject* tpl = K_ARRAY::buildArray(*fnew, varStringOut, 
                                        ni, nj, nk);
    delete fnew; 
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else
  {
    PyObject* tpl = K_ARRAY::buildArray(*fnew, varStringOut, 
                                        *c, -1, eltType);
    delete fnew; 
    RELEASESHAREDU(array, f, c);
    return tpl;
  }
}
//-----------------------------------------------------------------------------
/* Extrait de la chaine vars0 les variables a calculer. Une verification 
   est effectuee sur les noms de variables. La chaine varStringOut est 
   aussi construite pour l'array de sortie 
   IN: vars0: chaine contenant les variables a extraire
   OUT: vars: vecteur contenant les variables a calculer
   OUT: varStringOut: chaine de variables calculees pour l'array de sortie
   retourne 0 si aucune variable n a ete trouvee.
*/
//-----------------------------------------------------------------------------
short K_POST::checkAndExtractVariables(char* vars0, vector<char*>& vars, 
                                       char* varStringOut)
{
  vector<char*> tmpvars;
  //extraction de toutes les variables
  K_ARRAY::extractVars(vars0, tmpvars);

  E_Int maxStringLength = 80;

  E_Int first = 1; // test si premier element de vars
  // verification
  E_Int tmpvarsSize = tmpvars.size();
  short ok = 0;

  for (E_Int v = 0 ; v < tmpvarsSize; v++)
  {
    if (K_STRING::cmp(tmpvars[v], "VelocityX") == 0 ||
        K_STRING::cmp(tmpvars[v], "VelocityY") == 0 ||
        K_STRING::cmp(tmpvars[v], "VelocityZ") == 0 ||
        K_STRING::cmp(tmpvars[v], "VelocityMagnitude") == 0 ||
        K_STRING::cmp(tmpvars[v], "Pressure") == 0 ||
        K_STRING::cmp(tmpvars[v], "Temperature") == 0 ||
        K_STRING::cmp(tmpvars[v], "Entropy") == 0 ||
        K_STRING::cmp(tmpvars[v], "Enthalpy") == 0 ||
        K_STRING::cmp(tmpvars[v], "Mach") == 0 ||
        K_STRING::cmp(tmpvars[v], "ViscosityMolecular") == 0 ||
        K_STRING::cmp(tmpvars[v], "PressureStagnation") == 0  ||
        K_STRING::cmp(tmpvars[v], "TemperatureStagnation") == 0 ||
        K_STRING::cmp(tmpvars[v], "PressureDynamic") == 0 
//         K_STRING::cmp(tmpvars[v],"RotatingVelocityX") == 0 ||
//         K_STRING::cmp(tmpvars[v],"RotatingVelocityY") == 0 ||
//         K_STRING::cmp(tmpvars[v],"RotatingVelocityZ") == 0 ||
//         strcmp(tmpvars[v],"RotatingVelocityMagnitude") == 0
      )
    {
      char* tmp = new char[maxStringLength];
      strcpy(tmp, tmpvars[v]);
      vars.push_back(tmp);
      if (first == 1)
      {
        strcpy(varStringOut, tmpvars[v]);
        first = 0;
        ok = 1;
      }
      else 
      {
        strcat(varStringOut, ",");
        strcat(varStringOut, tmpvars[v]);
      }
    }
    else 
    {
      printf("Warning: computeVariables: %s is not a valid string name. Skipped...\n", tmpvars[v]);
    }
  }
  for (E_Int v = 0; v < tmpvarsSize; v++) delete [] tmpvars[v];
  
  return ok;
}

//-----------------------------------------------------------------------------
// Calcule les variables composees (a partir des variables conservatives)
//-----------------------------------------------------------------------------
E_Int K_POST::computeCompVariables(const FldArrayF& f, const E_Int posro,
                                   const E_Int posrou, const E_Int posrov,
                                   const E_Int posrow, const E_Int posroe,
                                   const E_Float gamma, const E_Float rgp,
                                   const E_Float s0,  
                                   const E_Float betas, const E_Float Cs,
                                   vector<char*>& vars,
                                   FldArrayF& fnew)
{ 
  E_Int npts = f.getSize();
  E_Float *rop, *roup, *rovp, *rowp, *roEp;

  E_Int varsSize = vars.size();
  bool fail = false;
  for (E_Int v = 0 ; v < varsSize; v++)
  {
    E_Float* fnewv = fnew.begin(v+1);

    if (K_STRING::cmp(vars[v], "VelocityX") == 0) //vitesse absolue vx 
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rou, vx;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          ro = rop[i];
          rou = roup[i];
          VELOCITYX;
          fnewv[i] = vx;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "VelocityY") == 0) //vitesse absolue vy
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rov, vy;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          ro = rop[i];
          rov = rovp[i];
          VELOCITYY;
          fnewv[i] = vy;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "VelocityZ") == 0) //vitesse absolue vz
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, row, vz;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          ro = rop[i];
          row = rowp[i];
          VELOCITYZ;
          fnewv[i] = vz;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "VelocityMagnitude") == 0)// vitesse abs: module
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          fnewv[i] = MAGNITUDE(vx,vy,vz);
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Pressure") == 0) // pression statique 
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          fnewv[i] = p;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Temperature") == 0) //temperature statique
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p, t;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          TEMPERATURE;
          fnewv[i] = t;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Entropy") == 0)
    {
      GETPTRS;
      GAM5;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p, t, s;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          TEMPERATURE;
          ENTROPY;
          fnewv[i] = s;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Enthalpy") == 0)
    {
      GETPTRS;
      GAM6;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p, h;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          ENTHALPY;
          fnewv[i] = h;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Mach") == 0)
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p, mach;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          MACH;
          fnewv[i] = mach;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "ViscosityMolecular") == 0) //viscosite du fluide 
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p, t, mu;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          TEMPERATURE;
          MU;
          fnewv[i] = mu;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "PressureStagnation") == 0) //pression d'arret
    {
      GETPTRS;
      GAM4;
      GAM6;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p, mach;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          MACH;
          fnewv[i] = p*pow(1.+gam4*mach*mach, gam6);
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "TemperatureStagnation") == 0) //temperature d'arret
    {
      GETPTRS;
      GAM4;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p, t, mach;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          TEMPERATURE;
          MACH;
          fnewv[i] = t*(1.+gam4*mach*mach);
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "PressureDynamic") == 0) //pression dynamique 
    {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rou, rov, row, roE, vx, vy, vz, p, mach;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETCONS;
          VELOCITY;
          PRESSURE;
          MACH;
          fnewv[i] = 0.5*gamma*p*mach*mach;
        }
      }
    }
    else fail = true;
  }
  if (fail == true) return 0;
  else return 1;
}

//-----------------------------------------------------------------------------
// Calcule les variables composees (a partir des variables ro, u, T)
//-----------------------------------------------------------------------------
E_Int K_POST::computeCompVariables2(const FldArrayF& f, const E_Int posro,
                                    const E_Int posu, const E_Int posv,
                                    const E_Int posw, const E_Int post,
                                    const E_Float gamma, const E_Float rgp,
                                    const E_Float s0,
                                    const E_Float betas, const E_Float Cs,
                                    vector<char*>& vars,
                                    FldArrayF& fnew)
{ 
  E_Int npts = f.getSize();
  E_Float *rop, *up, *vp, *wp, *tp;
  bool fail = false;

  E_Int varsSize = vars.size();
  for (E_Int v = 0 ; v < varsSize; v++)
  {
    E_Float* fnewv = fnew.begin(v+1);

    if (K_STRING::cmp(vars[v], "VelocityX") == 0)
    {
      GETPTRS2;
#pragma omp parallel
      {
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          fnewv[i] = up[i];
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "VelocityY") == 0)
    {
      GETPTRS2;
#pragma omp parallel
      {
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          fnewv[i] = vp[i];
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "VelocityZ") == 0)
    {
      GETPTRS2;
#pragma omp parallel
      {
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          fnewv[i] = wp[i];
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "VelocityMagnitude") == 0)// vitesse abs: module
    {
      GETPTRS2;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, t;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          fnewv[i] = MAGNITUDE(vx,vy,vz);
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Pressure") == 0) // pression statique 
    {
      GETPTRS2;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, p, t;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          PRESSURE2;
          fnewv[i] = p;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Temperature") == 0) //temperature statique
    {
      GETPTRS2;
#pragma omp parallel
      {
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          fnewv[i] = tp[i];
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Entropy") == 0)
    {
      GETPTRS2;
      GAM5;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, p, t, s;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          PRESSURE2;
          ENTROPY;
          fnewv[i] = s;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Enthalpy") == 0)
    {
      GETPTRS2;
      GAM6;
#pragma omp parallel
      {
        E_Float ro, p, t, vx, vy, vz, h;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          PRESSURE2;
          ENTHALPY;
          fnewv[i] = h;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "Mach") == 0)
    {
      GETPTRS2;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, t, p, mach;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          PRESSURE2;
          MACH;
          fnewv[i] = mach;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "ViscosityMolecular") == 0) //viscosite du fluide 
    {
      GETPTRS2;
#pragma omp parallel
      {
        E_Float t, mu;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          t = tp[i];
          MU;
          fnewv[i] = mu;
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "PressureStagnation") == 0) //pression d'arret
    {
      GETPTRS2;
      GAM4;
      GAM6;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, p, t, mach;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          PRESSURE2;
          MACH;
          fnewv[i] = p*pow(1.+gam4*mach*mach, gam6);
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "TemperatureStagnation") == 0) //temperature d'arret
    {
      GETPTRS2;
      GAM4;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, p, t, mach;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          PRESSURE2;
          MACH;
          fnewv[i] = t*(1.+gam4*mach*mach);
        }
      }
    }
    else if (K_STRING::cmp(vars[v], "PressureDynamic") == 0) //pression dynamique 
    {
      GETPTRS2;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, p, t, mach;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          PRESSURE2;
          MACH;
          fnewv[i] = 0.5*gamma*p*mach*mach;
        }
      }
    }
    else fail = true;
  }
  if (fail == true) return 0;
  else return 1;
}
//=============================================================================
// Calcul de la vitesse absolue
// IN: f: champ (doit contenir les variables conservatives)
// OUT: velo: champ de vitesse
//=============================================================================
void K_POST::computeVelocity(const FldArrayF& f, 
                             const E_Int posro, const E_Int posrou, 
                             const E_Int posrov, const E_Int posrow,
                             FldArrayF& velo)
{
  if (velo.getSize() != 0) return; // deja calcule
  E_Int npts = f.getSize();

  velo.malloc(npts, 3);
  computeVelocity__(npts, f.begin(posro), f.begin(posrou), f.begin(posrov), f.begin(posrow),
    velo.begin(1), velo.begin(2), velo.begin(3));
}

//=============================================================================
// Calcul de la pression statique
// IN: f: champ (doit contenir les variables conservatives)
// IN: gamma: constante gaz parfait
// OUT: velo: champ de vitesse
// OUT: press: pression statique pour un gaz parfait
//=============================================================================
void K_POST::computePressure(const FldArrayF& f, 
                             const E_Int posro, const E_Int posrou, 
                             const E_Int posrov, const E_Int posrow,
                             const E_Int posroe, const E_Float gamma,
                             FldArrayF& velo, FldArrayF& press)
{
  if (press.getSize() != 0) return; // deja calcule
  
  E_Int npts = f.getSize();
  press.malloc(npts);
    
  computeVelocity(f, posro, posrou, posrov, posrow, velo);
  
  // calcul de p = gam1 * roe, e energie interne
  computePressure__(npts, gamma, velo.begin(1), velo.begin(2), velo.begin(3),
    f.begin(posro), f.begin(posroe), press.begin());
}

//=============================================================================
// Calcul de la temperature statique
// IN: f: champ (doit contenir les variables conservatives)
// IN: gamma, rgp
// OUT: velo, press, temp
//=============================================================================
void K_POST::computeTemperature(const FldArrayF& f, 
                                const E_Int posro, const E_Int posrou, 
                                const E_Int posrov, const E_Int posrow,
                                const E_Int posroe, const E_Float gamma,
                                const E_Float rgp,
                                FldArrayF& velo, FldArrayF& press,
                                FldArrayF& temp)
{
  if (temp.getSize() != 0) return; // deja calcule
  
  E_Int npts = f.getSize();
  temp.malloc(npts);
  computePressure(f, posro, posrou, posrov, posrow, posroe,
                  gamma, velo, press); 
  const E_Float* rop = f.begin(posro);
  E_Float* tempp = temp.begin();
  E_Float* pressp = press.begin();
  for (E_Int ind = 0; ind < npts; ind++)
    tempp[ind] = pressp[ind] / (rop[ind]*rgp);
}

//=============================================================================
// Calcul de la viscosite du fluide
//=============================================================================
void K_POST::computeMu(const FldArrayF& f, 
                       const E_Int posro, const E_Int posrou, 
                       const E_Int posrov, const E_Int posrow,
                       const E_Int posroe, 
                       const E_Float gamma, const E_Float rgp,
                       const E_Float betas, const E_Float Cs,
                       FldArrayF& velo, FldArrayF& press,
                       FldArrayF& temp, FldArrayF& mu)
{
  if (mu.getSize() != 0) return; // deja calcule
  
  E_Float t, inv;
  E_Int npts = f.getSize();
  mu.malloc(npts);
  computeTemperature(f, posro, posrou, posrov, posrow, posroe,
                     gamma, rgp, velo, press, temp);
  E_Float* tempp = temp.begin();
  E_Float* mup = mu.begin();
  for (E_Int ind = 0 ; ind < npts; ind++)
  {
    t = tempp[ind];
    inv = 1.+Cs/t;      
    mup[ind] = betas * sqrt(t) / inv;
  }
}

//=============================================================================
/* Calcul de l'entropie: s = s0 + rgp*gam/(gam-1)*log(T) - rgp*log(P)
   ou s0 = sref -  rgp*gam/(gam-1)*log(Tref) + rgp*log(Pref)
   IN: f: champ (doit contenir les variables conservatives)
   IN: gamma, rgp: constantes des gaz parfaits
*/
//=============================================================================
void K_POST::computeEntropy(const FldArrayF& f, 
                            const E_Int posro, const E_Int posrou, 
                            const E_Int posrov, const E_Int posrow,
                            const E_Int posroe, const E_Float gamma,
                            const E_Float rgp, const E_Float s0,
                            FldArrayF& velo,
                            FldArrayF& temp, FldArrayF& press,
                            FldArrayF& s)
{ 
  E_Int npts = f.getSize();
  E_Float cp = rgp * gamma / (gamma-1.);
  
  //calcul de la pression si pas deja calculee
  computePressure(f, posro, posrou, posrov, posrow, posroe,
                  gamma, velo, press); 
  //calcul de la temperature si pas deja calculee
  computeTemperature(f, posro, posrou, posrov, posrow, posroe,
                    gamma, rgp, velo, press, temp);
  E_Float* sp = s.begin();
  E_Float* tempp = temp.begin();
  E_Float* pressp = press.begin();
  for (E_Int ind = 0 ; ind < npts; ind++)
    sp[ind] = s0 + cp * log(tempp[ind]) - rgp * log(pressp[ind]);
}
//=============================================================================
// calcul de l'enthalpie
//=============================================================================
void K_POST::computeEnthalpy(const FldArrayF& f, 
                             const E_Int posro, const E_Int posrou, 
                             const E_Int posrov, const E_Int posrow,
                             const E_Int posroe, const E_Float gamma,
                             FldArrayF& velo, FldArrayF& press, FldArrayF& h)
{
  E_Int npts = f.getSize();

  // calcul de la pression si pas deja calculee
  computePressure(f, posro, posrou, posrov, posrow, posroe,
                  gamma, velo, press);
  E_Float gam6 = gamma/(gamma-1.);
  const E_Float* rop = f.begin(posro);
  E_Float* hp = h.begin();
  E_Float* pressp = press.begin();

  for (E_Int i = 0; i < npts; i++) hp[i] = gam6 * pressp[i]/rop[i];
}

//=============================================================================
// Calcul du mach 
//=============================================================================
void K_POST::computeMach(const FldArrayF& f, 
                         const E_Int posro, const E_Int posrou, 
                         const E_Int posrov, const E_Int posrow,
                         const E_Int posroe, const E_Float gamma,
                         FldArrayF& velo, FldArrayF& press, FldArrayF& mach)
{
  if (mach.getSize() != 0) return; // deja calcule
  
  E_Int npts = f.getSize(); mach.malloc(npts);
  computePressure(f, posro, posrou, posrov, posrow, posroe,
                  gamma, velo, press);
  E_Float u2, ainv;
  E_Float gammai = 1./gamma;
  E_Float* vx = velo.begin(1);
  E_Float* vy = velo.begin(2);
  E_Float* vz = velo.begin(3);
  const E_Float* rop = f.begin(posro);
  E_Float* pressp = press.begin();
  E_Float* machp = mach.begin();
  for (E_Int ind = 0; ind < npts; ind++)
  {
    u2 = vx[ind]*vx[ind] + vy[ind]*vy[ind] + vz[ind]*vz[ind];   
    ainv = gammai * rop[ind] / pressp[ind];
    machp[ind] = sqrt(u2 * ainv);
  }
}
