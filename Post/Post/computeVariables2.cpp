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

// Recupere les pointeurs si grandeurs conservatives
#define GETPTRS \
  rop  = (E_Float*)f.begin(posro);  \
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

# include "post.h"

using namespace K_FLD;
using namespace std;

// =================================================================================
// Compute variables - in place version 
//==================================================================================
PyObject* K_POST::computeVariables2(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars0;
  E_Float gamma, rgp, s0, betas, Cs, mus, Ts;

  if (!PYPARSETUPLEF(args,
                    "OOddddddd", "OOfffffff",
                    &array, &vars0, &gamma, &rgp, &s0, &betas, &Cs, &mus, &Ts))
  {
      return NULL;
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Check array
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, "computeVariable2: invalid array.");
    return NULL;
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Extrait les variables a calculer de la chaine vars0. 
  // Insere dans vars uniquement celles qui seront effectivement calculees
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  char varStringOut[K_ARRAY::VARSTRINGLENGTH]; varStringOut[0] = '\0';
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
        char* str = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0));
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
  else // PyList_Check(vars0) == 0
  {
    if (PyString_Check(vars0))
    {
      char* str = PyString_AsString(vars0);
      checkAndExtractVariables(str, vars, varStringOut);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(vars0)) 
    {
      char* str = PyBytes_AsString(PyUnicode_AsUTF8String(vars0));
      checkAndExtractVariables(str, vars, varStringOut); 
    }
#endif
    else printf("Warning: computeVariables: varname must be a string. Skipped...\n");
  }
  E_Int nvarout = vars.size(); // variables a calculer
  if (nvarout == 0)
  {
    RELEASESHAREDB(res, array, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "computeVariables2: no variable computed.");
    return NULL;
  }

  RELEASESHAREDB(res, array, f, cn);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Determine les variables a notre disposition
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  E_Int typeVars = -1; // 0: conservatives, 1: ro, u, T
  E_Int posro, posrou, posrov, posrow, posroe;
  E_Int posu, posv, posw, post;

  posro  = K_ARRAY::isDensityPresent(varString); 
  posrou = K_ARRAY::isMomentumXPresent(varString); 
  posrov = K_ARRAY::isMomentumYPresent(varString); 
  posrow = K_ARRAY::isMomentumZPresent(varString); 
  posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);

  // Test ro, rou, roe 
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
      printf("Warning: computeVariables2: one conservative or primitive field was not found. Variables are not computed.\n");
      Py_INCREF(Py_None);
      return Py_None;
    }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Evaluation des nouveaux champs 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  E_Int varsSize = vars.size();
  for (E_Int v = 0; v < varsSize; v++) 
  { 
      E_Int posvar = K_ARRAY::isNamePresent(vars[v],varString);

      if (posvar == -1) 
      {
          K_ARRAY::addFieldInArray(array, vars[v]);
      }

      E_Int res2 = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
					    cn, eltType); // utile si posvar > 0 ??

      posvar = K_ARRAY::isNamePresent(vars[v],varString);
      posvar++;

      if  (typeVars == 0) // cons.  
          computeCompVars(*f, posvar, vars[v], posro, posrou, posrov, 
   	      	              posrow, posroe, gamma, rgp, s0, betas, Cs);
      else if (typeVars == 1) // ro,u,T 
          computeCompVars2(*f, posvar, vars[v], posro, posu, posv, 
   	      	               posw, post, gamma, rgp, s0, betas, Cs);

      RELEASESHAREDB(res2, array, f, cn);
  }
  for (E_Int v = 0; v < varsSize; v++) delete [] vars[v];
  Py_INCREF(Py_None);
  return Py_None;
}

// =================================================================================
// Calcule les variables composees (a partir des variables cons.)
//==================================================================================
E_Int K_POST::computeCompVars(const FldArrayF& f,  const E_Int posnew, 
                                    char* varnew,  const E_Int posro,
                              const E_Int posrou,  const E_Int posrov,
                              const E_Int posrow,  const E_Int posroe,
                              const E_Float gamma, const E_Float rgp,
                              const E_Float s0,    const E_Float betas, 
                              const E_Float Cs)

{

  E_Float *rop, *roup, *rovp, *rowp, *roEp;
  E_Float *varp = (E_Float*)f.begin(posnew);
  E_Int    npts = f.getSize();
  bool     fail = false;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // VelocityX
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (K_STRING::cmp(varnew, "VelocityX") == 0) //vitesse absolue vx 
  {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rou, vx;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          ro  = rop[i];
          rou = roup[i];
          VELOCITYX;
          varp[i] = vx;
        }
      }    
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // VelocityY
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "VelocityY") == 0) //vitesse absolue vy
  {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, rov, vy;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          ro  = rop[i];
          rov = rovp[i];
          VELOCITYY;
          varp[i] = vy;
        }
      }    
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // VelocityZ
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "VelocityZ") == 0) //vitesse absolue vz
  {
      GETPTRS;
#pragma omp parallel
      {
        E_Float ro, row, vz;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          ro  = rop[i];
          row = rowp[i];
          VELOCITYZ;
          varp[i] = vz;
        }
      }    
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // VelocityMagnitude
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "VelocityMagnitude") == 0)// vitesse abs: module
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
          varp[i] = MAGNITUDE(vx,vy,vz);
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Pression statique
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "Pressure") == 0) // pression statique 
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
          varp[i] = p;
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Temperature
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "Temperature") == 0) //temperature statique
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
          varp[i] = t;
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Entropy
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "Entropy") == 0)
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
          varp[i] = s;
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Enthalpy
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "Enthalpy") == 0)
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
          varp[i] = h;
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Mach
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "Mach") == 0)
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
          varp[i] = mach;
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Viscosité Moléculaire 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "ViscosityMolecular") == 0) //viscosite du fluide 
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
          varp[i] = mu;
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Pression d'arrêt
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "PressureStagnation") == 0) //pression d'arret
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
          varp[i] = p*pow(1.+gam4*mach*mach, gam6);
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Température d'arrêt 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "TemperatureStagnation") == 0) //temperature d'arret
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
          varp[i] = t*(1.+gam4*mach*mach);
        }
      }
    }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Pression dynamique 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if (K_STRING::cmp(varnew, "PressureDynamic") == 0) //pression dynamique 
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
          varp[i] = 0.5*gamma*p*mach*mach;
        }
      }
    }
    else fail = true;
  
  if (fail == true) return 0;
  else return 1;
}
// =================================================================================
// Calcule les variables composées (à partir des variables ro, u, T)
// =================================================================================
E_Int K_POST::computeCompVars2(const FldArrayF& f,    const E_Int posnew,
			             char* varnew,    const E_Int posro,
                               const E_Int posu,      const E_Int posv,
                               const E_Int posw,      const E_Int post,
                               const E_Float gamma,   const E_Float rgp,
                               const E_Float s0,      const E_Float betas, 
                               const E_Float Cs)
{ 
  E_Float *rop, *up, *vp, *wp, *tp;
  E_Float *varp = (E_Float*)f.begin(posnew);
  E_Int    npts = f.getSize();
  bool     fail = false;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // VelocityX
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (K_STRING::cmp(varnew, "VelocityX") == 0)
  {
      GETPTRS2;
#pragma omp parallel
      {
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          varp[i] = up[i];
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // VelocityY
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "VelocityY") == 0)
  {
      GETPTRS2;
#pragma omp parallel
      {
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          varp[i] = vp[i];
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // VelocityZ
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "VelocityZ") == 0)
  {
      GETPTRS2;
#pragma omp parallel
      {
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          varp[i] = wp[i];
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // VelocityMagnitude
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "VelocityMagnitude") == 0)// vitesse abs: module
  {
      GETPTRS2;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, t;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          varp[i] = MAGNITUDE(vx,vy,vz);
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Pressure
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "Pressure") == 0) // pression statique 
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
          varp[i] = p;
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Temperature
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "Temperature") == 0) //temperature statique
  {
      GETPTRS2;
#pragma omp parallel
      {
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          varp[i] = tp[i];
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Entropy
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "Entropy") == 0)
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
          varp[i] = s;
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Enthalpy
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "Enthalpy") == 0)
  {
      GETPTRS2;
      GAM6;
#pragma omp parallel
      {
        E_Float ro, vx, vy, vz, p, t, h;
#pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          GETPRIM;
          PRESSURE2;
          ENTHALPY;
          varp[i] = h;
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Mach
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "Mach") == 0)
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
          varp[i] = mach;
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Viscosité du fluide
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "ViscosityMolecular") == 0) //viscosite du fluide 
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
          varp[i] = mu;
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Pression d'arrêt 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "PressureStagnation") == 0) //pression d'arret
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
          varp[i] = p*pow(1.+gam4*mach*mach, gam6);
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Température d'arrêt 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "TemperatureStagnation") == 0) //temperature d'arret
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
          varp[i] = t*(1.+gam4*mach*mach);
        }
      }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Pression dynamique 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (K_STRING::cmp(varnew, "PressureDynamic") == 0) //pression dynamique 
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
          varp[i] = 0.5*gamma*p*mach*mach;
        }
      }
  }
  else fail = true;
  
  if (fail == true) return 0;
  else return 1;
}
