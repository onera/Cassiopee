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

//===========================================================================
PyObject* K_CONNECTOR::_computeFrictionVelocityIBM(PyObject* self, PyObject* args)
{
    PyObject *zone;
    char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
    if (!PYPARSETUPLE_(args, O_ SSS_,
                       &zone, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
        return NULL;

    vector<PyArrayObject*> hook;
    E_Int im, jm, km, cnSize, cnNfld;
    char* varString; char* eltType;
    vector<E_Float*> fields; vector<E_Int> locs;
    vector<E_Int*> cn;
    E_Int xyz = 1; E_Int loc = 0;
    E_Int res = K_PYTREE::getFromZone(zone, xyz, loc, varString, fields, locs,
                                      im, jm, km,
                                      cn, cnSize, cnNfld, eltType, hook,
                                      GridCoordinates,
                                      FlowSolutionNodes, FlowSolutionCenters);
    E_Int npts;
    if (res == 1) npts = im*jm*km;
    else npts = im;

    E_Int posDens = K_ARRAY::isNamePresent("Density",varString);
    E_Float* densPtr = fields[posDens];
    E_Int posu = K_ARRAY::isNamePresent("VelocityX",varString);
    E_Float* uPtr = fields[posu];
    E_Int posv = K_ARRAY::isNamePresent("VelocityY",varString);
    E_Float* vPtr = fields[posv];
    E_Int posw = K_ARRAY::isNamePresent("VelocityZ",varString);
    E_Float* wPtr = fields[posw];
    E_Int posPress = K_ARRAY::isNamePresent("Pressure",varString);
    E_Float* pressPtr = fields[posPress];
    E_Int posVisc = K_ARRAY::isNamePresent("ViscosityMolecular",varString);
    E_Float* viscPtr = fields[posVisc];
    E_Int posuTau = K_ARRAY::isNamePresent("utau",varString);
    E_Float* utauPtr = fields[posuTau];
    E_Int posyplus = K_ARRAY::isNamePresent("yplus",varString);
    E_Float* yplusPtr = fields[posyplus];

    E_Int posxPI = K_ARRAY::isNamePresent("CoordinateX_PI",varString);
    E_Float* xPI = fields[posxPI];
    E_Int posyPI = K_ARRAY::isNamePresent("CoordinateY_PI",varString);
    E_Float* yPI = fields[posyPI];
    E_Int poszPI = K_ARRAY::isNamePresent("CoordinateZ_PI",varString);
    E_Float* zPI = fields[poszPI];

    E_Int posxPW = K_ARRAY::isNamePresent("CoordinateX_PW",varString);
    E_Float* xPW = fields[posxPW];
    E_Int posyPW = K_ARRAY::isNamePresent("CoordinateY_PW",varString);
    E_Float* yPW = fields[posyPW];
    E_Int poszPW = K_ARRAY::isNamePresent("CoordinateZ_PW",varString);
    E_Float* zPW = fields[poszPW];

    if ( posDens == -1 || posu == -1 || posv == -1 || posw == -1 || posPress == -1 || posVisc == -1 || posuTau == -1 ||
     posyplus == -1 || posxPI== -1 || posyPI == -1 || poszPI == -1 || posxPW == -1 || posyPW == -1 || poszPW == -1 )
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeFrictionVelocityIBM: Some required quantities cannot be extracted from zone.");
      if (res == 2) delete [] eltType;
      RELEASESHAREDZ(hook, varString, (char*)NULL);
      return NULL;
    }

    E_Float ax, l1, l2, l3;
    E_Float roext, uext, pext, muext, yext, rowall, muwall;
    E_Float uscaln, un, vn, wn, ut, vt, wt, utau0;
    E_Float u, v, w, normb;
    E_Float b0, b1, b2, n0, n1, n2;

    //Lois de paroi : criteres d arret pour estimer le frottement par Newton
    E_Float newtoneps = 1.e-7; // critere d arret pour u+
    E_Float newtonepsprime = 1.e-12;// critere d arret pour la derivee
    E_Float fp, tp;

    E_Int ideb=0, ifin=npts;
    E_Int err=0, skip=0;

    E_Float* utau_vec    = utauPtr;
    E_Float* ro_vec      = densPtr;
    E_Float* yplus_vec   = yplusPtr;
    E_Float* press_vec   = pressPtr;
    FldArrayF aa_vec(npts), uext_vec(npts), nutcible_vec(npts);
    FldArrayF ut_vec(npts), vt_vec(npts), wt_vec(npts), mu_vec(npts), alpha_vec(npts), utauv_vec(npts);

    for (E_Int noind = 0; noind < npts; noind++)
    {
      roext = densPtr[noind];   // Densite du point interpole.
      pext  = pressPtr[noind];  // Pression du point interpole.
      muext  = viscPtr[noind];  // Viscosite du point interpole.
      // vitesse du pt ext
      u = uPtr[noind];
      v = vPtr[noind];
      w = wPtr[noind];
      # include "IBC/commonMuskerLaw_init_constantVars.h"
      // out= utau  et err
    }
    
    // Newton pour utau
    # include "IBC/commonMuskerLaw_Newton.h"
    // Compute the correct yplus
    for (E_Int noind = 0; noind < npts; noind++)
    {
      yplus_vec[noind] = aa_vec[noind]*utau_vec[noind];
    }
    
    if (res == 2) delete [] eltType;
    RELEASESHAREDZ(hook, varString, (char*)NULL);
    Py_INCREF(Py_None);
    return Py_None;
}
