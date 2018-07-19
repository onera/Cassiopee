/*    
    Copyright 2013-2018 Onera.

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

#define RELEASEDATA\
  RELEASESHAREDN(interpPtsCoordX, coordX);\
  RELEASESHAREDN(interpPtsCoordY, coordY);\
  RELEASESHAREDN(interpPtsCoordZ, coordZ);\
  RELEASESHAREDZ(hook, varStringD, eltTypeD);\

//=============================================================================
/* Recherche les donnees d interpolation et calcule les formules d interpolation
   directement pour une liste de points a partir d une zone donnee 
   Similaire a un extractMesh
   PENALISATION DU VOLUME DONNEUR : 
   INTERPOLATION DEPUIS UN BORD : vol += 1e3
   EXTRAPOLATION : vol += 1e6
   ORPHAN : vol+=1e12
*/
//=============================================================================
PyObject* K_CONNECTOR::transferFields(PyObject* self, PyObject* args)
{
  E_Float penaltyExtrap = 1e6;
  E_Float penaltyOrphan = 1e12;

  E_Int locDnr = 0;//localisation des champs dans la zone donneuse - en noeuds

  PyObject *interpPtsCoordX, *interpPtsCoordY, *interpPtsCoordZ; 
  PyObject *zoneD;//donor zone
  PyObject* hookADT;
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  PyObject* pyVariables;
  E_Int interporder;
  E_Int nature;// O: produit des cellN=0 -> donneur invalide; 1: cellN=0 ou 2 -> donneur invalide
  E_Int penalty;//1 : penalite sur le volume des pts ou cellules frontieres
  E_Float constraint;
  E_Int extrapOrder = 1;
  if (!PYPARSETUPLE(args, "OOOOllldOOsss", "OOOOiiidOOsss","OOOOlllfOOsss","OOOOiiifOOsss",
                    &zoneD, &interpPtsCoordX, &interpPtsCoordY, &interpPtsCoordZ,
                    &interporder, &nature, &penalty, &constraint, &hookADT, &pyVariables,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters)) return NULL;

  // Interpolation type
  K_INTERP::InterpAdt::InterpolationType interpType;
  E_Int nindi, ncfmax;
  switch (interporder)
  {
    case 2:
      interpType = K_INTERP::InterpAdt::O2CF;
      ncfmax = 8; nindi = 1;
      break;

    case 3:
      interpType = K_INTERP::InterpAdt::O3ABC;
      ncfmax = 9; nindi = 1;
      break;

    case 5:
      interpType = K_INTERP::InterpAdt::O5ABC;
      ncfmax = 15; nindi = 1;
      break;
        
    default:
      printf("Warning: transferVariables: unknown interpolation order.");
      printf(" Set to 2nd order.\n");
      interpType = K_INTERP::InterpAdt::O2CF;
      ncfmax = 8; nindi = 1;
      break;
  } 

  /*--------------------------------------------------*/
  /* Extraction des infos sur les points a interpoler */
  /*--------------------------------------------------*/
  FldArrayF* coordX;
  E_Int resn = K_NUMPY::getFromNumpyArray(interpPtsCoordX, coordX, true);
  if (resn == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "transferFields: 2nd arg must be a numpy of floats.");
    return NULL;
  }
  FldArrayF* coordY;
  resn = K_NUMPY::getFromNumpyArray(interpPtsCoordY, coordY, true);  
  if (resn == 0)
  {
    RELEASESHAREDN(interpPtsCoordX, coordX);
    PyErr_SetString(PyExc_TypeError, 
                    "transferFields: 3rd arg must be a numpy of floats.");
    return NULL;
  }
  FldArrayF* coordZ;
  resn = K_NUMPY::getFromNumpyArray(interpPtsCoordZ, coordZ, true);  
  if (resn == 0)
  {
    RELEASESHAREDN(interpPtsCoordX, coordX);
    RELEASESHAREDN(interpPtsCoordY, coordY);
    PyErr_SetString(PyExc_TypeError, 
                    "transferFields: 4th arg must be a numpy of floats.");
    return NULL;
  }

  /*--------------------------------------------------*/
  /* recupere les champs de la zone donneuse (nodes)  */
  /*--------------------------------------------------*/
  char* varStringD;
  vector<E_Int> locsD;
  vector<E_Int*> cnd;
  vector<E_Float*> fieldsD;
  E_Int imd, jmd, kmd, cnSizeD, cnNfldD;
  char* eltTypeD; vector<PyArrayObject*> hook;
  E_Int typeZoneD = K_PYTREE::getFromZone(zoneD, 1, locDnr, varStringD,
                                          fieldsD, locsD, imd, jmd, kmd,
                                          cnd, cnSizeD, cnNfldD, 
                                          eltTypeD, hook,
                                          GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

  if ( typeZoneD == 0)
  {
    RELEASESHAREDN(interpPtsCoordX, coordX);
    RELEASESHAREDN(interpPtsCoordY, coordY);
    RELEASESHAREDN(interpPtsCoordZ, coordZ);
    PyErr_SetString(PyExc_TypeError, 
                    "transferFields: 1st arg is not a valid zone.");
    return NULL;
  }   // cas seulement non structure : ncf a 4 (minimum)
  else if (typeZoneD==2)
  {
    RELEASEDATA;
    PyErr_SetString(PyExc_TypeError, 
                    "transferFields: only valid for structured donor zones.");
    return NULL;   
  }
  E_Int nfld = fieldsD.size();
  E_Int npts = imd*jmd*kmd;
  E_Float** rakeFieldsD = new E_Float*[nfld];
  for (E_Int i =0; i < nfld; i++)
    rakeFieldsD[i] = fieldsD[i];
  // FldArrayF* donorFields = new FldArrayF(npts, nfld, rakeFieldsD, true, true);
  FldArrayF donorFields(npts, nfld, rakeFieldsD, true, true);

  //recup de l ADT
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    void** packet = (void**) PyCObject_AsVoidPtr(hookADT);
#else
    void** packet = (void**) PyCapsule_GetPointer(hookADT, NULL);
#endif

  E_Int* typeHookp = (E_Int*)packet[0];
  E_Int typeHook = typeHookp[0];
  if  ( typeHook != 1 ) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "transferFields: 5th arg must be a hook on an ADT.");
    RELEASEDATA;
    return NULL;
  }
  E_Int s1 = typeHookp[1];
  K_INTERP::InterpAdt* adt;
  for (E_Int i = 0; i < s1; i++)
    adt = (K_INTERP::InterpAdt*)(packet[i+1]);

  E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringD);
  E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringD);
  E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringD);
  E_Int poscd = K_ARRAY::isCellNatureField2Present(varStringD);
  if ( posxd == -1 || posyd == -1 || poszd==-1)
  {
    PyErr_SetString(PyExc_TypeError, "transferFields: no coordinates found for donor zone.");
    RELEASEDATA;
    return NULL;
  }
  if ( poscd == -1)
  {
    PyErr_SetString(PyExc_TypeError, "transferFields: no cellN found in donor zone.");
    RELEASEDATA;
    return NULL;
  }
  posxd++; posyd++; poszd++; poscd++;
  /*---------------------------------------------------*/
  /*  Extrait les positions des variables a transferer */
  /*---------------------------------------------------*/
  vector<E_Int> posvars0;
  vector<char*> listOfVars;
  E_Int lenVarString = 0;
  if (PyList_Check(pyVariables) != 0)
  {
    int nvariables = PyList_Size(pyVariables);
    if ( nvariables > 0 )
    {
      for (int i = 0; i < nvariables; i++)
      {
        PyObject* tpl0 = PyList_GetItem(pyVariables, i);
        if (PyString_Check(tpl0) == 0)
          PyErr_Warn(PyExc_Warning, "transferFields: variable must be a string. Skipped.");
        else 
        {
          char* varname = PyString_AsString(tpl0);        
          E_Int posvd = K_ARRAY::isNamePresent(varname, varStringD);      
          if (posvd != -1) 
          {
            posvars0.push_back(posvd+1);            
            listOfVars.push_back(varname);
            lenVarString += strlen(varname);
          }
        }
      }
    }
    else
    {
      PyErr_SetString(PyExc_TypeError, "transferFields: variables to be transfered is an empty list.");
      RELEASEDATA;
      return NULL;
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "transferFields: variables to be transfered must be defined by a list.");
    RELEASEDATA;
    return NULL;
  }
  E_Int nvars = listOfVars.size();
  lenVarString += nvars+strlen("donorVol");// les virgules
  char* varStringOut = new char[lenVarString];

  for (E_Int nov = 0; nov < listOfVars.size(); nov++)
  {
    char*& varname = listOfVars[nov];
    if ( nov >0) 
    {
      strcat(varStringOut,",");      
      strcat(varStringOut,varname);
    }
    else strcpy(varStringOut,varname);
  }
  listOfVars.clear();
  strcat(varStringOut,",donorVol");

  /*---------------------------------------------------*/
  E_Float* xr = coordX->begin();
  E_Float* yr = coordY->begin();
  E_Float* zr = coordZ->begin();
  E_Int nbInterpPts = coordX->getSize();
  void* a2 = &imd;
  void* a3 = &jmd;
  void* a4 = &kmd;
  void* a5 = NULL;

  E_Int type, noblk;
  E_Float volD;

  //-----------------------//
  // Champs a interpoler   //
  //-----------------------//
  // LE VOLUME DONNEUR EST LE DERNIER CHAMP DE interpolatedFields
  // N EST PAS INTERPOLE 
  E_Int nfldD = posvars0.size()+1;// variables a transferer + le volume donneur 
  PyObject* tpl = K_ARRAY::buildArray(nfldD, varStringOut, nbInterpPts, 1, 1);
  E_Float* foutp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF interpolatedFields(nbInterpPts,nfldD, foutp, true);
  interpolatedFields.setAllValuesAtNull();
  E_Float* ptrVol = interpolatedFields.begin(nfldD);

#pragma omp parallel default(shared)  private(type, noblk, volD) if (nbInterpPts > 50)
{
  FldArrayI indi(nindi); FldArrayF cf(ncfmax);
  #pragma omp for
  for (E_Int noind = 0; noind < nbInterpPts; noind++)
  {
    E_Float x = xr[noind];
    E_Float y = yr[noind];
    E_Float z = zr[noind];
    volD = 0.;
    short ok = K_INTERP::getInterpolationCell(x, y, z, adt, &donorFields,
                                              a2, a3, a4, a5, posxd, posyd, poszd, poscd,
                                              volD, indi, cf, type, noblk, interpType, 
                                              nature, penalty);   
    if (ok != 1)
    {
      ok = K_INTERP::getExtrapolationCell(x, y, z, adt, &donorFields,
                                          a2, a3, a4, a5, posxd, posyd, poszd, poscd,
                                          volD, indi, cf, type, noblk, interpType, 
                                          nature, penalty, constraint, extrapOrder);
      if (noblk > 0) //extrapolated
        volD += penaltyExtrap;
      else //orphan
        volD = penaltyOrphan;
    }
    if (ok == 1)//formule d interpolation
    {
      // on n interpole que les variables dans posvars0, qui sont les 1eres variables ds interpolatedFields
      K_INTERP::compInterpolatedValues(indi.begin(), cf, donorFields, a2, a3, a4, 
                                       noind, type, interpolatedFields, posvars0);      
    }
    ptrVol[noind] = volD;
  }
  }
  RELEASEDATA;
  delete [] varStringOut;
  return tpl; 
}
