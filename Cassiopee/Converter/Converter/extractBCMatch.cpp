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
#include "converter.h"

#define signVar(a) (a < 0 ? -1 : 1)

using namespace std;
using namespace K_FLD;

#include <map>
//=============================================================================
PyObject* K_CONVERTER::extractBCMatchNG(PyObject* self, PyObject* args )
{
  // Return index of boundary faces in receiver zone and associated fields 
  // extracted from donor zone
  
  PyObject *zone, *pyIndices, *pyVariables;
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  
  if (!PYPARSETUPLE_(args, OOO_ SSS_, &zone, &pyIndices, &pyVariables, 
               &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters )) 
     return NULL;
  
  // Zone
  // ~~~~
  E_Int ni, nj, nk, cnSize, cnNfld ; 
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  vector<PyArrayObject*> hook;

  E_Int zoneType = K_PYTREE::getFromZone(zone, 0, 1, varString, fields, locs, ni, nj, nk, 
                                         cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
                                         FlowSolutionNodes, FlowSolutionCenters);

  E_Int nfld0 = fields.size();
  if (nfld0 == 0) 
  {
    RELEASESHAREDZ(hook, varString, eltType);
    PyErr_SetString(PyExc_TypeError,
                    "extractBCMatchNG: no field to perform computation.");
    return NULL;
  }
    
  if (zoneType == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "extractBCMatchNG: not a valid zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }

  // Parent Elements 
  // ~~~~~~~~~~~~~~~
  E_Int* PE = NULL;
  if (zoneType == 2)
  {
    if (cn.size() < 3)//PE does not exist
    {
      PyErr_SetString(PyExc_TypeError, "extractBCMatchNG: ParentElements node must be defined in zone.");
      RELEASESHAREDZ(hook, varString, eltType);
      return NULL;  
    }
    else
    {
      PE = cn[2];
    }
  }
  
  // Positions des variables a extraire 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vector<E_Int> posvars;
  E_Int posvar;   
  char* varStringOut = new char[K_ARRAY::VARSTRINGLENGTH];
  varStringOut[0] = '\0';
    
  if (PyList_Check(pyVariables) != 0)
  {
    int nvariables = PyList_Size(pyVariables);
    if (nvariables > 0)
    {
      for (int i = 0; i < nvariables; i++)
      {
        PyObject* tpl0 = PyList_GetItem(pyVariables, i);
        if (PyString_Check(tpl0)) 
        {
           char* varname = PyString_AsString(tpl0); 
           if (varStringOut[0] == '\0' ) strcpy(varStringOut, varname);
           else
           {
             strcat(varStringOut, ","); strcat(varStringOut, varname);
           }
           posvar = K_ARRAY::isNamePresent(varname, varString);  
           if (posvar != -1 ) posvars.push_back(posvar);
        }
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(tpl0))
        {
          const char* varname = PyUnicode_AsUTF8(tpl0);
          if (varStringOut[0] == '\0' ) strcpy(varStringOut, varname);
          else
          {
             strcat(varStringOut, ","); strcat(varStringOut, varname);
          }
          posvar = K_ARRAY::isNamePresent(varname, varString);
          if (posvar != -1 ) posvars.push_back(posvar);
        }
#endif
        else
        {
          PyErr_Warn(PyExc_Warning, "extractBCMatchNG: variable must be a string. Skipped.");
        }
      }
    }
  }
  
  // Indices des faces 
  // ~~~~~~~~~~~~~~~~~
  FldArrayI* ind;
  E_Int res = K_NUMPY::getFromPointList(pyIndices, ind);

  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, "extractBCMatchNG: not a valid numpy for indices.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;   
  }

  E_Int* ptrInd = ind->begin();

  // Tableau des champs 
  // ~~~~~~~~~~~~~~~~~~
  int nfld = PyList_Size(pyVariables);
  int nint = ind->getSize();
  PyObject* pyFldD = K_ARRAY::buildArray3(nfld, varStringOut, nint, 1, 1, 3); 

  delete [] varStringOut;

  FldArrayF* fldD; FldArrayI* cn2;
  E_Int ni2, nj2, nk2;
  char* varStringTmp;
  K_ARRAY::getFromArray3(pyFldD, varStringTmp, fldD, ni2, nj2, nk2, cn2, eltType);

  // Extrapolation
  // ~~~~~~~~~~~~~
  for (E_Int novar = 0; novar < nfld; novar++)     
  {
      E_Int posv       = posvars[novar];
      E_Float* fieldV  = fields[posv];
      E_Float* ptrFldD = fldD->begin(novar+1); 

      for (E_Int noint = 0; noint < nint; noint++)
      {
        E_Int indint   = ptrInd[noint]-1;
        E_Int indcell  = PE[indint]-1;
        ptrFldD[noint] = fieldV[indcell];
      }
  }

  RELEASESHAREDZ(hook, varString, eltType);
  RELEASESHAREDS(pyFldD, fldD);
  RELEASESHAREDN(pyIndices, ind);

  return pyFldD; 
}


//=============================================================================
PyObject* K_CONVERTER::extractBCMatchStruct(PyObject* self, PyObject* args )
{
  // Return index of boundary faces in receiver zone and associated fields 
  // extracted from donor zone

  PyObject *fields;

  E_Int niD, njD, nkD;       // dim zone donneuse
  E_Int niR, njR, nkR;       // dim zone receveuse 

  E_Int iminD, jminD, kminD; // indices fenetre donneuse
  E_Int imaxD, jmaxD, kmaxD; // indices fenetre donneuse
  E_Int iminR, jminR, kminR; // indices fenetre receveuse
  E_Int imaxR, jmaxR, kmaxR; // indices fenetre receveuse 

  E_Int triI, triJ, triK;   // transform (issu du GC de la zone "receveuse")

  if (!PYPARSETUPLE_(args, O_ "(" IIII_ II_ ")(" IIII_ II_ ")" TIII_ TIII_,
                     &fields, &iminD, &jminD, &kminD, &imaxD, &jmaxD, &kmaxD,
                     &iminR, &jminR, &kminR, &imaxR, &jmaxR, &kmaxR, 
                     &niR, &njR, &nkR, 
                     &triI, &triJ, &triK )) return NULL;

 
  // Check array
  // ===========
  FldArrayF* FCenter; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(fields, varString, FCenter, niD, njD, nkD, 
                                     cn, eltType); 

  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError, "extractBCMatchStruct: array must be structured."); 
    if (res == 2) RELEASESHAREDS(fields, FCenter);
    return NULL; 
  }


  // 
  E_Int dim       = 3;
  E_Int noindint  = 0;
  E_Int ind ;
  //E_Int nbIntID   = (niD+1)*njD*K_FUNC::E_max(1,nkD) ;
  //E_Int nbIntJD   = (njD+1)*niD*K_FUNC::E_max(1,nkD) ;
  E_Int nbIntIR   = (niR+1)*njR*K_FUNC::E_max(1,nkR) ;
  E_Int nbIntJR   = (njR+1)*niR*K_FUNC::E_max(1,nkR) ;

  E_Int ifaceR ;
  E_Int jfaceR ;
  E_Int kfaceR ;
  E_Int shift  ;

  // compute dim 
  // Check match indice before looking for dimension (necessary when using 'addMXZones+depth=1')
  // if ((iminD == imaxD) and (njD == 1) or (nkD ==1))
  // {
  //   dim = 2; 
  // }

  // if ((jminD == jmaxD) and (niD == 1) or (nkD ==1))
  // {
  //   dim = 2; 
  // }

  // if ((kminD == kmaxD) and (niD == 1) or (njD ==1))
  // {
  //   dim = 2; 
  // }

  if ( (niR == 0) or (njR == 0) or (nkR == 0) )
  {
    dim = 2 ;
  }

  // build output arrays 
  // ===================
  E_Int nfld = FCenter->getNfld();
  E_Int nint = max(E_Int(1),(imaxD-iminD))*max(E_Int(1),(jmaxD-jminD))*max(E_Int(1),(kmaxD-kminD)); 

  // 1. tableau des indices 
  // ~~~~~~~~~~~~~~~~~~~~~~
  // 1. a Indices des faces de la zone receveuse 
  // -------------------------------------------
  PyObject* indFaceR = K_NUMPY::buildNumpyArray(nint,1,1);
  E_Int* ptrIndFaceR = K_NUMPY::getNumpyPtrI(indFaceR);


  // 2. tableau des champs 
  // ~~~~~~~~~~~~~~~~~~~~~
  PyObject* pyFldD = K_ARRAY::buildArray3(nfld, varString, nint, 1, 1, 3); 

  FldArrayF*  fBC; 
  FldArrayI* cnBC;
  E_Int ni2, nj2, nk2;
  K_ARRAY::getFromArray3(pyFldD, varString, fBC, ni2, nj2, nk2, cnBC, eltType);
  
  // Cas 2D
  // ======
  if (dim == 2)
  {
    // Face donneuse en i 
    // ******************
    if (iminD == imaxD) 
    { 
      // printf("Frontiere en i \n");
      // 1. tableau des indices
      // ----------------------

      // 1.a. face receveuse en i
      if (iminR == imaxR) 
      {
        noindint = 0;

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
        { 
          if (triJ > 0)
      {
        jfaceR =  jface + jminR - jminD ;
      }
      else
      {
        jfaceR = (jmaxR-2) - (jface-jminD+1); 
      }

          ptrIndFaceR[noindint] = iminR - 1 + jfaceR*(niR+1) ;          noindint++;
    }
      }
      // 1.b. face receveuse en j
      else if (jminR == jmaxR) 
      {
    noindint = 0;

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
        { 
          if (triI > 0)
      {
        ifaceR = jface + iminR - jminD ;
      }
      else
      {
        ifaceR = (imaxR-2) - (jface-jminD+1)  ; 
      }

          shift  = (jminR-1)*niR + nbIntIR ; 
          ptrIndFaceR[noindint] = shift + ifaceR ;
          noindint++;
    }
      } // jminR=jmaxR

      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0 ;

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
        {
          if (iminD==1) { ind = jface*niD         ; } 
          else          { ind = jface*niD + niD-1 ; } 

          fld[noindint] = fce[ind] ; 
          noindint++;
        } 
      } // var loop 
    }
    // Si frontiere en j
    // *****************
    else if (jminD == jmaxD) 
    {      
      // 1.a. face receveuse en i
      if (iminR == imaxR) 
      {
        noindint = 0;

        for (E_Int iface = iminD - 1 ; iface < imaxD-1 ; iface ++)
        { 
          if (triJ > 0)
          {
            jfaceR =  iface + jminR - iminD ; 
          }
          else
          {
            jfaceR = (jmaxR-2) - (iface-iminD+1);
          }

          ptrIndFaceR[noindint] = iminR - 1 + jfaceR*(niR+1) ;
          noindint++;
        }
      } // iminR==imaxR 

      // 1.b. face receveuse en j
      else if (jminR == jmaxR) 
      {
        // printf("Frontiere receveuse en j \n");
        noindint = 0;

        for (E_Int iface = iminD - 1 ; iface < imaxD-1 ; iface ++)
        { 
          E_Int shift  = (jminR-1)*niR + nbIntIR ;

          if (triI > 0)
          {
            ifaceR = iface + iminR - iminD ;
          }
          else
          {
            ifaceR = (imaxR-2) - (iface-iminD+1) ;
          }

          ptrIndFaceR[noindint] = shift + ifaceR ;
          // printf(" %d ", ptrIndFaceR[noindint]);
          noindint++;
        }
      }
  
      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0 ;

        for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
        {
          if (jminD==1) { ind = iface               ; } 
          else          { ind = iface + niD*(njD-1) ; } 

          fld[noindint] = fce[ind] ; 
          noindint++;
        }
      } // var loop
    }// (si frontiere en j)
 
  } //(si dim=2)

  // Cas 3D
  // ======
  else if (dim == 3)
  {
    // ********************
    // Frontiere donneuse i 
    // ********************
    if (iminD == imaxD)
    {
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      // Frontiere receveuse en i 
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (iminR == imaxR)
      {
        noindint = 0 ;

        if (abs(triJ)==2) // kD <-> kR et  jD <-> jR
        {
          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triK > 0) { kfaceR = kface + kminR - kminD ;     }
            else          { kfaceR = (kmaxR-2)-(kface-kminD+1) ; }

            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
            {
              if (triJ > 0){ jfaceR = jface + jminR - jminD ;    }
              else         { jfaceR = (jmaxR-2)-(jface-jminD+1); }

              ptrIndFaceR[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              // printf("%d ", ptrIndFaceR[noindint]);
              noindint++;
            }
          }
        } // triJ=2

        else if (abs(triJ)==3) // kD <-> jR et  jD <-> kR
        {
          // printf("Face receveuse jD=kR, kD=jR) \n");
          // printf("indR : ");

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)  
            { 
              if (triK > 0) { kfaceR = jface + kminR - jminD ;     }
              else          { kfaceR = (kmaxR-2)-(jface-jminD+1) ; }

              if (triJ > 0) { jfaceR = kface + jminR - kminD ;     }
              else          { jfaceR = (jmaxR-2)-(kface-kminD+1);  }
      
              ptrIndFaceR[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              noindint++;
            }
          } 
        } // triJ=3
      }
     
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      // Frontiere receveuse en j 
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      else if (jminR == jmaxR)
      {
        noindint = 0 ;

        if (abs(triI)==2) // jD <-> iR et  kD <-> kR
        {
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triK > 0) { kfaceR = kface + kminR - kminD ;     }
            else          { kfaceR = (kmaxR-2)-(kface-kminD+1) ; }

            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
            {
              if (triI > 0){ ifaceR = jface + iminR - jminD ;    }
              else         { ifaceR = (imaxR-2)-(jface-jminD+1); }

              ptrIndFaceR[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              noindint++;
            }
          } 
        } // triI=1

        if (abs(triI)==3) // jD <-> kR et  kD <-> iR
        {
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)  
            { 
              if (triK > 0) { kfaceR = jface + kminR - jminD ;     }
              else          { kfaceR = (kmaxR-2)-(jface-jminD+1) ; }

              if (triI > 0) { ifaceR = kface + iminR - kminD ;     }
              else          { ifaceR = (imaxR-2)-(kface-kminD+1);  }
      
              ptrIndFaceR[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              noindint++;
            }
          }
        } // triJ=3
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~
      // Frontiere receveuse en k
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      else if (kminR == kmaxR)
      {
        noindint = 0 ;

        if (abs(triI)==2) // jD <-> iR et  kD <-> jR
        {
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triJ > 0) { jfaceR = kface + jminR - kminD ;     }
            else          { jfaceR = (jmaxR-2)-(kface-kminD+1) ; }

            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
            {
              if (triI > 0){ ifaceR = jface + iminR - jminD ;    }
              else         { ifaceR = (imaxR-2)-(jface-jminD+1); }

              ptrIndFaceR[noindint] = shift + ifaceR + jfaceR*niR ;
              noindint++;
            }
          } 
        } // triI=1

        else if (abs(triI)==3) // jD <-> jR et  kD <-> iR
        {
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            if (triI > 0) { ifaceR = kface + iminR - kminD ;     }
            else          { ifaceR = (imaxR-2)-(kface-kminD+1);  }

            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)  
            { 
              if (triJ > 0) { jfaceR = jface + jminR - jminD ;     }
              else          { jfaceR = (jmaxR-2)-(jface-jminD+1) ; }
      
            ptrIndFaceR[noindint] = shift + ifaceR + jfaceR*niR ;
            noindint++;
            }
          } 
        } // triJ=3
      }

      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0;

        for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
        {
          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
          {
            if (iminD==1) { ind = iminD-1 + jface*niD + kface*njD*niD; } 
            else          { ind = niD-1    + jface*niD + kface*njD*niD; } 

            fld[noindint] = fce[ind]; 
            noindint++;
          }
        }

      } // var loop 
    }
    // Si frontiere en j
    // *****************
    else if (jminD == jmaxD)
    {

    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en i 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (iminR == imaxR)
      {
        noindint = 0 ;

        if (abs(triJ)==1) // kD <-> kR et  iD <-> jR
        {
          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triK > 0) { kfaceR = kface + kminR - kminD ;     }
            else          { kfaceR = (kmaxR-2)-(kface-kminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triJ > 0){ jfaceR = iface + jminR - iminD ;    }
              else         { jfaceR = (jmaxR-2)-(iface-iminD+1); }

              ptrIndFaceR[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              noindint++;
            }
          } 
        } 

        else if (abs(triJ)==3) // kD <-> jR et  iD <-> kR
        {
          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triK > 0) { kfaceR = iface + kminR - iminD ;     }
              else          { kfaceR = (kmaxR-2)-(iface-iminD+1) ; }

              if (triJ > 0) { jfaceR = kface + jminR - kminD ;     }
              else          { jfaceR = (jmaxR-2)-(kface-kminD+1);  }
      
            ptrIndFaceR[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
            noindint++;
            }
          } 
        } // triJ=3
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~
      // Frontiere receveuse en j 
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      else if (jminR == jmaxR)
      {
        noindint = 0 ;

        if (abs(triI)==1) // iD <-> iR et  kD <-> kR
        {
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triK > 0) { kfaceR = kface + kminR - kminD ;     }
            else          { kfaceR = (kmaxR-2)-(kface-kminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triI > 0){ ifaceR = iface + iminR - iminD ;    }
              else         { ifaceR = (imaxR-2)-(iface-iminD+1); }

              ptrIndFaceR[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              // printf("%d ", ptrIndFaceR[noindint]);
              noindint++;
            }
          } 
        } // triI=1

        else if (abs(triI)==3) // iD <-> kR et  kD <-> iR
        {
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triK > 0) { kfaceR = iface + kminR - iminD ;     }
              else          { kfaceR = (kmaxR-2)-(iface-iminD+1) ; }

              if (triI > 0) { ifaceR = kface + iminR - kminD ;     }
              else          { ifaceR = (imaxR-2)-(kface-kminD+1);  }
      
            ptrIndFaceR[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
            noindint++;
            }
          } 
        } // triJ=3 
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~
      // Frontiere receveuse en k 
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      else if (kminR == kmaxR)
      {
        noindint = 0 ;

        if (abs(triI)==1) // iD <-> iR et  kD <-> jR
        {
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triJ > 0) { jfaceR = kface + jminR - kminD ;     }
            else          { jfaceR = (jmaxR-2)-(kface-kminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triI > 0){ ifaceR = iface + iminR - iminD ;    }
              else         { ifaceR = (imaxR-2)-(iface-iminD+1); }

              ptrIndFaceR[noindint] = shift + ifaceR + jfaceR*niR ;
              // printf("%d ", ptrIndFaceR[noindint]);
              noindint++;
            }
          } 
        } // triI=1

        else if (abs(triI)==3) // iD <-> jR et  kD <-> iR
        {
          // printf("Face receveuse iD=jR, kD=iR) \n");
          // printf("indR : ");
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            if (triI > 0) { ifaceR = kface + iminR - kminD ;     }
            else          { ifaceR = (imaxR-2)-(kface-kminD+1);  }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triJ > 0) { jfaceR = iface + jminR - iminD ;     }
              else          { jfaceR = (jmaxR-2)-(iface-iminD+1) ; }
      
              ptrIndFaceR[noindint] = shift + ifaceR + jfaceR*niR ;
              // printf("%d ", ptrIndFaceR[noindint]);
              noindint++;
            }
          } 
        } // triJ=3
      }

      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0 ;

        for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
        {
          for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
          {
            if (jminD==1) { ind = iface + kface*njD*niD            ; } 
            else         { ind = iface + kface*njD*niD + (njD-1)*niD; } 

            fld[noindint] = fce[ind] ; 
            noindint++;
          }
        }

      } // var loop
    }
    // Si frontiere en k
    // *****************
    else if (kminD == kmaxD)
    {

      // ~~~~~~~~~~~~~~~~~~~~~~~~
      // Frontiere receveuse en i 
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (iminR == imaxR)
      {
        noindint = 0 ;

        if (abs(triJ)==2) // iD <-> kR et  jD <-> jR
        {
          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
          {
            if (triJ > 0){ jfaceR = jface + jminR - jminD ;    }
            else         { jfaceR = (jmaxR-2)-(jface-jminD+1); }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triK > 0) { kfaceR = iface + kminR - iminD ;     }
              else          { kfaceR = (kmaxR-2)-(iface-iminD+1) ; }

              ptrIndFaceR[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              noindint++;
            }
          } 
        } // triJ=2

        else if (abs(triJ)==1) // iD <-> jR et  jD <-> kR
        {
          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)
          {
            if (triK > 0) { kfaceR = jface + kminR - jminD ;     }
            else          { kfaceR = (kmaxR-2)-(jface-jminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triJ > 0) { jfaceR = iface + jminR - iminD ;     }
              else          { jfaceR = (jmaxR-2)-(iface-iminD+1);  }

              ptrIndFaceR[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              noindint++;
            }
          } 
        } // triJ=3
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~
      // Frontiere receveuse en j 
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      else if (jminR == jmaxR)
      {
        noindint = 0 ;

        if (abs(triI)==2) // iD <-> kR et  jD <-> iR
        {
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
          {
            if (triI > 0){ ifaceR = jface + iminR - jminD ;    }
            else         { ifaceR = (imaxR-2)-(jface-jminD+1); }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triK > 0) { kfaceR = iface + kminR - iminD ;     }
              else          { kfaceR = (kmaxR-2)-(iface-iminD+1) ; }

              ptrIndFaceR[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              noindint++;
            }
          } 
        } // triJ=2

        else if (abs(triI)==1) // iD <-> iR et  jD <-> kR
        {
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)
          {
            if (triK > 0) { kfaceR = jface + kminR - jminD ;     }
            else          { kfaceR = (kmaxR-2)-(jface-jminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triI > 0) { ifaceR = iface + iminR - iminD ;     }
              else          { ifaceR = (imaxR-2)-(iface-iminD+1);  }

              ptrIndFaceR[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              noindint++;
            }
          } 
        } // triJ=3
      }
     
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      // Frontiere receveuse en k
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      else if (kminR == kmaxR)
      {
        noindint = 0 ;

        if (abs(triI)==1) // iD <-> iR et  jD <-> jR
        {
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
          {
            if (triJ > 0){ jfaceR = jface + jminR - jminD ;    }
            else         { jfaceR = (jmaxR-2)-(jface-jminD+1); }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triI > 0) { ifaceR = iface + iminR - iminD ;     }
              else          { ifaceR = (imaxR-2)-(iface-iminD+1) ; }

              ptrIndFaceR[noindint] = shift + ifaceR + jfaceR*niR ;
              noindint++;
            }
          } 
        } // triI=1

        else if (abs(triI)==2) // iD <-> jR et  jD <-> iR
        {
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)
          {
            if (triI > 0) { ifaceR = jface + iminR - jminD ;     }
            else          { ifaceR = (imaxR-2)-(jface-jminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triJ > 0) { jfaceR = iface + jminR - iminD ;     }
              else          { jfaceR = (jmaxR-2)-(iface-iminD+1);  }

              ptrIndFaceR[noindint] = shift + ifaceR + jfaceR*niR ;
              noindint++;
            }
          } 
        } // triJ=3
      }

      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0;

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
        {
          for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
          {
            if (kminD==1) { ind = iface + jface*niD            ; } 
            else          { ind = iface + jface*niD + (nkD-1)*niD*njD; } 

            fld[noindint] = fce[ind]; 
            noindint++;
          }
        }
      } // var loop
    }
  }

  PyObject* tplOut;
  tplOut = Py_BuildValue("[OO]", indFaceR, pyFldD);

  RELEASESHAREDS(fields, FCenter);
  RELEASESHAREDS(pyFldD, fBC);

  Py_DECREF(indFaceR);
  Py_DECREF(pyFldD);

  return tplOut; 
}

//=============================================================================
//=============================================================================
void K_CONVERTER::indface2index(E_Int indface, E_Int ni, E_Int nj, E_Int nk, E_Int& ind)
{
  // Return cell index 'ind' given a face index 'indFace' 
  // Warning: only valid for boundary faces imin/imax, jmin/jmax, kmin/kmax 
  //          because for general case 1 face connect with 2 cells
  //          here information of 'min or max' enable to pick a unique (i,j,k) 
  
  // printf("indface : %d \n", indface);
  // printf("ni : %d, nj : %d, nk : %d \n", ni,nj,nk);

  E_Int i,j,k,res;

  E_Int nbIntI = (ni+1)*nj*K_FUNC::E_max(1,nk);
  E_Int nbIntJ = (nj+1)*ni*K_FUNC::E_max(1,nk);

  i = 10;
  j = 20;
  k = 30; 

  // printf("nbIntI : %d , nbIntJ : %d \n",nbIntI, nbIntJ); 

  if ( indface < nbIntI )
  {
    k   = indface/( (ni+1)*nj ) + 1;
    res = indface - (k-1)*(ni+1)*nj;
    j   = res/(ni+1) + 1;
    i   = res - (j-1)*(ni+1) + 1;

    if ( i==ni+1) { i = ni; }
  }
  else if ( indface < nbIntI + nbIntJ )
  {
    res = indface - nbIntI;
    k   = res/(ni*(nj+1)) + 1;
    res = indface - nbIntI - (k-1)*ni*(nj+1); 
    j   = res/ni + 1 ;
    i   = res - (j-1)*ni + 1; 

    if ( j==nj+1) { j = nj; }
  }
  else
  {
    res = indface - nbIntI - nbIntJ;
    k   = res/(ni*nj) + 1;
    res = res - (k-1)*ni*nj;
    j   = res/ni + 1;
    i   = res - (j-1)*ni + 1;

    if ( k==nk+1) { k = nk; }
    
  }

  ind = i-1 + (j-1)*ni + (k-1)*nj*ni;
  // printf("iface: %d, i: %d, j: %d, k: %d, ind: %d \n",indface,i,j,k,ind);

  return;
}

//=============================================================================
// compute fld = 0.5(fldD+flR)
//=============================================================================
PyObject* K_CONVERTER::buildBCMatchFieldNG(PyObject* self, PyObject* args )
{

  PyObject *zone, *pyIndR, *pyFldD, *pyVariables ; 
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  
  if (!PYPARSETUPLE_(args, OOOO_ SSS_, &zone, &pyIndR, &pyFldD, 
                     &pyVariables, &GridCoordinates, &FlowSolutionNodes, 
                     &FlowSolutionCenters )) return NULL;

  // Zone
  // ~~~~
  E_Int ni, nj, nk, cnSize, cnNfld ; 
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  vector<PyArrayObject*> hook;

  E_Int zoneType = K_PYTREE::getFromZone(zone, 0, 1, varString, fields, locs, ni, nj, nk, 
                                         cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
                                         FlowSolutionNodes, FlowSolutionCenters);

  if (zoneType == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldNG: not a valid zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }

  // Parent Elements 
  // ~~~~~~~~~~~~~~~
  E_Int* PE = NULL;
  if (zoneType == 2)
  {
    if (cn.size() < 3) //PE does not exist
    {
      PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldNG: ParentElements node must be defined in zone.");
      RELEASESHAREDZ(hook, varString, eltType);
      return NULL;  
    }
    else PE = cn[2];
  }
  
  // Champs de la zone donneuse
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~
  E_Int ni2, nj2, nk2;
  FldArrayF* fldD;
  FldArrayI* cn2;
  char* varStringOut;
  E_Int res2 = K_ARRAY::getFromArray3(pyFldD, varStringOut, fldD, ni2, nj2, nk2, 
                                      cn2, eltType); 

  if (res2 != 1)
  {
    PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldNG: wrong array."); 
    RELEASESHAREDS(pyFldD, fldD);
    return NULL; 
  }

  // Positions des variables a extraire 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vector<E_Int> posvars;
  E_Int posvar;   

  if (PyList_Check(pyVariables) != 0)
  {
    int nvariables = PyList_Size(pyVariables);
    if (nvariables > 0)
    {
      for (int i = 0; i < nvariables; i++)
      {
        PyObject* tpl0 = PyList_GetItem(pyVariables, i);
        if (PyString_Check(tpl0)) 
        {
          char* varname = PyString_AsString(tpl0); 
          // Verif. presence variables a extraire dans le dict.
          E_Int verif = K_ARRAY::isNamePresent(varname, varStringOut);  
          if (verif == -1) 
          {
            PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldNG: Variable not found in dictionary allMatch.");
          }
          posvar = K_ARRAY::isNamePresent(varname, varString);  
          if (posvar != -1 ) posvars.push_back(posvar);
        }
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(tpl0))
        {
          const char* varname = PyUnicode_AsUTF8(tpl0); 
          // Verif. presence variables a extraire dans le dict.
          E_Int verif = K_ARRAY::isNamePresent(varname, varStringOut);  
          if (verif == -1) 
          {
            PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldNG: Variable not found in dictionary allMatch.");
          }
          posvar = K_ARRAY::isNamePresent(varname, varString);  
          if (posvar != -1 ) posvars.push_back(posvar); 
        }
#endif
        else
        {
          PyErr_Warn(PyExc_Warning, "buildBCMatchFieldNG: variable must be a string. Skipped.");
        }  
      }
    }
  }

  // Indices des faces 
  // ~~~~~~~~~~~~~~~~~
  FldArrayI* indR;
  E_Int res = K_NUMPY::getFromPointList(pyIndR, indR);

  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldNG: not a valid numpy for indR.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;   
  }

  // Tableau des champs (output)
  // ~~~~~~~~~~~~~~~~~~
  E_Int nfld = fldD->getNfld();
  E_Int nind = indR->getSize();
  PyObject* pyFld = K_ARRAY::buildArray3(nfld, varStringOut, nind, 1, 1, 3); 

  FldArrayF* fld; FldArrayI* cn3;
  char* varStringTmp;
  K_ARRAY::getFromArray3(pyFld, varStringTmp, fld, ni2, nj2, nk2, cn3, eltType);


  // Build 0.5(fldD+fldR) array on boundary faces
  // ============================================
  E_Int* ptrIndR = indR->begin();

  for (E_Int var = 1; var <= nfld; var++)
  {
    E_Int posv          = posvars[var-1];
    E_Float* fieldV     = fields[posv];
    E_Float* ptrFldD    = fldD->begin(var);
    E_Float* ptrFld     = fld->begin(var);

    for (E_Int noindint = 0 ; noindint < nind ; noindint++)
    {
      E_Int indFace    = ptrIndR[noindint]-1;
      E_Int indcell    = PE[indFace]-1;
      ptrFld[noindint] = 0.5*( fieldV[indcell]+ptrFldD[noindint] );    
    }
  }

  RELEASESHAREDS(pyFldD, fldD);
  RELEASESHAREDN(pyIndR, indR);
  RELEASESHAREDZ(hook, varString, eltType); 
  RELEASESHAREDS(pyFld, fld);

  return pyFld;  
}


//=============================================================================
//=============================================================================
PyObject* K_CONVERTER::buildBCMatchFieldStruct(PyObject* self, PyObject* args )
{
//   // indR fldD fldR >> fld = 0.5(fldD+flR)

  PyObject *pyFieldsR, *pyIndR, *pyFldD, *pyNcnt; 
  
  if (!PYPARSETUPLE_(args, OOOO_, &pyFieldsR, &pyIndR, &pyFldD, &pyNcnt )) return NULL;

  // Get current zone fields (in volume)
  // ===================================
  E_Int ni, nj, nk ;
  FldArrayF* fieldsR; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(pyFieldsR, varString, fieldsR, ni, nj, nk, 
                                     cn, eltType); 

  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldStruct: array must be structured."); 
    if (res == 2) RELEASESHAREDS(pyFieldsR, fieldsR);
    return NULL; 
  }

  // Get donor zone fields (on BCMatch)
  // ==================================
  E_Int ni2, nj2, nk2 ;
  FldArrayF* fldD;
  E_Int res2 = K_ARRAY::getFromArray3(pyFldD, varString, fldD, ni2, nj2, nk2, 
                                      cn, eltType); 

  // printf("ni2 : %d, nj2 : %d, nk2 : %d \n", ni2, nj2, nk2);

  if (res2 != 1)
  {
    PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldStruct: array must be structured."); 
    if (res2 == 2) RELEASESHAREDS(pyFldD, fldD);
    return NULL; 
  }

  E_Int nfld = fldD->getNfld();


  // Get index of boundary faces in current zone
  // ============================================
  FldArrayI* indR;
  E_Int resi = K_NUMPY::getFromNumpyArray(pyIndR, indR);
  if ( resi == 0)
  {
     PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldStruct: not a valid numpy for indices of BC (indR).");
     RELEASESHAREDN(pyIndR, indR);
     return NULL;   
  }

  E_Int  nind    = indR->getSize();

  // Get ncount array if supplied (used for near-match or TNC match)
  // ================================================================
  E_Bool needcount = true;
  FldArrayI* ncount = NULL;
    
  if (pyNcnt != Py_None)
  {
    E_Int resi = K_NUMPY::getFromNumpyArray(pyNcnt, ncount);
    if (resi == 0)
    {
      PyErr_SetString(PyExc_TypeError, "buildBCMatchFieldStruct: not a valid numpy for ncount array.");
      RELEASESHAREDN(pyNcnt, ncount);
      return NULL;   
    }
    
    // E_Int* ptrNcnt = ncount->begin();

    // for (E_Int ko=0 ; ko< ncount->getSize(); ko++)
    // std::cout << "ncout[" << ko << "]= " << ptrNcnt[ko] << std::endl;
  }
  else
  {
    needcount = false;
  }

  // Create output array 
  // ===================
  E_Int nn;
  PyObject* pyFld = K_ARRAY::buildArray3(nfld, varString, nind, 1, 1, 3); 
  // printf("nfld : %d, nind : %d \n",nfld,nind);
  FldArrayF* fld;
  K_ARRAY::getFromArray3(pyFld, varString, fld, nind, nn, nn, cn, eltType);


  // Build 0.5(fldD+fldR) array on boundary faces
  // ============================================
  E_Int  ind,indFace;
  E_Int* ptrIndR = indR->begin();

  E_Int* ptrNcnt = NULL;
  if (needcount) ptrNcnt = ncount->begin();
    
  for (E_Int noindint = 0 ; noindint < nind ; noindint++)
  {
    indFace = ptrIndR[noindint]; 
    indface2index(indFace,ni,nj,nk,ind);
    // printf("indFace: %d, ind: %d \n", indFace,ind);

    for (E_Int var = 1; var <= nfld; var++)
    {
      E_Float* ptrFieldsR = fieldsR->begin(var);
      E_Float* ptrFldD    = fldD->begin(var);
      E_Float* ptrFld     = fld->begin(var);
      if (needcount)
      {
        std::cout << "ncout[" << noindint << "]= " << ptrNcnt[noindint] << std::endl;
        ptrFld[noindint]  = 0.5*( ptrFieldsR[ind]/ptrNcnt[noindint]+ptrFldD[noindint] );   
      }
      else
      {
        ptrFld[noindint]  = 0.5*( ptrFieldsR[ind]+ptrFldD[noindint] );
      }
    }
  }

  RELEASESHAREDN(pyIndR, indR);
  RELEASESHAREDS(pyFldD, fldD);
  RELEASESHAREDS(pyFld , fld );
  RELEASESHAREDS(pyFieldsR, fieldsR);
  if (needcount) RELEASESHAREDN(pyNcnt, ncount);

  return pyFld;
}
