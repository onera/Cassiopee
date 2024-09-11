/*    
    Copyright 2013-2024 Onera.

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
# include <stdio.h>
#include <map>
# include "connector.h"
using namespace std;
using namespace K_FLD;

#define BLOCK1\
  for (size_t v2 =0; v2 < pyRcvIndices.size(); v2++)\
    RELEASESHAREDN(pyRcvIndices[v2], FldRcvIndices[v2]);\

#define BLOCK2\
  BLOCK1;  \
  for (size_t v2 =0; v2 < pyDnrIndices.size(); v2++)\
    RELEASESHAREDN(pyDnrIndices[v2], FldDnrIndices[v2]);     \

#define BLOCK3\
  BLOCK2; \
  for (size_t v2 =0; v2 < pyEXIndir.size(); v2++)\
    RELEASESHAREDN(pyEXIndir[v2], FldEXIndir[v2]);\

# define BLOCK4 \
  BLOCK3; \
  for (size_t v2 = 0; v2 < pyDnrCoefs.size(); v2++)\
    RELEASESHAREDN(pyDnrCoefs[v2], FldDnrCoefs[v2]);\

# define BLOCK5\
  BLOCK4;\
  for (size_t v2 =0; v2 < pyTypes.size(); v2++)\
    RELEASESHAREDN(pyTypes[v2], FldTypes[v2]); \

#define BLOCK6\
  BLOCK5; \
  for (size_t v2 = 0; v2 < pyCellN.size(); v2++)\
    RELEASESHAREDS(pyCellN[v2], FldCellN[v2]);\

//=============================================================================
/* Ecriture des coefficients d'interpolation dans un fichier relisible 
   par elsA */
//=============================================================================
PyObject* K_CONNECTOR::writeCoefs(PyObject* self, PyObject* args)
{
  PyObject *RcvIndexMap, *DonorIndexMap;  // indices des pts interpoles et des cellules donneuses
  PyObject *EXDirectionMap;
  PyObject *DonorInterpCoefMap;           // coefficients d interpolation
  PyObject *DonorCellNMap;                // champs cellN des zones d interpolations
  PyObject *InterpTypesMap;               // types des interpolations
  PyObject *ZoneDimMap;                   // dimensions des grilles aux centres des cellules pour chaque zone
  PyObject *BlockRcvIdMap;                // Id des zones interpolees
  PyObject *Npts;                         // nombre total de points interpoles (dans tous les blocs) par ne zone d interpolation
  char* PrefixFile;                       // prefixe pour le nommage des fichiers elsA
  E_Int isEX = 0;                      // isEX = 0 ou 1 : Chimere avec 2 ou 1 rangees de cellules fictives
  E_Int NZones;                        // nombre de points interpoles
  E_Int Solver;                        // solveur pour lequel les fichiers sont ecrits
  E_Int NGhostCells;
  if (!PYPARSETUPLE_(args, I_ OOOO_ OOOO_ O_ S_ III_,
                    &NZones, &BlockRcvIdMap, &RcvIndexMap, 
                    &EXDirectionMap, &DonorIndexMap, &DonorInterpCoefMap, &InterpTypesMap,
                    &DonorCellNMap, &ZoneDimMap,
                    &Npts, &PrefixFile, &isEX, &Solver, &NGhostCells))
  {
    return NULL;
  }

  /*-------------------*/
  /* Variables locales */
  /*-------------------*/
  E_Int intKey;

  /*-------------------------------*/
  /* Extraction du solveur */
  /*-------------------------------*/
  E_Int solver = Solver; //1: elsA, 2: Cassiopee

  /*-------------------------------*/
  /* Extraction du nombre de zones */
  /*-------------------------------*/
  E_Int nzones = NZones;

  /*-------------------------------*/
  /* Nb de ghost cells */
  /*-------------------------------*/
  E_Int nbOfGc = NGhostCells;

  /*-------------------------------------------------------------*/
  /* Extraction du nombre de points interpoles par zone donneuse */
  /*-------------------------------------------------------------*/
  vector<E_Int> nbInterpCells;
  PyObject* tpl;
  E_Int ind;
  E_Int size = PyList_Size(Npts);
  for (E_Int v  = 0 ; v < size; v++)
  {
    tpl = PyList_GetItem(Npts, v);
    if (PyInt_Check(tpl) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable " SF_D_ ". Skipped...\n", v);
    }
    else
    {
      ind = PyInt_AsLong(tpl);
      nbInterpCells.push_back(ind);
    }
  }

  /*------------------------------------------------*/
  /* Extraction des listes d'Id des blocs receveurs */
  /*------------------------------------------------*/
  map<E_Int, vector<E_Int> > blockRcvIdMap;
  PyObject* pyListValues;
  if (PyDict_Check (BlockRcvIdMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: BlockRcvIdMap must be a dictionary.");
    return NULL;
  }
  PyObject * key_list = PyDict_Keys(BlockRcvIdMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable " SF_D_ ". Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int>
    vector<E_Int> tmpListValue;
    pyListValues = PyDict_GetItem(BlockRcvIdMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: BlockRcvIdMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* pyIntValue = PyList_GetItem(pyListValues, v);
      if (PyInt_Check(pyIntValue) == 0)
      {
        printf("Warning: writeCoefs: invalid int value in  BlockRcvIdMap\n");
      }
      else
      {
        E_Int intValue = PyInt_AsLong(pyIntValue);
        tmpListValue.push_back(intValue);
      }
    }
   blockRcvIdMap[intKey] = tmpListValue;
  }
  /*--------------------------------------------------------*/
  /* Extraction des listes de tableaux d indices            */
  /*--------------------------------------------------------*/
  map<E_Int,vector<E_Int*> > rcvIndexMap;
  map<E_Int,vector<E_Int*> > donorIndexMap;
  vector<PyObject*> pyRcvIndices;
  vector<FldArrayI*> FldRcvIndices;
  vector<PyObject*> pyDnrIndices;
  vector<FldArrayI*> FldDnrIndices;
  // rcvIndexMap : indices du bloc interpole
  if (PyDict_Check (RcvIndexMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: RcvIndexMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(RcvIndexMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable " SF_D_ ". Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int*>
    vector<E_Int*> tmpListValue;
    pyListValues = PyDict_GetItem(RcvIndexMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: RcvIndexMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* intArray = PyList_GetItem(pyListValues, v);
      pyRcvIndices.push_back(intArray);

      FldArrayI* indArrayI;
      E_Int ok = K_NUMPY::getFromNumpyArray(intArray, indArrayI, true);
      if ( ok == 1)
      {
        E_Int* indArray = indArrayI->begin();
        FldRcvIndices.push_back(indArrayI);
        tmpListValue.push_back(indArray);
      }
      else 
      {
        BLOCK1;
        PyErr_SetString(PyExc_TypeError,"writeCoefs: invalid array in RcvIndexMap");
        return NULL;
      }
    }
    rcvIndexMap[intKey] = tmpListValue;
  }

  // DonorIndexArray  : indices du bloc donneur
  if (PyDict_Check (DonorIndexMap) == 0)
  {
    BLOCK1;
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: DonorIndexMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(DonorIndexMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable " SF_D_ ". Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int*>
    vector<E_Int*> tmpListValue;
    pyListValues = PyDict_GetItem(DonorIndexMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      BLOCK1;
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: DonorIndexMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* intArray = PyList_GetItem(pyListValues, v);
      pyDnrIndices.push_back(intArray);
      FldArrayI* indArrayI;
      E_Int ok = K_NUMPY::getFromNumpyArray(intArray, indArrayI, true);
      if ( ok == 1)
      {
        E_Int* indArray = indArrayI->begin();
        FldDnrIndices.push_back(indArrayI);
        tmpListValue.push_back(indArray);
      }
      else 
      {
        BLOCK2;
        PyErr_SetString(PyExc_TypeError,"writeCoefs: invalid array in DonorIndexMap");
        return NULL;
      }
    }
    donorIndexMap[intKey] = tmpListValue;
  }

  /*----------------------------------------------------------*/
  /* Extraction des tableaux d indirection pour les points EX */
  /*----------------------------------------------------------*/
  map<E_Int,vector<E_Int*> >directionEXMap;
  vector<PyObject*> pyEXIndir;
  vector<FldArrayI*> FldEXIndir;
  if (isEX)
  {    
    if (PyDict_Check (EXDirectionMap) == 0)
    {
      BLOCK2;
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: EXDirectionMap must be a dictionary.");
      return NULL;
    }
    key_list = PyDict_Keys(EXDirectionMap);
    size = PyList_Size(key_list);
    //Get keys and corresponding values from the dictionary
    for (E_Int i  = 0 ; i < size; i++)
    {
      PyObject * pyKey = 0;
      pyKey = PyList_GetItem(key_list, i);
      if (PyInt_Check(pyKey) == 0)
      {
        printf("Warning: writeCoefs: invalid int for variable " SF_D_ ". Skipped...\n", i);
      }
      else
        intKey = PyInt_AsLong(pyKey);
      //Convert to a C++ vector<E_Int*>
      vector<E_Int*> tmpListValue;
      pyListValues = PyDict_GetItem(EXDirectionMap, pyKey);
      // check if pyListValue is a list
      if (PyList_Check (pyListValues) == 0)
      {
        BLOCK2;
        PyErr_SetString(PyExc_TypeError, 
                        "writeCoefs: EXDirectionMap must contain lists.");
        return NULL;
      }
      // fill C++ map
      E_Int tmpListValueSize = PyList_Size(pyListValues);
      for (E_Int v  = 0 ; v < tmpListValueSize; v++)
      {
        PyObject* intArray = PyList_GetItem(pyListValues, v);
        pyEXIndir.push_back(intArray);
        FldArrayI* indArrayI;
        E_Int ok = K_NUMPY::getFromNumpyArray(intArray, indArrayI, true);
        if ( ok == 1)
        {
          E_Int* indArray = indArrayI->begin();
          FldEXIndir.push_back(indArrayI);
          tmpListValue.push_back(indArray);
        }
        else 
        {
          BLOCK3;
          PyErr_SetString(PyExc_TypeError,"writeCoefs: invalid array in EXIndirMap");
          return NULL;
        }
      }
      directionEXMap[intKey] = tmpListValue;
    }
  }

  /*--------------------------------------------------------*/
  /* Extraction des listes de coefficients d interpolation  */
  /*--------------------------------------------------------*/
  // DonorInterpCoef : coefficient d'interpolation
  map<E_Int,vector<FldArrayF> > donorCoefMap;
  vector<PyObject*> pyDnrCoefs;
  vector<FldArrayF*> FldDnrCoefs;
  if (PyDict_Check (DonorInterpCoefMap) == 0)
  {
    BLOCK3;
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: DonorInterpCoefMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(DonorInterpCoefMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable " SF_D_ ". Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<FldArrayF>
    vector<FldArrayF> tmpListValue;
    pyListValues = PyDict_GetItem(DonorInterpCoefMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      BLOCK3;
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: DonorInterpCoefMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* pyArrayValue = PyList_GetItem(pyListValues, v);
      pyDnrCoefs.push_back(pyArrayValue);
      FldArrayF* dnrCoefs;
      E_Int ok = K_NUMPY::getFromNumpyArray(pyArrayValue, dnrCoefs, true);
      if (ok == 1)
      {
        FldDnrCoefs.push_back(dnrCoefs);
        tmpListValue.push_back(*dnrCoefs);
      }
      else
      {
        BLOCK4;
        PyErr_SetString(PyExc_TypeError,"writeCoefs: invalid array in DonorInterpCoefMap");
        return NULL;
      }      
    }
    donorCoefMap[intKey] = tmpListValue;
  }

  /*--------------------------------------------------------*/
  /* Extraction des types d interpolation                   */
  /*--------------------------------------------------------*/
  map<E_Int,vector<E_Int*> > interpTypesMap;
  vector<PyObject*> pyTypes;
  vector<FldArrayI*> FldTypes;

  if (PyDict_Check (InterpTypesMap) == 0)
  {
    BLOCK4;
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: InterpTypesMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(InterpTypesMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i=0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
      printf("Warning: writeCoefs: invalid integer for variable " SF_D_ ". Skipped...\n", i);
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int*>
    vector<E_Int*> tmpListValue;
    pyListValues = PyDict_GetItem(InterpTypesMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      BLOCK4;
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: InterpTypesMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* intArray = PyList_GetItem(pyListValues, v);
      pyTypes.push_back(intArray);
      FldArrayI* indArrayI;
      E_Int ok = K_NUMPY::getFromNumpyArray(intArray, indArrayI, true);
      if ( ok == 1)
      {
        FldTypes.push_back(indArrayI);
        tmpListValue.push_back(indArrayI->begin());
      }
      else 
      {
        BLOCK5;
        PyErr_SetString(PyExc_TypeError,"writeCoefs: invalid array in InterpTypesMap");
        return NULL;
      }
    }
    interpTypesMap[intKey] = tmpListValue;
  }

  /*--------------------------------------------------------*/
  /* Extraction des listes de champs cellN                  */
  /*--------------------------------------------------------*/
  // DonorCellN : champs cellNatureField
  map<E_Int,FldArrayF> lDonorCellN;
  vector<PyObject*> pyCellN;
  vector<FldArrayF*> FldCellN;

  if (PyDict_Check (DonorCellNMap) == 0)
  {
    BLOCK5;
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: DonorCellNMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(DonorCellNMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable " SF_D_ ". Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    // fill C++ map
    PyObject* pyArrayValue = PyDict_GetItem(DonorCellNMap, pyKey);
    E_Int im, jm, km;
    FldArrayF* cellNF; FldArrayI* cn;
    char* varString; char* eltType;
    E_Int res = K_ARRAY::getFromArray2(pyArrayValue, varString, 
                                       cellNF, im, jm, km, cn, eltType); 
    if ( res != 1 )
    {
      BLOCK6;
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: DonorCellNMap must be a structured array.");
      return NULL;
    }
    E_Int posc = K_ARRAY::isCellNatureField2Present(varString);
    if ( posc != 0 )
    {
      BLOCK6;
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: DonorCellNMap must contain only cellN variable.");
      return NULL;
    }

    pyCellN.push_back(pyArrayValue);
    FldCellN.push_back(cellNF);
    lDonorCellN[intKey]=*cellNF;
  }

  /*--------------------------------------------------------*/
  /* Extraction des listes de dimensions des blocs donneurs */
  /*--------------------------------------------------------*/
  map<E_Int, vector<E_Int> > zoneDimMap;
  if (PyDict_Check (ZoneDimMap) == 0)
  {
    BLOCK6;
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: ZoneDimMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(ZoneDimMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable " SF_D_ ". Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int>
    vector<E_Int> tmpListValue;
    pyListValues = PyDict_GetItem(ZoneDimMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      BLOCK6;
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: ZoneDimMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* pyIntValue = PyList_GetItem(pyListValues, v);
      if (PyInt_Check(pyIntValue) == 0)
      {
        printf("Warning: writeCoefs: invalid int value in  ZoneDimMap\n");
      }
      else
      {
        E_Int intValue = PyInt_AsLong(pyIntValue);
        tmpListValue.push_back(intValue);
      }
    }
   zoneDimMap[intKey] = tmpListValue;
  }

  /*--------------------------------------------------------*/
  /* Ecriture des fichiers d interpolation                  */
  /*--------------------------------------------------------*/
  E_Int typeInterp;
  // Parcours des zones d interpolation
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    // Id du bloc donneur
    E_Int BlockDonorId = noz;

    // Donnees pour la zone noz
    E_Int nbRcvZones = rcvIndexMap[noz].size();

    // 1- Fichier de cellN (necessaire pour une relecture de fichiers des coefficients dans elsA)
    // REMARQUES : - le fichier cellN doit contenir les cellules fictives!
    //             - on peut mettre cellN a 1 sur ces cellules fictives, car la valeur sera ensuite  
    //             ecrasee dans le traitement des raccords
    if (!isEX)
    {
      E_Int ni = zoneDimMap[noz][0]; E_Int nj = zoneDimMap[noz][1]; E_Int nk = zoneDimMap[noz][2];
      E_Int nig = ni+2*nbOfGc; E_Int njg = nj+2*nbOfGc; E_Int nkg = nk+2*nbOfGc;
      E_Int nignjg = nig*njg; E_Int ninj = ni*nj;
      FldArrayI cellNgc;
      E_Float eps = K_CONST::E_GEOM_CUTOFF;
      if (solver == 1) // elsA : prise en compte des cellules fictives pour l indicage
      {
        cellNgc.malloc(nig*njg*nkg); cellNgc.setAllValuesAt(1);
        for (E_Int k=0; k<nk; k++)
          for (E_Int j=0; j<nj; j++)
            for (E_Int i=0; i<ni; i++)
            {
              E_Int ind = k*ninj + j*ni +i;
              E_Int indg = (k+nbOfGc)*nignjg + (j+nbOfGc)*nig +i+nbOfGc;
              if (K_FUNC::fEqualZero(lDonorCellN[noz][ind],eps) == true)
                cellNgc[indg] = -1;
              else if (K_FUNC::fEqualZero(lDonorCellN[noz][ind] - 2.,eps) == true)
                cellNgc[indg] = 0;
            }
      }
      else // Cassiopee : pas de cellules fictives
      {
        cellNgc.malloc(ni*nj*nk); cellNgc.setAllValuesAt(1);
         for (E_Int k=0; k<nk; k++)
          for (E_Int j=0; j<nj; j++)
            for (E_Int i=0; i<ni; i++)
            {
              E_Int ind = k*ninj + j*ni +i;
              if (K_FUNC::fEqualZero(lDonorCellN[noz][ind],eps) == true)
                cellNgc[ind] = -1;
              else if (K_FUNC::fEqualZero(lDonorCellN[noz][ind] - 2.,eps) == true)
                cellNgc[ind] = 0;
            }       
      }
      // Ouverture du fichier
      FILE* ptr_file = NULL;
      char* file = new char[K_ARRAY::VARSTRINGLENGTH];
      strcpy(file,PrefixFile);
      char* strId = new char[K_ARRAY::VARSTRINGLENGTH];
      #if defined E_DOUBLEINT
        #if defined _WIN32
        sprintf(strId,"%04lld",BlockDonorId);
        #else
        sprintf(strId,"%04ld",BlockDonorId);
        #endif
      #else
      sprintf(strId,"%04d",BlockDonorId);
      #endif
      strcat(file,strId);
      strcat(file,"_Blanking");
      ptr_file = fopen(file, "w");
      printf("Open file %s\n",file);fflush(stdout);
      // Ecriture du nombre de points du domaine d interpolation
      E_Int nptsInterp = cellNgc.getSize();
      fprintf(ptr_file, SF_D_ "\n",nptsInterp);
      for (E_Int np = 0; np < nptsInterp; np++)
      {
        fprintf(ptr_file, SF_D_ " ", cellNgc[np]);
      }
      fclose(ptr_file);
      delete [] file; delete [] strId;
    }

    // 1- Fichier d'interpolations
    // Ouverture du fichier
    FILE* ptr_file = NULL;
    char* file = new char[K_ARRAY::VARSTRINGLENGTH];
    strcpy(file, PrefixFile);
    char* strId = new char[K_ARRAY::VARSTRINGLENGTH];
    #if defined E_DOUBLEINT
      #if defined _WIN32
      sprintf(strId,"%04lld",BlockDonorId);
      #else
      sprintf(strId,"%04ld",BlockDonorId);
      #endif
    #else
    sprintf(strId,"%04d",BlockDonorId);
    #endif
    strcat(file, strId);
    if (isEX) strcat(file,"_Int");
    ptr_file = fopen(file, "w");
    printf("Open file %s\n",file);fflush(stdout);
    // Ecriture du nombre de points d interpolations
    E_Int npts = nbInterpCells[noz];
    fprintf(ptr_file, SF_D_ "\n",npts);

    for (E_Int n = 0; n < nbRcvZones; n++)
    {
      E_Int nbOfDatas = donorCoefMap[noz][n].getSize();

      FldArrayF cellN = lDonorCellN[noz];
      
      E_Int blockRcvId = blockRcvIdMap[noz][n];
      E_Int RcvId; 
      // pointeurs sur le champ donorCoefMap de coefficients d'interpolation
      E_Float* donorCoefMap1 = donorCoefMap[noz][n].begin(1); E_Float* donorCoefMap2 = donorCoefMap[noz][n].begin(2); E_Float* donorCoefMap3 = donorCoefMap[noz][n].begin(3); 
      E_Float* donorCoefMap4 = donorCoefMap[noz][n].begin(4); E_Float* donorCoefMap5 = donorCoefMap[noz][n].begin(5); E_Float* donorCoefMap6 = donorCoefMap[noz][n].begin(6); 
      E_Float* donorCoefMap7 = donorCoefMap[noz][n].begin(7); E_Float* donorCoefMap8 = donorCoefMap[noz][n].begin(8);
      
      // Ecriture des donnees sous la forme :
      // indice_pt_interp indice_cell_donneuse id_bloc_rcv coef_interp_1 coef_interp_2 ... coef_interp_8
      for (E_Int np = 0; np < nbOfDatas; np++)
      {
        // indice pour la cellule donneuse exprimee dans le maillage en centres etendus
        E_Int indDonor = (E_Int) donorIndexMap[noz][n][np];
        // indice pour la cellule interpolee
        E_Int indRcv;
        E_Int ni = zoneDimMap[blockRcvId][0]; E_Int nj = zoneDimMap[blockRcvId][1]; E_Int nk = zoneDimMap[blockRcvId][2];
        E_Int ninj = ni*nj;
        E_Int kd  = rcvIndexMap[noz][n][np]/ninj; E_Int val = rcvIndexMap[noz][n][np] - kd*ninj;
        E_Int jd =  val/ni;
        E_Int id = val - jd*ni;
        E_Int nig = ni+2*nbOfGc; E_Int njg = nj+2*nbOfGc; E_Int nkg = nk+2*nbOfGc;
        E_Int nignjg = nig*njg; 
        E_Int kg = kd + nbOfGc; E_Int jg = jd + nbOfGc; E_Int ig = id + nbOfGc;
        typeInterp = (E_Int)interpTypesMap[noz][n][np];
        if (solver == 1) // elsA : prise en compte des cellules fictives pour l indicage
        {
          indRcv = kg*nignjg + jg*nig +ig;
        }
        else // Cassiopee : pas de cellules fictives
        {
          indRcv = kd*ni*nj + jd*ni +id;
        }
        if (solver == 1)  RcvId = blockRcvId; // elsA : les Id des blocs commencent a 0
        else RcvId = blockRcvId+1;// Cassiopee : les Id des blocs commencent a 1
        if (!isEX)
        {
          if (solver == 1) // elsA/Kernel
            fprintf(ptr_file, SF_D4_ " %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",
                    indDonor, indRcv, RcvId, typeInterp, 
                    donorCoefMap1[np], donorCoefMap2[np],donorCoefMap3[np], donorCoefMap4[np], 
                    donorCoefMap5[np], donorCoefMap6[np], donorCoefMap7[np], donorCoefMap8[np]);
          else //Cassiopee/Kernel
            fprintf(ptr_file, SF_D3_ " %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",
                    indDonor, indRcv, RcvId, donorCoefMap1[np], donorCoefMap2[np],donorCoefMap3[np], 
                    donorCoefMap4[np], donorCoefMap5[np], donorCoefMap6[np], donorCoefMap7[np], donorCoefMap8[np]);
        }
        else
        {
          // Calcul de l indice de l interface pour le point EX
          E_Int indIntEX=0;
          if (solver == 1) // elsA : prise en compte des cellules fictives pour l indicage
          {
            RcvId = blockRcvId; // elsA : les Id des blocs commencent a 0
            E_Int kg = kd + nbOfGc; E_Int jg = jd + nbOfGc; E_Int ig = id + nbOfGc;
            E_Int nbIntByDir = nignjg*nkg;
            if      (directionEXMap[noz][n][np] == 0) indIntEX = kg*nignjg+jg*nig+ig+1;                  // interface de frontiere Imax
            else if (directionEXMap[noz][n][np] == 1) indIntEX = kg*nignjg+jg*nig+ig;                    // interface de frontiere Imin
            else if (directionEXMap[noz][n][np] == 2) indIntEX = kg*nignjg+(jg+1)*nig+ig +   nbIntByDir; // interface de frontiere Jmax
            else if (directionEXMap[noz][n][np] == 3) indIntEX = kg*nignjg+jg*nig+ig     +   nbIntByDir; // interface de frontiere Jmin
            else if (directionEXMap[noz][n][np] == 4) indIntEX = (kg+1)*nignjg+jg*nig+ig + 2*nbIntByDir; // interface de frontiere Kmax
            else if (directionEXMap[noz][n][np] == 5) indIntEX = kg*nignjg+jg*nig+ig     + 2*nbIntByDir; // interface de frontiere Kmin
            fprintf(ptr_file, SF_D5_ " %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",
                    indDonor, indIntEX, directionEXMap[noz][n][np], RcvId, typeInterp, 
                    donorCoefMap1[np], donorCoefMap2[np],donorCoefMap3[np], donorCoefMap4[np], 
                    donorCoefMap5[np], donorCoefMap6[np], donorCoefMap7[np], donorCoefMap8[np]);
          }
          else // Cassiopee/Kernel : pas de cellules fictives
          {
            RcvId = blockRcvId+1;// Cassiopee : les Id des blocs commencent a 1
            E_Int  nbIntByDiri = (ni+1)*nj*nk;
            E_Int nbIntByDirj = ni*(nj+1)*nk;
            if      (directionEXMap[noz][n][np] == 0) indIntEX = kd*(ni+1)*nj+jd*(ni+1)+id+1;                               // interface de frontiere Imax
            else if (directionEXMap[noz][n][np] == 1) indIntEX = kd*(ni+1)*nj+jd*(ni+1)+id;                                 // interface de frontiere Imin
            else if (directionEXMap[noz][n][np] == 2) indIntEX = kd*ni*(nj+1)+(jd+1)*ni+id + nbIntByDiri;               // interface de frontiere Jmax
            else if (directionEXMap[noz][n][np] == 3) indIntEX = kd*ni*(nj+1)+jd*ni+id     + nbIntByDiri;               // interface de frontiere Jmin
            else if (directionEXMap[noz][n][np] == 4) indIntEX = (kd+1)*ninj+jd*ni+id + nbIntByDiri + nbIntByDirj; // interface de frontiere Kmax
            else if (directionEXMap[noz][n][np] == 5) indIntEX = kd*ninj+jd*ni+id     + nbIntByDiri + nbIntByDirj; // interface de frontiere Kmin      
            
            fprintf(ptr_file, SF_D4_ " %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",
                    indDonor, indIntEX, directionEXMap[noz][n][np], RcvId, donorCoefMap1[np], donorCoefMap2[np],donorCoefMap3[np], 
                    donorCoefMap4[np], donorCoefMap5[np], donorCoefMap6[np], donorCoefMap7[np], donorCoefMap8[np]);
          }          
        }
      }
    }
    fclose(ptr_file);
    delete [] file; delete [] strId;
  }

  /*-------------------------------*/
  /* Destruction des objets crees  */
  /*-------------------------------*/
//   // destruction de rcvIndexMap
//   map<E_Int,vector<E_Int*> >::iterator itr;
//   for (itr=rcvIndexMap.begin(); itr != rcvIndexMap.end();itr++)
//   {
//     vector<E_Int*> vect = (*itr).second; 
//     for (E_Int w = 0; w < vect.size(); w++) delete [] vect[w]; 
//   }
//   // destruction de donorIndexMap
//   for (itr=donorIndexMap.begin(); itr != donorIndexMap.end();itr++)
//   {
//     vector<E_Int*> vect = (*itr).second; 
//     for (E_Int w = 0; w < vect.size(); w++) delete [] vect[w]; 
//   }

  BLOCK6;
  Py_INCREF(Py_None);
  return Py_None;
}
