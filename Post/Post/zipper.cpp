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

// Surface extractor for Chimera computation

# include "stdio.h"
# include <vector>

# include "post.h"
# include "zipper/zipper.h"

using namespace K_FLD;
using namespace std; 
 
// ============================================================================
/* Surface extractor for Chimera computation */
// ============================================================================
PyObject* K_POST::zipperF(PyObject* self, PyObject* args)
{ 
  PyObject* listFields;
  PyObject* options;
  if (!PyArg_ParseTuple(args, "OO", &listFields, &options))
  {
    return NULL;
  }
  
  // Check every arrays 
  if (PyList_Check(listFields) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "zipper: first argument must be a list of arrays.");
    return NULL;
  }
  
  if (PyList_Check(options) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "zipper: option argument must be a list.");
    return NULL;
  }
  
  // read list of options
  E_Float overlapTol;
  E_Float matchTol;
  readZipperOptions(options, overlapTol, matchTol);

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; 
  vector<E_Int> njt; 
  vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true; 
  E_Boolean skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    listFields, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType,  objst, objut,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = structF.size();
  
  E_Int nfld = 0;
  if ( nzones != 0 ) 
    nfld = structF[0]->getNfld();

  if ( isOk == -1 || nfld < 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "zipper: invalid list of arrays.");
    for (unsigned int nos = 0; nos < objst.size(); nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    return NULL;
  }
  vector<StructBlock*> structBlocks;

  E_Int c = 0;
  char* varString = structVarString[0];
  // coordonnees deja verifiees dans getFromArrays
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); 
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); 
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;

  E_Int posc = K_ARRAY::isCellNatureField1Present(varString);  
  posc++;  
  if ( posc == 0)
  {
    printf("Warning: zipper: no cellnaturefield defined in arrays. Set to 1.\n");
    nfld = nfld+1;

    E_Int structFSize = structF.size();
    for (E_Int v = 0; v < structFSize; v++)
    {
      posc = nfld;
      FldArrayF& field = *structF[v];
      E_Int npts = field.getSize();
      FldArrayF celln0(npts);
      celln0.setAllValuesAt(1.);
      field.setOneField(celln0, posc, 1);
    }
  }
  else if ( posc != nfld )
  {
    PyErr_SetString(PyExc_TypeError,
                    "zipper: celln must be the last variable in all the arrays.");
    for (unsigned int nos = 0; nos < objst.size(); nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    return NULL;
  }
  
  for (E_Int i = 0; i < nzones; i++)
  {
    E_Int nil = nit[i]; E_Int njl = njt[i]; E_Int nkl = nkt[i];
    // check surfaces
    if ( nil == 1 ) 
    {
      nil = njl;
      njl = nkl;
    }
    else if ( njl == 1)
    {
      njl = nkl;
      nkl = 1;
    }
    else if ( nkl == 1)
    {;}
    else
    {
      PyErr_SetString(PyExc_TypeError, 
                      "zipper : array must define a surface." );
      return NULL;
    }
    if ( nil*njl*nkl > 0 )
    {
      FldArrayF& f = *structF[i];
      StructBlock* s = new StructBlock(c, nil, njl, nkl, 
                                       posx, posy, posz,
                                       overlapTol, matchTol, f);
      structBlocks.push_back(s);
      c++;
    }
    else 
      printf("Warning: zipper: one array is empty. Skipped.\n");
    
  }
  for (unsigned int nos = 0; nos < objst.size(); nos++)
    RELEASESHAREDS(objst[nos], structF[nos]);
  if (c == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "zipper: nothing to zip.");
    return NULL;
  }
  //Compute the blanked cells
  computeIBlank(structBlocks);

  // Compute stringing
  vector<CString*> strings;
  computeStringing(structBlocks, strings);

  // Compute pairs of strings
  vector<SegmentPair*> segPairs;
  computeMatchingSegments(strings, segPairs);

  // Zip in between pairs
  vector<FldArrayF*> field;
  vector<FldArrayI*> triConnect;
  // isZipped : retourne True si la paire de segments a ete triangulee
  FldArrayB isZipped(segPairs.size());
  
  zipInbetweenPairs(segPairs,
                    field,
                    triConnect, isZipped);

  // Close remaining pockets
  closePockets(strings, segPairs, field, triConnect);

  E_Int stringsSize = strings.size();
  for (E_Int v = 0; v < stringsSize; v++)
    delete strings[v];
  strings.clear();
  E_Int segPairsSize = segPairs.size();
  for (E_Int v = 0; v < segPairsSize; v++)
    delete segPairs[v];
  segPairs.clear();
  
  // Merge all zones in one
  vector<FldArrayI*> idgT;
  FldArrayI idgG;
  FldArrayF* fieldG = new FldArrayF();
  FldArrayI* unsConnectENG = new FldArrayI();
  FldArrayI FirstPt(structBlocks.size());  
  mergeAllZonesInOne(structBlocks, field, idgT, idgG, *fieldG, FirstPt);

  // Build unstructured connectivity 
  computeConnectivity(structBlocks, triConnect, idgT, idgG, 
                      *unsConnectENG, FirstPt);
    
  for (E_Int v = 0; v < nzones; v++)
    delete structBlocks[v];
  
  structBlocks.clear();
  
  // Deleting unused nodes
  deletingUnusedNodes(*unsConnectENG, *fieldG);
  // Building numpy arrays for unstructured block    

  if ( unsConnectENG->getSize() == 0 || fieldG->getSize() == 0)  
  {
    if ( unsConnectENG->getSize() != 0 ) delete unsConnectENG; 
    if ( fieldG->getSize() != 0 ) delete fieldG;
    PyErr_SetString(PyExc_TypeError,
                    "zipper: no unstructured zone created.");
    return NULL;
  }

  // Build array of unstructured field
  PyObject* tpl2 = K_ARRAY::buildArray(*fieldG, varString,
                                       *unsConnectENG, -1, "TRI");
  delete fieldG; delete unsConnectENG;
  return tpl2;
}

//=============================================================================
// Lecture des options et retourne les valeurs associées: 
// geomTol: tolerance geometrique 
//=============================================================================
void K_POST::readZipperOptions(PyObject* optionList, E_Float& overlapTol,
                               E_Float& matchTol)
{
  //default values 
  matchTol = 1.e-7;
  overlapTol = 1.e-5;
 
  /* parameters reading*/
  E_Int nparam = PyList_Size(optionList);
  if (nparam == 0)
  {
    printf("Warning: zipper: default parameters are set : \n");
    printf("Overlap geometric tolerance is set to 1e-5 .\n");
    printf(" Matching boundaries tolerance is set to 1.e-6.\n");
    return;
  }
  //verification que la liste contient un nombre pair d'elements
  int nparam2 = int(nparam/2);  
  PyObject* tpl1;
  PyObject* tpl2;

  if (nparam*0.5 - nparam2 == 0)
  {
    for (int i = 0; i < nparam2; i++)
    {
      tpl1 = PyList_GetItem(optionList, 2*i);
      tpl2 = PyList_GetItem(optionList, 2*i+1);
    
      char* s = PyString_AsString(tpl1);
      
      // geomTol 
      if ( strcmp(s,"overlapTol") == 0 ) 
      {
        if ( PyFloat_Check(tpl2) == 0)
        {
          printf("Warning: zipper: must be a float. Set to default value 1e-5.\n");
        }
        else overlapTol =  PyFloat_AsDouble(tpl2);
      }
      else if ( strcmp(s, "matchTol") == 0 ) 
      {
        if ( PyFloat_Check(tpl2) == 0)
        {
          printf("Warning: zipper: matchTol must be a float. Set to default value: 1.e-6.\n");
        }
        else matchTol =  PyFloat_AsDouble(tpl2);
      } 
      else //autres 
      {
        printf("Warning: zipper: parameter %s is unknown.\n", s);
      }
    }
  }
}
