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

//=============================================================================
/* gatherMatchingNGon */
//=============================================================================
PyObject* K_CONNECTOR::gatherMatchingNGon(PyObject* self, PyObject* args)
{
  // INPUT: arrays of tag1 (number of the opposite zone)
  //                  tag2 (number of the opposite element)
  //        list of indices for all exterior faces
  PyObject* AllTags;
  PyObject* OriginalExteriorFaceIndices;
  if (!PYPARSETUPLE_(args, OO_, &AllTags, &OriginalExteriorFaceIndices))
  {
    return NULL;
  }

  /* Extract info on tag1 and tag2 */
  vector<E_Int> rest;
  vector<char*> structVarStringt; vector<char*> unstrVarStringt;
  vector<FldArrayF*> structFt; vector<FldArrayF*> unstrFt;
  vector<E_Int> nitt; vector<E_Int> njtt; vector<E_Int> nktt;
  vector<FldArrayI*> cntt;
  vector<char*> eltTypett;
  vector<PyObject*> objstt, objutt;
  bool skipNoCoord = false;
  bool skipStructured = true;
  bool skipUnstructured = false;
  bool skipDiffVars = false;

  K_ARRAY::getFromArrays(
    AllTags, rest, structVarStringt, unstrVarStringt,
    structFt, unstrFt, nitt, njtt, nktt, cntt, eltTypett, objstt, objutt, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = unstrFt.size();
  vector<E_Int> post1; vector<E_Int> post2;
  E_Int postag1, postag2;
  for (E_Int i = 0; i < nzones; i++)
  {   
    postag1 = K_ARRAY::isNamePresent("tag1", unstrVarStringt[i]); postag1++;
    postag2 = K_ARRAY::isNamePresent("tag2", unstrVarStringt[i]); postag2++;
    if (postag1 < 1 || postag2 < 1) 
    {
      PyErr_SetString(PyExc_TypeError, "gatherMatchingNGon: tag1 and tag2 must be defined in listOfAllWins.");
      for (E_Int is = 0; is < nzones; is++)
      {
        RELEASESHAREDU(objutt[is], unstrFt[is], cntt[is]);
      }
      return NULL;
    } 
    post1.push_back(postag1); post2.push_back(postag2);
  }

  // Original face numbering in volume mesh : extract numpys
  vector<FldArrayI*> vectOfOrigIndicesFaces;
  vector<PyObject*> vectOfIndicesObjects;
  int sizeOfOrigIndices = 0;
  if (PyList_Check(OriginalExteriorFaceIndices) != 0)
  {   
    sizeOfOrigIndices = PyList_Size(OriginalExteriorFaceIndices);
    for (int i = 0; i < sizeOfOrigIndices; i++)
    {
      PyObject* tpl = PyList_GetItem(OriginalExteriorFaceIndices, i);
      FldArrayI* indicesFacesL;
      K_NUMPY::getFromNumpyArray(tpl, indicesFacesL);
      vectOfOrigIndicesFaces.push_back(indicesFacesL);
      vectOfIndicesObjects.push_back(tpl);
    }
  }
  
  if (nzones != sizeOfOrigIndices)
  {
    PyErr_SetString(PyExc_TypeError,
                    "gatherMatchingNGon: 1st and 2nd args must be of same length.");
    for (E_Int is = 0; is < nzones; is++)
      RELEASESHAREDU(objutt[is], unstrFt[is], cntt[is]);
    return NULL;
  }
  // Parcours des tags
  E_Int nfacesExt, nozopp, nofopp, nofaceorig, nofaceorigopp;
 
  PyObject* allPyOrigFacesR = PyList_New(0);  
  PyObject* allPyOrigFacesD = PyList_New(0);  

  vector<E_Int> noZR; vector<E_Int> noZD;

  for (E_Int noz = 0; noz < nzones; noz++)
  {
    nfacesExt = unstrFt[noz]->getSize(); 
    E_Int* origFaces = vectOfOrigIndicesFaces[noz]->begin();
    E_Float* oppZones = unstrFt[noz]->begin(post1[noz]);
    E_Float* oppFaces = unstrFt[noz]->begin(post2[noz]);
    vector< vector<E_Int> > origFacesR(nzones);
    vector< vector<E_Int> > origFacesD(nzones);
    
    for (E_Int nof = 0; nof < nfacesExt; nof++)
    {
      nozopp = E_Int(oppZones[nof]);
      nofopp = E_Int(oppFaces[nof]);
      if (nozopp>-1 && nofopp>-1)
      {
        nofaceorig = origFaces[nof];
        //vector<E_Int> origFacesRL = origFacesR[nozopp];
        //vector<E_Int> origFacesDL = origFacesD[nozopp];
        E_Int* origFacesOpp = vectOfOrigIndicesFaces[nozopp]->begin();
        nofaceorigopp = origFacesOpp[nofopp];
        origFacesR[nozopp].push_back(nofaceorig);
        origFacesD[nozopp].push_back(nofaceorigopp);
        //oppZones[nof] = -1.; 
        oppFaces[nof] = -1.;
        //E_Float* currZones = unstrFt[nozopp]->begin(post1[nozopp]);
        E_Float* currFaces = unstrFt[nozopp]->begin(post2[nozopp]);
        //currZones[nofopp] = -1.; 
        currFaces[nofopp] = -1.;
      }
    }
    // on cree les tableaux 
    for (E_Int nozopp = 0; nozopp < nzones; nozopp++)
    {
      vector<E_Int> origFacesRL = origFacesR[nozopp];
      vector<E_Int> origFacesDL = origFacesD[nozopp];
      nfacesExt = origFacesRL.size();
      if ( nfacesExt > 0 )
      {
        PyObject* PyOrigFacesR = K_NUMPY::buildNumpyArray(nfacesExt, 1, 1, 0); 
        PyObject* PyOrigFacesD = K_NUMPY::buildNumpyArray(nfacesExt, 1, 1, 0); 
        E_Int* oneOrigFacesR = K_NUMPY::getNumpyPtrI(PyOrigFacesR);
        E_Int* oneOrigFacesD = K_NUMPY::getNumpyPtrI(PyOrigFacesD);
        for (E_Int nof = 0; nof < nfacesExt; nof++)
        {
          oneOrigFacesR[nof] = origFacesRL[nof];
          oneOrigFacesD[nof] = origFacesDL[nof];
        }
        PyList_Append(allPyOrigFacesR, PyOrigFacesR); Py_DECREF(PyOrigFacesR);
        PyList_Append(allPyOrigFacesD, PyOrigFacesD); Py_DECREF(PyOrigFacesD);
        origFacesRL.clear(); origFacesDL.clear();
        noZR.push_back(noz); noZD.push_back(nozopp);
      }      
    }
  }
  E_Int sizeTot = noZR.size();
  PyObject* allPyZonesD = K_NUMPY::buildNumpyArray(sizeTot,1, 1, 0);
  PyObject* allPyZonesR = K_NUMPY::buildNumpyArray(sizeTot,1, 1, 0);
  E_Int* allZonesR = K_NUMPY::getNumpyPtrI(allPyZonesR);
  E_Int* allZonesD = K_NUMPY::getNumpyPtrI(allPyZonesD);
  for (E_Int no = 0; no < sizeTot; no++)
  { 
    allZonesR[no] = noZR[no];
    allZonesD[no] = noZD[no];
  }
  
  // Cleaning
  for (E_Int is = 0; is < nzones; is++)
  {
    RELEASESHAREDU(objutt[is], unstrFt[is], cntt[is]);
  }
  for (size_t i = 0; i < vectOfOrigIndicesFaces.size(); i++)
  {
    RELEASESHAREDN(vectOfIndicesObjects[i], vectOfOrigIndicesFaces[i]);
  }

  // Sortie
  PyObject* tpl = Py_BuildValue("[OOOO]", allPyZonesR, allPyZonesD, allPyOrigFacesR, allPyOrigFacesD);
  Py_DECREF(allPyZonesD); Py_DECREF(allPyZonesR); Py_DECREF(allPyOrigFacesR);Py_DECREF(allPyOrigFacesD);  
  return tpl;
}
 
