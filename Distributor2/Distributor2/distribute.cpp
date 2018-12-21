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

# include "distributor2.h"
using namespace std;

//============================================================================
/* Distribute blocks among N processors */
//============================================================================
PyObject* K_DISTRIBUTOR2::distribute(PyObject* self, PyObject* args)
{
  PyObject* nbPts; PyObject* setBlocks;
  PyObject* perfProcs; PyObject* weight;
  PyObject* com; E_Int NProc;
  char* algorithm;
  if (!PYPARSETUPLEI(args,"OOOOOls","OOOOOis", 
                     &nbPts, &setBlocks, &perfProcs, &weight, 
                     &com, &NProc, &algorithm)) return NULL;

  // Check algorithm
  E_Int algo, param;
  if (strcmp(algorithm, "gradient") == 0) {algo=0; param=0;}
  else if (strcmp(algorithm, "gradient0") == 0) {algo=0; param=0;}
  else if (strcmp(algorithm, "gradient1") == 0) {algo=0; param=1;}
  else if (strcmp(algorithm, "genetic") == 0) {algo=1; param=0;}
  else if (strcmp(algorithm, "fast") == 0) {algo=1; param=1;}
  else if (NProc > 1 && strcmp(algorithm, "graph") == 0) {algo=2; param=0;}
  else if (NProc == 1 && strcmp(algorithm, "graph") == 0) {algo=1; param=1;}
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "distribute: algorithm is unknown.");
    return NULL;
  }

  // Check args
  if (PyList_Check(nbPts) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distribute: arrays must be a list.");
    return NULL;
  }
  // Nombre de blocs
  E_Int nb = PyList_Size(nbPts);

  // Analyse setBlocks
  if (PyList_Check(setBlocks) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distribute: set must be a list.");
    return NULL;
  }
  if (PyList_Size(setBlocks) != nb)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distribute: arrays and set must have the same size.");
    return NULL;
  }
 
  // Analyse weight
  if (PyList_Check(weight) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distribute: weight must be a list.");
    return NULL;
  }
  if (PyList_Size(weight) != nb)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distribute: arrays and weight must have the same size.");
    return NULL;
  }

  // Analyse de com
  IMPORTNUMPY;
  if (PyArray_Check(com) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distribute: com must a numpy array.");
    return NULL;
  }
  
  PyArrayObject* coma = (PyArrayObject*)com;
  int* volCom = (int*)PyArray_DATA(coma);
  Py_INCREF(coma);

  // Analyse perfProcs
  if (PyList_Check(perfProcs) == 0)
  {
    Py_DECREF(coma);
    PyErr_SetString(PyExc_TypeError,
                    "distribute: perfo must be a list.");
    return NULL;
  }
  if (PyList_Size(perfProcs) != NProc)
  {
    Py_DECREF(coma);
    PyErr_SetString(PyExc_TypeError,
                    "distribute: perfo must have NProc elements.");
    return NULL;
  }

  vector<E_Float> solver; vector<E_Float> latence; vector<E_Float> comSpeed;
  E_Float valf;
  for (E_Int i = 0; i < NProc; i++)
  {
    PyObject* o = PyList_GetItem(perfProcs, i);
    if (PyTuple_Check(o) == 0 || PyTuple_Size(o) != 3)
    {
      Py_DECREF(coma);
      PyErr_SetString(PyExc_TypeError,
                      "distribute: perfo must be tuples (solver,latency,comSpeed).");
      return NULL;
    }
    
    valf = PyFloat_AsDouble(PyTuple_GetItem(o, 0));
    solver.push_back(valf);
    valf = PyFloat_AsDouble(PyTuple_GetItem(o, 1));
    latence.push_back(valf);
    valf = PyFloat_AsDouble(PyTuple_GetItem(o, 2));
    comSpeed.push_back(valf);
  }
  
  // Recupere le vecteur des poids du solveur pour les blocs
  vector<E_Float> lweight;
  E_Float val;
  for (E_Int i = 0; i < nb; i++)
  {
    val = PyFloat_AsDouble(PyList_GetItem(weight, i));
    lweight.push_back(val);
  }

  // Recupere le vecteur du nombre de pts de chaque blocs dans une liste c++
  vector<E_Float> lnbPts;
  E_Int ival;
  for (E_Int i = 0; i < nb; i++)
  {
    ival = PyLong_AsLong(PyList_GetItem(nbPts, i));
    val = lweight[i] * double(ival);
    lnbPts.push_back(val);
  }

  // Recupere le vecteur des blocs imposes
  vector<E_Int> lsetBlocks;
  for (E_Int i = 0; i < nb; i++)
  {
    ival = PyLong_AsLong(PyList_GetItem(setBlocks, i));
    lsetBlocks.push_back(ival);
  }

  // Algo genetique ou gradient
  vector<E_Int> out(nb);
  E_Float meanPtsPerProc;
  E_Float varMin, varMax, varRMS, comRatio, bestAdapt;
  E_Int nptsCom; 
  if (algo == 1) // genetic
    K_DISTRIBUTOR2::genetic(lnbPts, lsetBlocks, NProc, volCom, 
                            solver, latence, comSpeed, param, out,
                            meanPtsPerProc, varMin,
                            varMax, varRMS, nptsCom, comRatio,
                            bestAdapt);
  else if (algo == 2) // graph
  {
    K_DISTRIBUTOR2::graph(lnbPts, lsetBlocks, NProc, volCom, 
                          solver, latence, comSpeed, param, out,
                          meanPtsPerProc, varMin,
                          varMax, varRMS, nptsCom, comRatio,
                          bestAdapt);
  }
  else // algo = 0
    K_DISTRIBUTOR2::gradient(lnbPts, lsetBlocks, NProc, volCom, 
                             solver, latence, comSpeed, param, out,
                             meanPtsPerProc, varMin,
                             varMax, varRMS, nptsCom, comRatio, 
                             bestAdapt);
    
  // Formation de la sortie 
  PyObject* tpl = PyList_New(0);
  E_Int size = out.size();
  PyObject* o;
  for (E_Int i = 0; i < size; i++)
  {
    o = Py_BuildValue("l", out[i]);
    PyList_Append(tpl, o); Py_DECREF(o);
  }
  Py_DECREF(coma);

  // Formation du dictionnaire de stats
  PyObject* stats = PyDict_New();
  PyDict_SetItemString(stats, "distrib", tpl); Py_DECREF(tpl);
  o = Py_BuildValue("d", meanPtsPerProc);
  PyDict_SetItemString(stats, "meanPtsPerProc", o); Py_DECREF(o);
  o = Py_BuildValue("d", varMin);
  PyDict_SetItemString(stats, "varMin", o); Py_DECREF(o);
  o = Py_BuildValue("d", varMax);
  PyDict_SetItemString(stats, "varMax", o); Py_DECREF(o);
  o = Py_BuildValue("d", varRMS);
  PyDict_SetItemString(stats, "varRMS", o); Py_DECREF(o);
  o = Py_BuildValue("l", nptsCom);
  PyDict_SetItemString(stats, "nptsCom", o); Py_DECREF(o);
  o = Py_BuildValue("d", comRatio);
  PyDict_SetItemString(stats, "comRatio", o); Py_DECREF(o);
  o = Py_BuildValue("d", bestAdapt);
  PyDict_SetItemString(stats, "adaptation", o); Py_DECREF(o);
  return stats;
}
