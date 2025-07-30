#include "template.h"

// Recupere une numpy de python au C 
PyObject* K_TEMPLATE::numpyExample(PyObject* self, PyObject* args)
{
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;

  /* More info in KCore/Numpy/Numpy.h */
  E_Float* f; E_Int size; E_Int nfld;
  E_Int ret = K_NUMPY::getFromNumpyArray(o, f, size, nfld);

  for (E_Int n = 0; n < nfld; n++)
    for (E_Int i = 0; i < size; i++) f[i+size*n] = 12.;

  Py_INCREF(Py_None);
  return Py_None;
}

// Recupere un numpy de python et le passe au fortran
extern "C" 
{
  void myfunction_(const E_Int& n, E_Float* a);
}
PyObject* K_TEMPLATE::numpyExample2(PyObject* self, PyObject* args)
{
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;

  /* More info in KCore/Numpy/Numpy.h */
  E_Float* f; E_Int size; E_Int nfld;
  E_Int ret = K_NUMPY::getFromNumpyArray(o, f, size, nfld);

  myfunction_(size*nfld, f);

  Py_INCREF(Py_None);
  return Py_None;
}
