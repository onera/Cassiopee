/* Exemple de passage d'entiers et de float de python au C */

#include "template.h"

PyObject* K_TEMPLATE::intFloatExample(PyObject* self, PyObject* args)
{
  E_Float f; E_Int n;
#if defined E_DOUBLEREAL && defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "dl", &f, &n)) return NULL;
#elif defined E_DOUBLEREAL && !defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "di", &f, &n)) return NULL;
#elif !defined E_DOUBLEREAL && defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "fl", &f, &n)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "fi", &f, &n)) return NULL;
#endif 
  printf("Read %f and %d.\n", f, n);
  
  // In python, returning NULL means an error
  // Nothing must be returned by Py_None;
  Py_INCREF(Py_None);
  return Py_None;
}
