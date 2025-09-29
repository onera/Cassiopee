/* Exemple de passage d'entiers et de float de python au C */

#include "template.h"

PyObject* K_TEMPLATE::intFloatExample(PyObject* self, PyObject* args)
{
  E_Float f; E_Int n;
  if (!PYPARSETUPLE_(args, R_ I_, &f, &n)) return NULL;
  printf("Read " SF_F_ " and " SF_D_ ".\n", f, n);

  // In python, returning NULL means an error
  // Nothing must be returned by Py_None;
  Py_INCREF(Py_None);
  return Py_None;
}
