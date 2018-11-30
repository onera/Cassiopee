#include "kcore.h"
PyObject* K_KCORE::activation(PyObject* self, PyObject* args)
{
  int date = activation();
  return Py_BuildValue("l", date);
}
int K_KCORE::activation(const char* name)
{
  return 0;
}
