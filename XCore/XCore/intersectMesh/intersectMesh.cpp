#include "../common/common.h"
#include "xcore.h"

PyObject *K_XCORE::intersectMesh(PyObject *self, PyObject *args)
{
  PyObject *MASTER, *SLAVE;
  
  if (!PYPARSETUPLE_(args, OO_, &MASTER, &SLAVE)) {
    RAISE("Bad input.");
    return NULL;
  }


  return Py_None;
}