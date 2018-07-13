// Distributor

# include "Python.h"
# include "distributor.h"
# include "Spl/Api/SplApi.h"

static PyObject*
PyregisterAsSplitter( PyObject* dp, PyObject* args )
{
  SplApi::setImplementation( SplDistributor() );
  Py_INCREF(Py_None);
  return Py_None;
}

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pydistributor [] =
{
  {"register", PyregisterAsSplitter, METH_VARARGS},
  {NULL, NULL}
};
// ============================================================================
/* Init of module */
// ============================================================================
extern"C"
{
  PyMODINIT_FUNC initdistributor();
  PyMODINIT_FUNC initdistributor()
  {
    Py_InitModule("distributor", 
                  Pydistributor);
  }
}

