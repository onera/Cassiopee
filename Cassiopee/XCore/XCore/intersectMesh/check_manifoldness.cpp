#include "mesh.h"
#include "common/Karray.h"

PyObject *K_XCORE::check_manifoldness(PyObject *self, PyObject *args)
{
    PyObject *MESH;
    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray array;
    if (Karray_parse_ngon(MESH, array) != 0) {
        RAISE("Input not NGon.");
        return NULL;
    }

    IMesh M(array);

    bool is_manifold = true;

    M.make_skin();
    Smesh Mf(M, M.skin, false);

    Karray_free_ngon(array);

    return PyBool_FromLong(is_manifold);
}