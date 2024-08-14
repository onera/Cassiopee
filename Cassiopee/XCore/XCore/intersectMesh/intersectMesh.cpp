/*    
    Copyright 2013-2024 Onera.

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
#include "xcore.h"
#include "karray.h"
#include "common/common.h"
#include "mesh.h"
#include "smesh.h"
#include "dcel.h"
#include "vertex.h"
#include "face.h"
#include "hedge.h"
#include "io.h"
#include "cycle.h"
#include "triangle.h"
#include "primitives.h"

PyObject *K_XCORE::intersectMesh(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVE, *MPATCH, *SPATCH;
  
    if (!PYPARSETUPLE_(args, OOOO_, &MASTER, &SLAVE, &MPATCH, &SPATCH)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;
    Karray sarray;

    Int ret;

    ret = Karray_parse_ngon(MASTER, marray);

    if (ret != 0) return NULL;

    ret = Karray_parse_ngon(SLAVE, sarray);

    if (ret != 0) {
        Karray_free_ngon(marray);
        return NULL;
    }

    // Init and orient master/slave meshes
    IMesh M(*marray.cn, marray.X, marray.Y, marray.Z, marray.npts);
    IMesh S(*sarray.cn, sarray.X, sarray.Y, sarray.Z, sarray.npts);

    // Check master intersection patch (zero-based)
    Int *mpatch = NULL;
    Int mpatch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(MPATCH, mpatch, mpatch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        RAISE("Bad master patch.");
        return NULL;
    }

    printf("Master patch: " SF_D_ " faces\n", mpatch_size);

    // Check slave intersection patch (zero-based)
    Int *spatch = NULL;
    Int spatch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(SPATCH, spatch, spatch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        RAISE("Bad slave patch.");
        return NULL;
    }

    printf("Slave patch: " SF_D_ " faces\n", spatch_size);

    for (Int i = 0; i < mpatch_size; i++) M.patch.insert(mpatch[i]);
    for (Int i = 0; i < spatch_size; i++) S.patch.insert(spatch[i]);

    M.orient_skin(OUT);
    S.orient_skin(IN);

    // Extract surface meshes
    Smesh Mf(M);
    Smesh Sf(S);

    Mf.write_ngon("Mf");
    Sf.write_ngon("Sf");

    Mf.make_point_edges();
    Sf.make_point_edges();

    Mf.make_point_faces_all();
    Sf.make_point_faces_all();

    Mf.make_pnormals();
    Sf.make_pnormals();

    Dcel D(Mf, Sf);

    D.locate_spoints(Mf, Sf);

    D.init_Cp(Mf, Sf);

    D.find_intersections_3D(Mf, Sf);

    return Py_None;
}
