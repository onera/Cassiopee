/*    
    Copyright 2013-2026 ONERA.

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
#include "icapsule.h"
#include "dcel.h"

PyObject *K_XCORE::icapsule_intersect2(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;

    if (!PYPARSETUPLE_(args, O_, &ICAPSULE)) {
        RAISE("Bad input");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICAPSULE")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICAPSULE");

    auto &M = icap->M;
    auto &Ss = icap->Ss;

    M.make_skin();
    M.orient_skin(OUT);
    M.Mf = Smesh(M, M.skin, false);
    auto &Mf = M.Mf;

    for (size_t i = 0; i < Ss.size(); i++) {
        printf("Intersecting slave %zu ...\n", i);

        Mf.make_fcenters();
        Mf.make_fnormals();
        Mf.make_pnormals();
        Mf.make_point_faces();
        Mf.make_BVH();

        auto &S = Ss[i];

        Smesh Sf = Smesh::Smesh_from_tagged_faces(S, true);
        Sf.make_fcenters();
        Sf.make_fnormals();
        Sf.make_pnormals();
        
        auto plocs = Mf.locate2(Sf);

        Dcel D = Dcel::intersect(Mf, Sf, plocs);

        D.reconstruct(Mf, Dcel::RED);
        D.reconstruct(Sf, Dcel::BLACK);

        Sf.reconstruct(S);

        // Tag Sf faces
        Sf.tag_faces(S);

        printf("Done.\n");
        fflush(stdout);
    }
    
    Mf.reconstruct(M);

    return Py_None;
}