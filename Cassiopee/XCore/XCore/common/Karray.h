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
#pragma once

#include "xcore.h"
#include "common/common.h"

struct Karray {
    // Reference to the python array
    PyObject *pyobject;

    K_FLD::FldArrayF *f;
    E_Int ni, nj, nk;

    K_FLD::FldArrayI *cn;

    E_Float *x, *y, *z;
    E_Int npts;

    inline E_Int npoints() const { return npts; }
    inline E_Int nfaces() const { return cn->getNFaces(); }
    inline E_Int ncells() const { return cn->getNElts(); }
    
    inline E_Int *indpg() const { return cn->getIndPG(); }
    inline E_Int *indph() const { return cn->getIndPH(); }

    inline E_Int *ngon() const { return cn->getNGon(); }
    inline E_Int *nface() const { return cn->getNFace(); }

    inline E_Int *get_face(E_Int face, E_Int &np) const
    { return cn->getFace(face, np, ngon(), indpg()); }

    inline E_Int *get_cell(E_Int cell, E_Int &nf) const
    { return cn->getElt(cell, nf, nface(), indph()); }

    inline E_Int orient_boundary()
    { return K_CONNECT::orient_boundary_ngon(x, y, z, *cn); }

    inline E_Int build_parent_elements(E_Int *owner, E_Int *neigh)
    { return K_CONNECT::build_parent_elements_ngon(*cn, owner, neigh); }
};

void Karray_free_ngon(Karray &karray);

E_Int Karray_parse_ngon(PyObject *pyobject, Karray &karray);

void Karray_free_structured(Karray &karray);

E_Int Karray_parse_structured(PyObject *pyobject, Karray &karray);
