#pragma once

#include "xcore.h"
#include "../common/common.h"

struct Karray {
    // Reference to the python array
    PyObject *pyobject;

    K_FLD::FldArrayF *f;
    Int ni, nj, nk;

    K_FLD::FldArrayI *cn;

    E_Float *x, *y, *z;
    Int npts;

    inline Int npoints() const { return npts; }
    inline Int nfaces() const { return cn->getNFaces(); }
    inline Int ncells() const { return cn->getNElts(); }
    
    inline Int *indpg() const { return cn->getIndPG(); }
    inline Int *indph() const { return cn->getIndPH(); }

    inline Int *ngon() const { return cn->getNGon(); }
    inline Int *nface() const { return cn->getNFace(); }

    inline Float *X() const { return x; }
    inline Float *Y() const { return y; }
    inline Float *Z() const { return z; }

    inline Int *get_face(Int face, Int &np) const
    { return cn->getFace(face, np, ngon(), indpg()); }

    inline Int *get_cell(Int cell, Int &nf) const
    { return cn->getElt(cell, nf, nface(), indph()); }

    inline Int orient_boundary()
    { return K_CONNECT::orient_boundary_ngon(X(), Y(), Z(), *cn); }

    inline Int build_parent_elements(Int *owner, Int *neigh)
    { return K_CONNECT::build_parent_elements_ngon(*cn, owner, neigh); }
};

void Karray_free_ngon(Karray &karray);

Int Karray_parse_ngon(PyObject *pyobject, Karray &karray);

void Karray_free_structured(Karray &karray);

Int Karray_parse_structured(PyObject *pyobject, Karray &karray);