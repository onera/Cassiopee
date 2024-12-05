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
#include "common/Karray.h"

PyObject *K_XCORE::write_bim_s(PyObject *self, PyObject *args)
{
    PyObject *ARRAY, *FNAME;
    if (!PYPARSETUPLE_(args, OO_, &ARRAY, &FNAME)) {
        RAISE("Bad input.");
        return NULL;
    }

    const char *fname = PyUnicode_AsUTF8(FNAME);
    if (!fname) {
        RAISE("Bad output file name");
        return NULL;
    }

    Karray k;
    E_Int ret = Karray_parse_structured(ARRAY, k);
    if (ret != 0) {
        RAISE("Bad structured mesh");
        return NULL;
    }

    FILE *fh = fopen(fname, "wb");
    if (!fh) {
        RAISE("Could not open output file.");
        Karray_free_ngon(k);
        return NULL;
    }

    fwrite(&k.ni, sizeof(E_Int), 1, fh);
    fwrite(&k.nj, sizeof(E_Int), 1, fh);
    fwrite(&k.nk, sizeof(E_Int), 1, fh);
    fwrite(&k.npts, sizeof(E_Int), 1, fh);
    fwrite(k.x, sizeof(E_Float), k.npts, fh);
    fwrite(k.y, sizeof(E_Float), k.npts, fh);
    fwrite(k.z, sizeof(E_Float), k.npts, fh);

    fclose(fh);

    Karray_free_structured(k);

    return Py_None;
}

PyObject *K_XCORE::write_bim(PyObject *self, PyObject *args)
{
    PyObject *ARRAY, *FNAME;
    if (!PYPARSETUPLE_(args, OO_, &ARRAY, &FNAME)) {
        RAISE("Bad input.");
        return NULL;
    }

    const char *fname = PyUnicode_AsUTF8(FNAME);
    if (!fname) {
        RAISE("Bad output file name.");
        return NULL;
    }

    Karray k;
    E_Int ret = Karray_parse_ngon(ARRAY, k);
    if (ret != 0) {
        RAISE("Bad NGon.");
        return NULL;
    }

    FILE *fh = fopen(fname, "wb");
    if (!fh) {
        RAISE("Could not open output file.");
        Karray_free_ngon(k);
        return NULL;
    }

    E_Int np = k.npoints();
    E_Int nf = k.nfaces();
    E_Int nc = k.ncells();
    E_Float *X = k.x, *Y = k.y, *Z = k.z;

    E_Int *INDPG = new E_Int[nf+1];
    INDPG[0] = 0;
    for (E_Int i = 0; i < nf; i++) {
        E_Int NP = -1;
        k.get_face(i, NP);
        INDPG[i+1] = INDPG[i] + NP;
    }

    E_Int *INDPH = new E_Int[nc+1];
    INDPH[0] = 0;
    for (E_Int i = 0; i < nc; i++) {
        E_Int NF = -1;
        k.get_cell(i, NF);
        INDPH[i+1] = INDPH[i] + NF;
    }

    E_Int *NGON = new E_Int[INDPG[nf]];
    E_Int *NFACE = new E_Int[INDPH[nc]];

    E_Int *ptr = NGON;

    for (E_Int i = 0; i < nf; i++) {
        E_Int NP = -1;
        E_Int *pn = k.get_face(i, NP);
        for (E_Int j = 0; j < NP; j++) *ptr++ = pn[j];
    }

    ptr = NFACE;

    for (E_Int i = 0; i < nc; i++) {
        E_Int NF = -1;
        E_Int *pf = k.get_cell(i, NF);
        for (E_Int j = 0; j < NF; j++) *ptr++ = pf[j];
    }

    E_Int dummy[1];

    // Points
    printf("np: %d\n", np);
    dummy[0] = np;
    fwrite(dummy, sizeof(E_Int), 1, fh);
    fwrite(X, sizeof(E_Float), np, fh);
    fwrite(Y, sizeof(E_Float), np, fh);
    fwrite(Z, sizeof(E_Float), np, fh);

    // IndPG
    printf("nf: %d\n", nf);
    dummy[0] = nf+1;
    fwrite(dummy, sizeof(E_Int), 1, fh);
    fwrite(INDPG, sizeof(E_Int), nf+1, fh);

    // Ngon
    printf("size ngon: %d\n", INDPG[nf]);
    fwrite(NGON, sizeof(E_Int), INDPG[nf], fh);

    // IndPH
    printf("nc: %d\n", nc);
    dummy[0] = nc+1;
    fwrite(dummy, sizeof(E_Int), 1, fh);
    fwrite(INDPH, sizeof(E_Int), nc+1, fh);

    // Ngon
    printf("size nface: %d\n", INDPH[nc]);
    fwrite(NFACE, sizeof(E_Int), INDPH[nc], fh);

    fclose(fh);

    Karray_free_ngon(k);

    delete [] INDPH;
    delete [] INDPG;
    delete [] NFACE;
    delete [] NGON;

    return Py_None;
}

PyObject *K_XCORE::write_im(PyObject *self, PyObject *args)
{
    PyObject *ARRAY, *FNAME;
    if (!PYPARSETUPLE_(args, OO_, &ARRAY, &FNAME)) {
        RAISE("Bad input.");
        return NULL;
    }

    const char *fname = PyUnicode_AsUTF8(FNAME);
    if (!fname) {
        RAISE("Bad output file name.");
        return NULL;
    }

    Karray k;
    E_Int ret = Karray_parse_ngon(ARRAY, k);
    if (ret != 0) {
        RAISE("Bad NGon.");
        return NULL;
    }

    FILE *fh = fopen(fname, "w");
    if (!fh) {
        RAISE("Could not open output file.");
        Karray_free_ngon(k);
        return NULL;
    }

    E_Int np = k.npoints();
    E_Int nf = k.nfaces();
    E_Int nc = k.ncells();
    E_Float *X = k.x, *Y = k.y, *Z = k.z;
    E_Int *indpg = k.indpg(), *indph = k.indph();
    E_Int *ngon = k.ngon(), *nface = k.nface();

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%d\n", np);
    for (E_Int i = 0; i < np; i++) {
        fprintf(fh, "%f %f %f ", X[i], Y[i], Z[i]);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPG\n");
    fprintf(fh, "%d\n", nf+1);
    for (E_Int i = 0; i < nf+1; i++) fprintf(fh, "%d ", indpg[i]);
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, "%d\n", indpg[nf]);
    for (E_Int i = 0; i < indpg[nf]; i++) fprintf(fh, "%d ", ngon[i]);
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, "%d\n", nc+1);
    for (E_Int i = 0; i < nc+1; i++) fprintf(fh, "%d ", indph[i]);
    fprintf(fh, "\n");

    fprintf(fh, "NFACE\n");
    fprintf(fh, "%d\n", indph[nc]);
    for (E_Int i = 0; i < indph[nc]; i++) fprintf(fh, "%d ", nface[i]);
    fprintf(fh, "\n");

    fclose(fh);

    Karray_free_ngon(k);

    return Py_None;
}
