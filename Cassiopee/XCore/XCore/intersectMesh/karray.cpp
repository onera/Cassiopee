#include "karray.h"
#include "../common/common.h"

void Karray_free_ngon(Karray &karray)
{
    RELEASESHAREDU(karray.pyobject, karray.f, karray.cn);
}

E_Int Karray_parse_ngon(PyObject *pyobject, Karray &karray)
{
    E_Int ret;

    char *varString;
    char *eltType;
    
    ret = K_ARRAY::getFromArray3(pyobject, varString, karray.f, karray.ni,
        karray.nj, karray.nk, karray.cn, eltType);
    
    if (ret <= 0) {
        RAISE("Bad input array.");
        return 1;
    }

    if (ret == 1) {
        RAISE("Mesh should be an NGon.");
        Karray_free_structured(karray);
        return 1;
    }
    
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1) {
        Karray_free_ngon(karray);
        RAISE("Coordinates not found.");
        return 1;
    }

    posx++; posy++; posz++;

    karray.X = karray.f->begin(posx);
    karray.Y = karray.f->begin(posy);
    karray.Z = karray.f->begin(posz);
    karray.npts = karray.f->getSize();

    karray.pyobject = pyobject;

    return 0;
}

void Karray_free_structured(Karray &karray)
{
    RELEASESHAREDS(karray.pyobject, karray.f);
}

E_Int Karray_parse_structured(PyObject *pyobject, Karray &karray)
{
    E_Int ret;

    char *varString;
    char *eltType;

    ret = K_ARRAY::getFromArray3(pyobject, varString, karray.f, karray.ni,
        karray.nj, karray.nk, karray.cn, eltType);

    if (ret <= 0) {
        RAISE("Bad input array");
        return 1;
    }

    if (ret != 1) {
        RAISE("Mesh should be structured.");
        RELEASESHAREDB(ret, pyobject, karray.f, karray.cn);
        return 1;
    }

    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1) {
        RAISE("Coordinates not found");
        Karray_free_structured(karray);
        return 1;
    }

    posx++; posy++; posz++;

    karray.X = karray.f->begin(posx);
    karray.Y = karray.f->begin(posy);
    karray.Z = karray.f->begin(posz);
    karray.npts = karray.f->getSize();

    karray.pyobject = pyobject;

    return 0;
}