/*    
    Copyright 2013-2025 Onera.

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
#include "Karray.h"
#include "common/common.h"

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
        merr("Bad input array.");
        return 1;
    }

    if (ret == 1) {
        merr("Mesh should be an NGon.");
        Karray_free_structured(karray);
        return 1;
    }
    
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1) 
    {
        Karray_free_ngon(karray);
        merr("Coordinates not found.");
        return 1;
    }

    posx++; posy++; posz++;

    karray.x = karray.f->begin(posx);
    karray.y = karray.f->begin(posy);
    karray.z = karray.f->begin(posz);
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

    if (ret <= 0) 
    {
        merr("Bad input array");
        return 1;
    }

    if (ret != 1) 
    {
        merr("Mesh should be structured.");
        RELEASESHAREDB(ret, pyobject, karray.f, karray.cn);
        return 1;
    }

    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1) 
    {
        merr("Coordinates not found");
        Karray_free_structured(karray);
        return 1;
    }

    posx++; posy++; posz++;

    karray.x = karray.f->begin(posx);
    karray.y = karray.f->begin(posy);
    karray.z = karray.f->begin(posz);
    karray.npts = karray.f->getSize();

    karray.pyobject = pyobject;

    return 0;
}