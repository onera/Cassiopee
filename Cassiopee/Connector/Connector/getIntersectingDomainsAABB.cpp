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
# include "connector.h"
# include "Nuga/include/BbTree.h"

using namespace K_FUNC;
using namespace K_CONST;
using namespace K_FLD;
using namespace std;

// IN: liste de arrays de bounding boxes
// ============================================================================
/* Return max length of a cell */
// ============================================================================
PyObject* K_CONNECTOR::getIntersectingDomainsAABB(PyObject* self, PyObject* args)
{
    PyObject* arrays; E_Float tol;
    if (!PYPARSETUPLE_(args, O_ R_, &arrays, &tol)) return NULL;

    E_Int na = PyList_Size(arrays);
    vector<FldArrayF*> fs(na);
    for (E_Int i = 0; i < na; i++)
    {
      E_Int nil, njl, nkl;
      FldArrayF* f; FldArrayI* cn;
      char* varString; char* eltType;
      PyObject* array = PyList_GetItem(arrays, i);
      E_Int ret = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                                       cn, eltType);
      if (ret != 1) 
      {
        RELEASESHAREDB(ret, array, f, cn);
        PyErr_SetString(PyExc_TypeError,
                        "getIntersectingDomainsAABB: need bounding boxes arrays.");
        return NULL;
      }
      fs[i] = f;
    }

    // Cree le BbTree
    typedef K_SEARCH::BoundingBox<3> BBox3DType;
    std::vector<BBox3DType*> boxes(na); // a detruire a la fin
#pragma omp parallel default(shared)
    {
      E_Float xmin, ymin, zmin;
      E_Float xmax, ymax, zmax;
      E_Float minB0[3];  E_Float maxB0[3];
#pragma omp for
       for (E_Int et = 0; et < na; et++)
       {
         E_Float* x = fs[et]->begin(1);
         E_Float* y = fs[et]->begin(2);
         E_Float* z = fs[et]->begin(3);

         xmin = K_FUNC::E_min(x[0], x[1], x[2]);
         xmin = K_FUNC::E_min(xmin, x[3], x[4]);
         xmin = K_FUNC::E_min(xmin, x[5]);
         ymin = K_FUNC::E_min(y[0], y[1], y[2]);
         ymin = K_FUNC::E_min(ymin, y[3], y[4]);
         ymin = K_FUNC::E_min(ymin, y[5]);
         zmin = K_FUNC::E_min(z[0], z[1], z[2]);
         zmin = K_FUNC::E_min(zmin, z[3], z[4]);
         zmin = K_FUNC::E_min(zmin, z[5]);
         xmax = K_FUNC::E_max(x[0], x[1], x[2]);
         xmax = K_FUNC::E_max(xmax, x[3], x[4]);
         xmax = K_FUNC::E_max(xmax, x[5]);
         ymax = K_FUNC::E_max(y[0], y[1], y[2]);
         ymax = K_FUNC::E_max(ymax, y[3], y[4]);
         ymax = K_FUNC::E_max(ymax, y[5]);
         zmax = K_FUNC::E_max(z[0], z[1], z[2]);
         zmax = K_FUNC::E_max(zmax, z[3], z[4]);
         zmax = K_FUNC::E_max(zmax, z[5]);

         minB0[0] = xmin-tol; minB0[1] = ymin-tol; minB0[2] = zmin-tol;
         maxB0[0] = xmax+tol; maxB0[1] = ymax+tol; maxB0[2] = zmax+tol;
         boxes[et] = new BBox3DType(minB0, maxB0);
       }
    }
    K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);

    // Search
    std::vector<E_Int>* indicesBB = new std::vector<E_Int> [na];

#pragma omp parallel default(shared)
    {
        E_Float minB0[3];  E_Float maxB0[3];
        E_Float xmin, ymin, zmin;
        E_Float xmax, ymax, zmax;

#pragma omp for
    for (E_Int i = 0; i < na; i++)
    {
        E_Float* x = fs[i]->begin(1);
        E_Float* y = fs[i]->begin(2);
        E_Float* z = fs[i]->begin(3);

        xmin = K_FUNC::E_min(x[0], x[1], x[2]);
        xmin = K_FUNC::E_min(xmin, x[3], x[4]);
        xmin = K_FUNC::E_min(xmin, x[5]);
        ymin = K_FUNC::E_min(y[0], y[1], y[2]);
        ymin = K_FUNC::E_min(ymin, y[3], y[4]);
        ymin = K_FUNC::E_min(ymin, y[5]);
        zmin = K_FUNC::E_min(z[0], z[1], z[2]);
        zmin = K_FUNC::E_min(zmin, z[3], z[4]);
        zmin = K_FUNC::E_min(zmin, z[5]);
        xmax = K_FUNC::E_max(x[0], x[1], x[2]);
        xmax = K_FUNC::E_max(xmax, x[3], x[4]);
        xmax = K_FUNC::E_max(xmax, x[5]);
        ymax = K_FUNC::E_max(y[0], y[1], y[2]);
        ymax = K_FUNC::E_max(ymax, y[3], y[4]);
        ymax = K_FUNC::E_max(ymax, y[5]);
        zmax = K_FUNC::E_max(z[0], z[1], z[2]);
        zmax = K_FUNC::E_max(zmax, z[3], z[4]);
        zmax = K_FUNC::E_max(zmax, z[5]);
        minB0[0] = xmin-tol; minB0[1] = ymin-tol; minB0[2] = zmin-tol;
        maxB0[0] = xmax+tol; maxB0[1] = ymax+tol; maxB0[2] = zmax+tol; 
        bbtree->getOverlappingBoxes(minB0, maxB0, indicesBB[i]);
    }
    }

    PyObject* tpl = PyList_New(0);
    for (E_Int i = 0; i < na; i++)
    {
        PyObject* l = PyList_New(0);
        E_Int no;
        for (size_t j = 0; j < indicesBB[i].size(); j++)
        {
            no = indicesBB[i][j];
            if (no != i) PyList_Append(l, Py_BuildValue("l", no));
        }
        PyList_Append(tpl, l); Py_DECREF(l);
    }
    delete [] indicesBB;
    
    /*
    E_Float minB0[3];  E_Float maxB0[3];
    E_Float xmin, ymin, zmin;
    E_Float xmax, ymax, zmax;
    std::vector<E_Int> indicesBB;
    PyObject* tpl = PyList_New(0);
    for (E_Int i = 0; i < na; i++)
    {
        E_Float* x = fs[i]->begin(1);
        E_Float* y = fs[i]->begin(2);
        E_Float* z = fs[i]->begin(3);

        xmin = K_FUNC::E_min(x[0], x[1], x[2]);
        xmin = K_FUNC::E_min(xmin, x[3], x[4]);
        xmin = K_FUNC::E_min(xmin, x[5]);
        ymin = K_FUNC::E_min(y[0], y[1], y[2]);
        ymin = K_FUNC::E_min(ymin, y[3], y[4]);
        ymin = K_FUNC::E_min(ymin, y[5]);
        zmin = K_FUNC::E_min(z[0], z[1], z[2]);
        zmin = K_FUNC::E_min(zmin, z[3], z[4]);
        zmin = K_FUNC::E_min(zmin, z[5]);
        xmax = K_FUNC::E_max(x[0], x[1], x[2]);
        xmax = K_FUNC::E_max(xmax, x[3], x[4]);
        xmax = K_FUNC::E_max(xmax, x[5]);
        ymax = K_FUNC::E_max(y[0], y[1], y[2]);
        ymax = K_FUNC::E_max(ymax, y[3], y[4]);
        ymax = K_FUNC::E_max(ymax, y[5]);
        zmax = K_FUNC::E_max(z[0], z[1], z[2]);
        zmax = K_FUNC::E_max(zmax, z[3], z[4]);
        zmax = K_FUNC::E_max(zmax, z[5]);

        minB0[0] = xmin-tol; minB0[1] = ymin-tol; minB0[2] = zmin-tol;
        maxB0[0] = xmax+tol; maxB0[1] = ymax+tol; maxB0[2] = zmax+tol; 
        bbtree->getOverlappingBoxes(minB0, maxB0, indicesBB);
        PyObject* l = PyList_New(0);
        E_Int no;
        for (size_t j = 0; j < indicesBB.size(); j++)
        {
            no = indicesBB[j];
            if (no != i) PyList_Append(l, Py_BuildValue("l", no));
        }
        PyList_Append(tpl, l); Py_DECREF(l);
        indicesBB.clear();
    }
    */

    for (E_Int i = 0; i < na; i++) delete boxes[i];
    delete bbtree;

    for (E_Int i = 0; i < na; i++)
    {
        PyObject* array = PyList_GetItem(arrays, i);
        RELEASESHAREDS(array, fs[i]);
    }

    return tpl;
}
