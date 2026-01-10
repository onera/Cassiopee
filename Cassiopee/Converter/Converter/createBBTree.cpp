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
// bbtree en python

# include "converter.h"
# include "Nuga/include/BbTree.h"

// Create a bbtree from a list of boxes mBB (min of bboxes), MBB (max of bboxes)
PyObject* K_CONVERTER::createBBTree(PyObject* self, PyObject* args)
{
    PyObject* pymBB;
    PyObject* pyMBB;
    if (!PYPARSETUPLE_(args, OO_, &pymBB, &pyMBB)) return NULL;

    E_Int nBB = PyList_Size(pyMBB);
    E_Float mBB[3];
    E_Float MBB[3];

    std::vector<K_SEARCH::BBox3D*> BBoxes(nBB);

    for (E_Int i = 0; i < nBB; i++)
    {
        PyObject* pyMin = PyList_GetItem(pymBB, i);
        PyObject* pyMax = PyList_GetItem(pyMBB, i);
        mBB[0] = PyFloat_AsDouble(PyList_GetItem(pyMin, 0));
        mBB[1] = PyFloat_AsDouble(PyList_GetItem(pyMin, 1));
        mBB[2] = PyFloat_AsDouble(PyList_GetItem(pyMin, 2));
        MBB[0] = PyFloat_AsDouble(PyList_GetItem(pyMax, 0));
        MBB[1] = PyFloat_AsDouble(PyList_GetItem(pyMax, 1));
        MBB[2] = PyFloat_AsDouble(PyList_GetItem(pyMax, 2));

        K_SEARCH::BBox3D* BBox = new K_SEARCH::BBox3D(mBB, MBB);
        BBoxes[i] = BBox;
    }

    K_SEARCH::BbTree3D* BBTree = new K_SEARCH::BbTree3D(BBoxes);

    PyObject* hook;
    void** packet = new void* [1];
    packet[0] = (void*)BBTree;

    #if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    hook = PyCObject_FromVoidPtr(packet, NULL);
    #else
    hook = PyCapsule_New(packet, NULL, NULL);
    #endif

    return hook;
}

// Return intersection of a box with hook stored bbtree
PyObject* K_CONVERTER::intersect(PyObject* self, PyObject* args)
{
    PyObject* pymBB; PyObject* pyMBB; PyObject* hook;
    if (!PYPARSETUPLE_(args, OOO_, &pymBB, &pyMBB, &hook)) return NULL;

    // Calcul de la BB de la zone
    E_Float mBB[3]; E_Float MBB[3];

    mBB[0] = PyFloat_AsDouble(PyList_GetItem(pymBB, 0));
    mBB[1] = PyFloat_AsDouble(PyList_GetItem(pymBB, 1));
    mBB[2] = PyFloat_AsDouble(PyList_GetItem(pymBB, 2));
    MBB[0] = PyFloat_AsDouble(PyList_GetItem(pyMBB, 0));
    MBB[1] = PyFloat_AsDouble(PyList_GetItem(pyMBB, 1));
    MBB[2] = PyFloat_AsDouble(PyList_GetItem(pyMBB, 2));

    K_SEARCH::BBox3D BBoxZone(mBB, MBB);
    
    // recupere le hook
    void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    packet = (void**) PyCObject_AsVoidPtr(hook);
#else
    packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
    
    K_SEARCH::BbTree3D* BBTree = (K_SEARCH::BbTree3D*)packet[0]; 

    //E_Int nBB = BBTree->_boxes.size();
    std::vector<E_Int> BBintersected; // vecteur d'indice des BB intersectees
    // recherche dans l'arbre les zones intersectees
    BBTree->getOverlappingBoxes(BBoxZone.minB, BBoxZone.maxB, BBintersected);

    E_Int nIntersect = BBintersected.size();

    PyObject* listOfIntersection = PyList_New(nIntersect);
    for (E_Int i = 0; i < nIntersect; i++)
    {
        PyList_SET_ITEM(listOfIntersection, i, PyLong_FromLong(BBintersected[i]));
    }

    return listOfIntersection;
}

// Return intersection of a list of boxes with hook stored bbtree
PyObject* K_CONVERTER::intersect2(PyObject* self, PyObject* args)
{
    PyObject* hook; PyObject* inBB;
    if (!PYPARSETUPLE_(args, OO_, &inBB, &hook)) return NULL;

    // recupere le hook
    void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    packet = (void**) PyCObject_AsVoidPtr(hook);
#else
    packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
    
    K_SEARCH::BbTree3D* BBTree = (K_SEARCH::BbTree3D*)packet[0]; 

    FldArrayF* f;
    K_NUMPY::getFromNumpyArray(inBB, f);
    E_Int size = f->getSize()/6;
    E_Float* fp = f->begin();

    //E_Int size = PyList_GET_SIZE(inBB);
    std::vector< std::vector<E_Int> > out(size);
    //printf("start %d\n", size); fflush(stdout);

#pragma omp parallel
{
    E_Float mBB[3];
    E_Float MBB[3];
    //PyObject* o;
    //PyObject* pymBB;
    //PyObject* pyMBB;
    std::vector<E_Int> BBintersected; // vecteur d'indice des BB intersectees
        
#pragma omp for
    for (E_Int i = 0; i < size; i++)
    {
        //o = PyList_GetItem(inBB, i);
        //pymBB = PyList_GetItem(o, 0);
        //pyMBB = PyList_GetItem(o, 1);
        //mBB[0] = PyFloat_AsDouble(PyList_GetItem(pymBB, 0));
        //mBB[1] = PyFloat_AsDouble(PyList_GetItem(pymBB, 1));
        //mBB[2] = PyFloat_AsDouble(PyList_GetItem(pymBB, 2));
        //MBB[0] = PyFloat_AsDouble(PyList_GetItem(pyMBB, 0));
        //MBB[1] = PyFloat_AsDouble(PyList_GetItem(pyMBB, 1));
        //MBB[2] = PyFloat_AsDouble(PyList_GetItem(pyMBB, 2));

        mBB[0] = fp[6*i+0];
        mBB[1] = fp[6*i+1];
        mBB[2] = fp[6*i+2];
        MBB[0] = fp[6*i+3];
        MBB[1] = fp[6*i+4];
        MBB[2] = fp[6*i+5];
        
        K_SEARCH::BBox3D BBoxZone(mBB, MBB);
        BBintersected.clear();

        // recherche dans l'arbre les zones intersectees
        BBTree->getOverlappingBoxes(BBoxZone.minB, BBoxZone.maxB, BBintersected);
        //E_Int nIntersect = BBintersected.size();
        out[i] = BBintersected; // copy
    }
}
    //printf("end\n"); fflush(stdout);
    
    PyObject* listOfIntersection = PyList_New(size);
    E_Int nbIntersect;
    for (E_Int i = 0; i < size; i++)
    {
        nbIntersect = out[i].size();
        PyObject* l = PyList_New(nbIntersect);
        PyList_SET_ITEM(listOfIntersection, i, l);
        
        for (E_Int j=0; j < nbIntersect; j++)
        {
            PyList_SET_ITEM(l, j, PyLong_FromLong(out[i][j]));
        }
    }
    RELEASESHAREDN(inBB, f);
    
    return listOfIntersection;
}

// Delete hook stored bbtree
PyObject* K_CONVERTER::deleteBBTree(PyObject* self, PyObject* args)
{
    PyObject* hook;
    if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;

    // recupere le hook
    void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    packet = (void**) PyCObject_AsVoidPtr(hook);
#else
    packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
    
    K_SEARCH::BbTree3D* BBTree = (K_SEARCH::BbTree3D*)packet[0]; 
    delete BBTree;
    Py_INCREF(Py_None);
    return Py_None;
}
