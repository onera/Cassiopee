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
#include "stream_line.hpp"
#include "structured_data_view.hpp"
#include "unstructured_data_view.hpp"
#include "ngon_data_view.hpp"

PyObject* K_POST::comp_stream_line(PyObject* self, PyObject* args)
{
    PyObject* arrays;
    PyObject* surfArray;
    PyObject* listOfPoints;
    E_Float   x0, y0, z0;
    PyObject* vectorNames;
    E_Int     nStreamPtsMax;
    E_Int     signe;

    if (!PYPARSETUPLE_(args, OOOO_ II_, &arrays, &surfArray, &listOfPoints,
                      &vectorNames, &signe, &nStreamPtsMax)) {
        return NULL;
    }
    std::vector<point3d> beg_nodes;
    if (!PyList_Check(listOfPoints))
    {
        if (!PyTuple_Check(listOfPoints))
        {
            PyErr_SetString(PyExc_TypeError, "streamLine: third argument must be a point or a list of point [(x0,y0,z0),(x1,y1,z1),...]");
            return NULL;
        }
        if (! PYPARSETUPLE_(listOfPoints, RRR_, &x0, &y0, &z0) )
        {
            PyErr_SetString(PyExc_TypeError, "streamLine: a point must be a triplet of doubles");
            return NULL;
        }
        beg_nodes.emplace_back(point3d{x0,y0,z0});
    }
    else
    {
        E_Int nb_points = PyList_Size(listOfPoints);
        //std::cout << "Nombre de points : " << nb_points << std::endl;
        beg_nodes.reserve(nb_points);
        for (E_Int i = 0; i < nb_points; ++i)
        {
            PyObject* tple = PyList_GetItem(listOfPoints, i);
            if (!PyTuple_Check(tple))
            {
                PyErr_SetString(PyExc_TypeError, "streamLine: a point inside the list of nodes must be a tuple of reals");
                return NULL;
            }
            //std::cout << "on depiote tple at " << (void*)tple << std::flush << std::endl;
            if (! PYPARSETUPLE_(tple, RRR_, &x0, &y0, &z0) )
            {
                PyErr_SetString(PyExc_TypeError, "streamLine: a point must be a triplet of reals");
                return NULL;
            }
            //std::cout << "Rajout point " << std::string(point3d{x0,y0,z0}) << std::flush << std::endl;
            beg_nodes.emplace_back(point3d{x0,y0,z0});
        }
    }
    // Check every array in arrays
    if (PyList_Check(arrays) == 0) {
        PyErr_SetString(PyExc_TypeError, "streamLine: first argument must be a list.");
        return NULL;
    }
    // extraction of the 3 components of the vector used in the streamline computation
    std::vector<char*> vnames;
    if (PyList_Check(vectorNames) == 0) {
        PyErr_SetString(PyExc_TypeError,
                        "streamLine: three variable inside a liste must be defined as component of the streamline vector.");
        return NULL;
    }
    E_Int sizeOfVector = PyList_Size(vectorNames);
    if (sizeOfVector != 3) {
        PyErr_SetString(PyExc_TypeError, "streamLine: vector must be defined by 3 components.");
        return NULL;
    }
    for (E_Int i = 0; i < PyList_Size(vectorNames); i++) {
        PyObject* tpl0 = PyList_GetItem(vectorNames, i);
        if (PyString_Check(tpl0)) {
            char* str = PyString_AsString(tpl0);
            vnames.push_back(str);
        }
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(tpl0)) {
            char* str = (char*)PyUnicode_AsUTF8(tpl0);
            vnames.push_back(str);
        }
#endif
        else {
            PyErr_SetString(PyExc_TypeError, "streamLine: vector component name must be a string.");
            return NULL;
        }
    }
    // Construction de la vue de chaque zone :
    // Extract infos from arrays
    std::vector<E_Int> resl;
    std::vector<char*> structVarString; std::vector<char*> unstrVarString;
    std::vector<FldArrayF*> structF; std::vector<FldArrayF*> unstrF;
    std::vector<E_Int> nit; std::vector<E_Int> njt; std::vector<E_Int> nkt;
    std::vector<FldArrayI*> cnt;
    std::vector<char*> eltType;
    std::vector<PyObject*> objs, obju;
    E_Boolean skipNoCoord = true;
    E_Boolean skipStructured = false;
    E_Boolean skipUnstructured = false; 
    E_Boolean skipDiffVars = true;
    //E_Int isOk = 
    K_ARRAY::getFromArrays( arrays, resl, 
                            structVarString, unstrVarString,
                            structF, unstrF, nit, njt, nkt, 
                            cnt, eltType, objs, obju, skipDiffVars,
                            skipNoCoord, skipStructured, 
                            skipUnstructured, true);

    char* varStringOut;
    if (structVarString.size() > 0)
    {
        varStringOut = new char [strlen(structVarString[0])+1];
        strcpy(varStringOut, structVarString[0]);
    }
    else if (unstrVarString.size() > 0)
    {
        varStringOut = new char [strlen(unstrVarString[0])+1];
        strcpy(varStringOut, unstrVarString[0]);
    }
    else
    {
        varStringOut = new char [2];
        varStringOut[0] = '\0';
    }
    E_Int nb_zones = structVarString.size() + unstrVarString.size();
    if (nb_zones == 0)
    {
        PyObject* tpl = PyList_New(0);
        return tpl;
    }

    // Build interpData 
    E_Int nzonesS = structF.size();
    E_Int nzonesU = unstrF.size();
    std::vector<E_Int> posxs1; std::vector<E_Int> posys1; std::vector<E_Int> poszs1; std::vector<E_Int> poscs1;
    std::vector<FldArrayF*> structF1;
    std::vector<E_Int> nis1; std::vector<E_Int> njs1; std::vector<E_Int> nks1;
    std::vector<char*> structVarStrings1;
    //E_Int isBuilt;
    std::vector<zone_data_view> zones;
    zones.reserve(nb_zones);
    for (E_Int no = 0; no < nzonesS; no++)
    {
        E_Int posx = K_ARRAY::isCoordinateXPresent(structVarString[no]); posx++;
        E_Int posy = K_ARRAY::isCoordinateYPresent(structVarString[no]); posy++;
        E_Int posz = K_ARRAY::isCoordinateZPresent(structVarString[no]); posz++;
        E_Int posc = K_ARRAY::isCellNatureField2Present(structVarString[no]); posc++;
        // On extrait la position du vecteur servant au streamline (a priori la vitesse!) du champs structF :
        E_Int posv1 = K_ARRAY::isNamePresent(vnames[0], structVarString[no]);
        E_Int posv2 = K_ARRAY::isNamePresent(vnames[1], structVarString[no]);
        E_Int posv3 = K_ARRAY::isNamePresent(vnames[2], structVarString[no]);
        // variables pas presentes dans l'array v
        if (posv1 == -1 || posv2 == -1 || posv3 == -1)
        {
            for (unsigned int nos = 0; nos < objs.size(); nos++) RELEASESHAREDS(objs[nos], structF[nos]);
            for (unsigned int nos = 0; nos < obju.size(); nos++) RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
            PyErr_SetString(PyExc_TypeError, 
                    "streamLine: Vector is missing.");
            return NULL;
        }
        posv1++; posv2++; posv3++;

        zones.push_back( structured_data_view( {nit[no],njt[no],nkt[no]}, structF[no], {posx, posy, posz}, 
                                               {posv1, posv2, posv3}, posc));
    }

    // InterpData non structuree (pour les n-gons, c'est aussi ici ?)
    //std::cout << "Construction des vues pour zones non structurees (" << nzonesU << ")." << std::flush << std::endl;
    for (E_Int no = 0; no < nzonesU; no++)
    {
      E_Int posx = K_ARRAY::isCoordinateXPresent(unstrVarString[no]); posx++;
      E_Int posy = K_ARRAY::isCoordinateYPresent(unstrVarString[no]); posy++;
      E_Int posz = K_ARRAY::isCoordinateZPresent(unstrVarString[no]); posz++;
      E_Int posc = K_ARRAY::isCellNatureField2Present(unstrVarString[no]); posc++;
      // On extrait la position du vecteur servant au streamline (a priori la vitesse!) du champs structF :
      E_Int posv1 = K_ARRAY::isNamePresent(vnames[0], unstrVarString[no]);
      E_Int posv2 = K_ARRAY::isNamePresent(vnames[1], unstrVarString[no]);
      E_Int posv3 = K_ARRAY::isNamePresent(vnames[2], unstrVarString[no]);
      // variables pas presentes dans l'array v
      if (posv1== -1 || posv2 == -1 || posv3 == -1)
      {
        for (unsigned int nos = 0; nos < objs.size(); nos++) RELEASESHAREDS(objs[nos], structF[nos]);
            for (unsigned int nos = 0; nos < obju.size(); nos++) RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
        PyErr_SetString(PyExc_TypeError, 
                "streamLine: Vector is missing.");
        return NULL;
      }
      posv1++; posv2++; posv3++;

      // On separe la connectivité par type d'elements :
      // cnt[no] => FldArrayI -> connectivité element vers sommets, je pense...
      // eltType[no] : Type d'élément
      // Vérifier si c'est un n-gon ou pas :
      if (K_STRING::cmp(eltType[no], "NGON") != 0 && K_STRING::cmp(eltType[no], "NGON*") != 0)
        zones.push_back( unstructured_data_view( eltType[no], cnt[no], unstrF[no], {posx, posy, posz}, {posv1,posv2,posv3}) );
      else
      { // Pour le n-gon, cnt[no] a une autre signification et possède le format suivant :
        // [ nb_faces, nb_entiers_pour_connectivité_face_vers_sommets, nb_sommets_face_1, index_sommet_1_1, index_sommet_1_2, ...,
        //   nb_sommets_face_2, ...., index_sommet_n_p,
        //   nb_elts, nb_entiers_pour_connectivité_elt_vers_faces, nbre_face_elt1, face1, face2, ..., nbre_face_elt2, face1, face2,
        //   ..., facep]
        E_Int nfaces = (*cnt[no])[0];
        E_Int ndata_faces = (*cnt[no])[1];
        E_Int nelts  = (*cnt[no])[2+ndata_faces];
        //E_Int ndata_elts  = (*cnt[no])[3+ndata_faces];
        zones.push_back( ngon_data_view( nfaces, {cnt[no]->begin()+2,cnt[no]->begin()+2+ndata_faces}, 
                                         nelts , {cnt[no]->begin()+4+ndata_faces, cnt[no]->end()}, 
                                         unstrF[no], {posx, posy, posz}, {posv1,posv2,posv3}) );
      }// Fin cas traitement ngon
    }// Fin boucle sur les zones non structurées
    //std::cout << "Je vais calculer une streamline" << std::flush << std::endl;
    /*
    if (beg_nodes.size() == 1)
    {
        PyObject* list_of_streamlines = PyList_New(1);
        streamline sline( {x0,y0,z0}, zones, nStreamPtsMax, (signe==2) );
        FldArrayF& field = sline.field();
        E_Int number_of_points = field.getSize();
        PyObject* tpl = K_ARRAY::buildArray(field, varStringOut, number_of_points, 1, 1);
        //delete [] varStringOut;
        PyList_SetItem(list_of_streamlines, 0, tpl);
        return list_of_streamlines;
    }
    */
    PyObject* list_of_streamlines = PyList_New(beg_nodes.size());
    std::vector<FldArrayF> arrFields(beg_nodes.size());
#   pragma omp parallel for schedule(dynamic,10)
    for (size_t i = 0; i < beg_nodes.size(); ++i)
    {
        //#pragma omp critical
        //std::cout << "Calcul streamline no" << i+1 << std::flush << std::endl;
        try
        {
            streamline sline(beg_nodes[i], zones, nStreamPtsMax, (signe==2));
            FldArrayF& field = sline.field();
            arrFields[i] = field;
        }
        catch (std::exception& err)
        {
            printf("Warning: streamLine: %s\n", err.what());
            Py_INCREF(Py_None);
            PyList_SetItem(list_of_streamlines, i, Py_None);
        }
    }
    for (size_t i = 0; i < beg_nodes.size(); ++i)
    {
        E_Int number_of_points = arrFields[i].getSize();
        /* Essai pour supprimer les streams avec 0 points */
        if (number_of_points > 0)
        {
            PyObject* tpl = K_ARRAY::buildArray(arrFields[i], varStringOut, number_of_points, 1, 1);
            PyList_SetItem(list_of_streamlines, i, tpl);
        }
        else 
        {
            Py_INCREF(Py_None);
            PyList_SetItem(list_of_streamlines, i, Py_None);
        }
    }
    
    // Compact - Essai pour enlever des streamlines qui auraient 0 points
    /*
    size_t size = 0;
    for (size_t i = 0; i < beg_nodes.size(); ++i)
    {
        PyObject* tpl = PyList_GetItem(list_of_streamlines, i);
        if (tpl != Py_None) size += 1; 
    }
    if (size < beg_nodes.size())
    {
        PyObject* list_of_streamlines2 = PyList_New(size);
        E_Int j = 0;
        for (size_t i = 0; i < beg_nodes.size(); ++i)
        {
            PyObject* tpl = PyList_GetItem(list_of_streamlines, i);
            if (tpl != Py_None) { PyList_SetItem(list_of_streamlines2, j, tpl); j += 1; }
        }
        Py_DECREF(list_of_streamlines);
        list_of_streamlines = list_of_streamlines2;
    }
    */
   
    delete [] varStringOut;
    for (unsigned int nos = 0; nos < objs.size(); nos++)
        RELEASESHAREDS(objs[nos], structF[nos]);
    for (unsigned int nos = 0; nos < obju.size(); nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

    return list_of_streamlines;
}
