/*    
    Copyright 2013-2019 Onera.

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
#include <string.h>
#include "Array/Array.h"
#include "String/kstring.h"

using namespace K_FLD;

//=====================================================================================
// Analyse an element string of the type "TRI_6" and return "TRI", nvpe, loc and typeId
//=====================================================================================
E_Int K_ARRAY::eltString2TypeId(char* eltString, char* eltType, E_Int& nvpe, 
  E_Int& loc, E_Int& typeId)
{
  E_Int ret = 1;
  char n[128];
  E_Int l = strlen(eltString);
  if (eltString[l-1] == '*') { loc = 1; l = l-1; }
  else loc = 0;
  typeId = 0;
  
  if (K_STRING::cmp(eltString, 5, "TETRA") == 0)
  {
    strncpy(eltType, eltString, 5); eltType[5] = '\0';
    if (l == 5) nvpe = 4;
    else
    {
      strncpy(n, eltString+6, l-6); n[l-6] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 4:
        typeId = 4; break;
      case 10:
        typeId = 14; break;
      case 16:
        typeId = 35; break;
      case 20:
        typeId = 36; break;
      default:
        ret = 0; break;
    }
  }
  else if (K_STRING::cmp(eltString, 4, "NODE") == 0)
  {
    strncpy(eltType, eltString, 4); eltType[4] = '\0';
    nvpe = 1;
    typeId = 0;
  }
  else if (K_STRING::cmp(eltString, 4, "QUAD") == 0)
  {
    if (l == 4) nvpe = 4;
    else
    {
      strncpy(eltType, eltString, 4); eltType[4] = '\0';
      strncpy(n, eltString+5, l-5); n[l-5] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 4:
        typeId = 3; break;
      case 8:
        typeId = 12; break;
      case 9:
        typeId = 13; break;
      case 12:
        typeId = 33; break;
      case 16:
        typeId = 34; break;
      case 25:
        typeId = 54; break;
      default:
        ret = 0; break;
    }
  }
  else if (K_STRING::cmp(eltString, 4, "HEXA") == 0)
  {
    if (l == 4) nvpe = 8;
    else
    {
      strncpy(eltType, eltString, 4); eltType[4] = '\0';
      strncpy(n, eltString+5, l-5); n[l-5] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 8:
        typeId = 7; break;
      case 20:
        typeId = 18; break;
      case 27:
        typeId = 19; break;
      case 32:
        typeId = 40; break;
      case 56:
        typeId = 41; break;
      case 64:
        typeId = 42; break;
      default:
        ret = 0; break;  
    }
  }
  else if (K_STRING::cmp(eltString, 4, "PYRA") == 0)
  {
    if (l == 4) nvpe = 5;
    else
    {
      strncpy(eltType, eltString, 4); eltType[4] = '\0';
      strncpy(n, eltString+5, l-5); n[l-5] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 5:
        typeId = 5; break;
      case 14:
        typeId = 15; break;
      case 13:
        typeId = 20; break;
      case 21:
        typeId = 37; break;
      case 29:
        typeId = 38; break;
      case 30:
        typeId = 39; break;
      default:
        ret = 0; break;   
    }
  }
  else if (K_STRING::cmp(eltString, 5, "PENTA") == 0)
  {
    if (l == 5) nvpe = 6;
    else
    {
      strncpy(eltType, eltString, 5); eltType[5] = '\0';
      strncpy(n, eltString+6, l-6); n[l-6] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 6:
        typeId = 6; break;
      case 15:
        typeId = 16; break;
      case 18:
        typeId = 17; break;
      default:
        ret = 0; break;
    }
  }
  else if (K_STRING::cmp(eltString, 4, "NGON") == 0)
  {
    strncpy(eltType, eltString, 4); eltType[4] = '\0';
    nvpe = 0;
    typeId = 8;
  }
  else if (K_STRING::cmp(eltString, 3, "TRI") == 0)
  {
    if (l == 3) nvpe = 3;
    else
    {
      strncpy(eltType, eltString, 3); eltType[3] = '\0';
      strncpy(n, eltString+4, l-4); n[l-4] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 3:
        typeId = 2; break;
      case 6:
        typeId = 11; break;
      case 9:
        typeId = 31; break;
      case 10:
        typeId = 32; break;
      case 12:
        typeId = 52; break;
      case 15:
        typeId = 53; break;
      default:
        ret = 0; break;  
    }
  }
  else if (K_STRING::cmp(eltString, 3, "BAR") == 0)
  {
    if (l == 3) nvpe = 2;
    else
    {
      strncpy(eltType, eltString, 3); eltType[3] = '\0';
      strncpy(n, eltString+4, l-4); n[l-4] = '\0';
      nvpe = atoi(n);
    }
    switch (nvpe)
    {
      case 2:
        typeId = 1; break;
      case 3:
        typeId = 10; break;
      case 4:
        typeId = 30; break;
      case 5:
        typeId = 51; break;
      default:
        ret = 0; break;    
    }
  }
  return ret;
}

// Analyse a typeId and loc and return eltString and nvpe
E_Int K_ARRAY::typeId2eltString(E_Int typeId, E_Int loc, char* eltString, E_Int& nvpe)
{  
  E_Int ret = 1;
  switch (typeId)
  {
    case 0:
      strcpy(eltString, "NODE"); nvpe = 1; break;
    
    case 1:
      strcpy(eltString, "BAR"); nvpe = 2; break;
      
    case 2:
      strcpy(eltString, "TRI"); nvpe = 3; break;
    
    case 3:
      strcpy(eltString, "QUAD"); nvpe = 4; break;
    
    case 4:
      strcpy(eltString, "TETRA"); nvpe = 4; break;
    
    case 5:
      strcpy(eltString, "PYRA"); nvpe = 5; break;
    
    case 6:
      strcpy(eltString, "PENTA"); nvpe = 6; break;
  
    case 7:
      strcpy(eltString, "HEXA"); nvpe = 8; break;
      
    case 8:
      strcpy(eltString, "NGON"); nvpe = 1; break;
      
    case 10:
      strcpy(eltString, "BAR_3"); nvpe = 3; break;
      
    case 11: 
      strcpy(eltString, "TRI_6"); nvpe = 6; break;
      
    case 12:
      strcpy(eltString, "QUAD_8"); nvpe = 8; break;
      
    case 13:
      strcpy(eltString, "QUAD_9"); nvpe = 9; break;
    
    case 14:
      strcpy(eltString, "TETRA_10"); nvpe = 10; break;
    
    case 15:
      strcpy(eltString, "PYRA_14"); nvpe = 14; break;
    
    case 16:
      strcpy(eltString, "PENTA_15"); nvpe = 15; break;
    
    case 17:
      strcpy(eltString, "PENTA_18"); nvpe = 18; break;
    
    case 18:
      strcpy(eltString, "HEXA_20"); nvpe = 20; break;
    
    case 19:
      strcpy(eltString, "HEXA_27"); nvpe = 27; break;
    
    case 20:
      strcpy(eltString, "PYRA_13"); nvpe = 13; break;
    
    case 30:
      strcpy(eltString, "BAR_4"); nvpe = 4; break;
    
    case 31:
      strcpy(eltString, "TRI_9"); nvpe = 9; break;
    
    case 32:
      strcpy(eltString, "TRI_10"); nvpe = 10; break;
    
    case 33:
      strcpy(eltString, "QUAD_12"); nvpe = 12; break;
    
    case 34:
      strcpy(eltString, "QUAD_16"); nvpe = 16; break;
    
    case 35:
      strcpy(eltString, "TETRA_16"); nvpe = 16; break;
    
    case 36:
      strcpy(eltString, "TETRA_20"); nvpe = 20; break;
    
    case 37:
      strcpy(eltString, "PYRA_21"); nvpe = 21; break;
    
    case 38:
      strcpy(eltString, "PYRA_29"); nvpe = 29; break;
    
    case 39:
      strcpy(eltString, "PYRA_30"); nvpe = 30; break;
    
    case 40:
      strcpy(eltString, "HEXA_32"); nvpe = 32; break;
    
    case 41:
      strcpy(eltString, "HEXA_56"); nvpe = 56; break;
    
    case 42:
      strcpy(eltString, "HEXA_64"); nvpe = 64; break;
      
    default:
      ret = 0; break;
  }
          
  if (loc == 1) strcat(eltString, "*");
  
  return ret;
}
//=============================================================================
// Extrait les donnees (shared) d'un objet python struct array
// defini par: Array1: [ 'vars', a, ni, nj, nk ]
//             Array2: [ 'vars', [a], ni, nj, nk ]
// ou d'un objet python unstruct array
// defini par: Array1: [ 'vars', a, c, "ELTTYPE"]
//             Array2: [ 'vars', [a], [c], "ELTTYPE"]
// ou ELTTYPE vaut: NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON
// avec ou sans star.
// return 1: valid struct array
//           f: field en stockage Fld (champs par variable)
//           ni, nj, nk: number of points
//           varString
// return 2: valid unstruct array
//           f: field en stockage Fld (champs par variable)
//           c: connectivity (champs des indices de noeuds commencants a 1)
//           eltType: type of elements
//           varString
// Dans tous les cas sauf NGON: c est la connectivite elt-noeuds 
// dimensionnee (elt, nbre de noeuds par elements)
// Dans le cas NGON: c est la connectivite Face-noeuds, suivi de la 
// connectivite elt-faces + 4 entiers indiquant les tailles
//
// return -1: given object is not a list.
// return -2: not a valid number of elts in list.
// return -3: first element is not a var string.
// return -4: a is not a valid numpy array.
// return -5: array is structured but ni, nj, nk unvalid.
// return -6: array is unstructured but connectivity is unvalid.
// return -7: array is unstructured but elt type is unknown.
// C'est la responsabilite de l'appelant de liberer la memoire de f et 
// eventuellement de c
//=============================================================================
E_Int K_ARRAY::getFromArray2(PyObject* o,
                             char*& varString,
                             FldArrayF*& f,
                             E_Int& ni, E_Int& nj, E_Int& nk,
                             FldArrayI*& c,
                             char*& eltType)
{
  PyObject* tpl; PyObject* p;
  PyArrayObject* a; PyArrayObject* ac;
  PyObject* ref; //PyObject* ref2;
  IMPORTNUMPY;

  // -- list --
  if (PyList_Check(o) == false)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
    return -1;
  }
  E_Int size = PyList_Size(o);
  if (size != 4 && size != 5)
  {
    PyErr_Warn(PyExc_Warning, 
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check number of elements in list.");
    return -2;
  }
 
  // -- varString --
  PyObject* l = PyList_GetItem(o,0);
  if (PyString_Check(l))
  {
    // pointeur sur la chaine python
    varString = PyString_AsString(PyList_GetItem(o,0));
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(l))
  {
    varString = PyBytes_AsString(PyUnicode_AsUTF8String(l)); 
  }
#endif
  else
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. First element must be a string.");
    return -3;
  }
  
  E_Int nvar = getNumberOfVariables(varString);

  // -- field --
  tpl = PyList_GetItem(o, 1);
  //ref = tpl; Py_INCREF(ref); // trick
  Py_INCREF(tpl); ref = tpl;

  if (PyArray_Check(tpl) == true) // -- Array1 --
  {
    a = (PyArrayObject*)tpl;

    if (PyArray_NDIM(a) != 2)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: field must have two dimensions.");
      Py_DECREF(tpl);
      return -4;
    }
  
    E_Int s = PyArray_DIMS(a)[1];
    E_Int nfld = PyArray_DIMS(a)[0];

    if (nfld != nvar)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: number of variables different in varString and field.");
      Py_DECREF(tpl);
      return -4;
    }    
    f = new FldArrayF(s, nfld, (E_Float*)PyArray_DATA(a), true, true);
  }
  else if (PyList_Check(tpl) == true) // -- Array2 --
  {
    E_Int nfld = PyList_Size(tpl);
    E_Float** acu = new E_Float* [nfld];
    E_Int s = 0;
    for (E_Int i = 0; i < nfld; i++)
    {
      p = PyList_GetItem(tpl, i);
      if (PyArray_Check(p) == true)
      {
        a = (PyArrayObject*)p; //Py_INCREF(a);
        s = PyArray_SIZE(a); // must be constant
        acu[i] = (E_Float*)PyArray_DATA(a);
      }
      else
      {
        PyErr_Warn(PyExc_Warning,
                   "getFromArray: field must a list of numpys.");
        Py_DECREF(tpl);
        return -4;
      }
    }
    f = new FldArrayF(s, nfld, acu, true, true);
    delete [] acu;
  }
  else
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: second arg in array must be a numpy array.");
    return -4;
  }
  
  if (size == 4) // unstruct array
  {
    // -- connectivity --
    tpl = PyList_GetItem(o, 2);
    //ref2 = tpl; Py_INCREF(ref2);
    Py_INCREF(tpl);
    if (PyArray_Check(tpl) == true) // -- Array1 --
    {
      ac = (PyArrayObject*)tpl;
  
      if (PyArray_NDIM(ac) > 2)
      {
        PyErr_Warn(PyExc_Warning,
                   "getFromArray: connectivity must have two dimensions max.");
        Py_DECREF(ref); Py_DECREF(tpl);
        return -4;
      }
      E_Int s, nfld;
      if (PyArray_NDIM(ac) == 2)
      { s = PyArray_DIMS(ac)[1]; nfld = PyArray_DIMS(ac)[0]; }
      else { s = PyArray_DIMS(ac)[0]; nfld = 1; }
      c = new FldArrayI(s, nfld, (E_Int*)PyArray_DATA(ac), true);
    }
    else if (PyList_Check(tpl) == true) // -- Array2 --
    {
      E_Int nc = PyList_Size(tpl);
      if (nc == 1) // BE => compact + stride
      {
        ac = (PyArrayObject*)PyList_GetItem(tpl,0);
        E_Int s, nfld;
        if (PyArray_NDIM(ac) == 2)
        { s = PyArray_DIMS(ac)[0]; nfld = PyArray_DIMS(ac)[1]; }
        else { s = PyArray_DIMS(ac)[0]; nfld = 1; }
        c = new FldArrayI(s, nfld, (E_Int*)PyArray_DATA(ac), true, false);
      }
      else // suppose NGON
      {
        if (nc == 4) // NGON/NFACE/indPG/indPH
        {
          PyArrayObject* ac1 = (PyArrayObject*)PyList_GetItem(tpl,0); // NGON
          PyArrayObject* ac2 = (PyArrayObject*)PyList_GetItem(tpl,1); // NFACE
          PyArrayObject* ac3 = (PyArrayObject*)PyList_GetItem(tpl,2); // indPG
          PyArrayObject* ac4 = (PyArrayObject*)PyList_GetItem(tpl,3); // indPH
          E_Int nfaces = PyArray_SIZE(ac3);
          E_Int nelts = PyArray_SIZE(ac4);
          E_Int sizeNGon = PyArray_SIZE(ac1);
          E_Int sizeNFace = PyArray_SIZE(ac2);
          c = new FldArrayI(nfaces, nelts, (E_Int*)PyArray_DATA(ac1), (E_Int*)PyArray_DATA(ac2),
                            (E_Int*)PyArray_DATA(ac3), (E_Int*)PyArray_DATA(ac4), sizeNGon, sizeNFace);
        }
      }
    }
    else
    {
      PyErr_Warn(PyExc_Warning, 
                 "getFromArray: third arg in array must be a numpy array.");
      Py_DECREF(ref); Py_DECREF(tpl);
      return -6;
    }

    // -- element type --
    PyObject* l = PyList_GetItem(o,3);
    if (PyString_Check(l))
    {
      eltType = PyString_AsString(l);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l))
    {
      eltType = PyBytes_AsString(PyUnicode_AsUTF8String(l)); 
    }
#endif
    else
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: an unstruct array must be of list of type ['vars', a, c, 'ELTTYPE']. Last element must be a string.");
      Py_DECREF(ref); Py_DECREF(tpl);
      return -7;
    }
    
    char st[256]; E_Int dummy;
    if (eltString2TypeId(eltType, st, dummy, dummy, dummy) == 0)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: element type unknown: %s. Must be in NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON or NODE*, BAR*, TRI*, QUAD*, TETRA*, PYRA*, PENTA*, HEXA*, NGON*.");
      Py_DECREF(ref); Py_DECREF(tpl);
      return -7;
    }
    return 2;
  }
  else // struct array
  {
    tpl = PyList_GetItem(o,2);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: third arg must be an integer.");
      Py_DECREF(ref);
      return -5;
    }
    ni = PyLong_AsLong(tpl);
    tpl = PyList_GetItem(o,3);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: fourth arg must be an integer.");
      Py_DECREF(ref);
      return -5;
    }
    nj = PyLong_AsLong(tpl);
    tpl = PyList_GetItem(o,4);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: fifth arg must be an integer.");
      Py_DECREF(ref);
      return -5;
    }
    nk = PyLong_AsLong(tpl);
    return 1;
  }
}

// Extrait sans ni,nj,nk
E_Int K_ARRAY::getFromArray2(PyObject* o,
                             char*& varString,
                             FldArrayF*& f,
                             FldArrayI*& c,
                             char*& eltType)
{
  E_Int ni, nj, nk;
  E_Int ret = getFromArray2(o, varString, f,
                            ni, nj, nk, c, eltType);
  return ret;
}

// Extrait sans ni,nj,nk,varString,eltType
E_Int K_ARRAY::getFromArray2(PyObject* o,
                             FldArrayF*& f,
                             FldArrayI*& c)
{
  E_Int ni, nj, nk; char* varString; char* eltType;
  E_Int ret = getFromArray2(o, varString, f,
                            ni, nj, nk, c, eltType);
  return ret;
}

//=============================================================================
// Extrait les donnees (shared) d'un objet python struct array
// defini par: Array1: [ 'vars', a, ni, nj, nk ]
//             Array2: [ 'vars', [a], ni, nj, nk ]
// ou d'un objet python unstruct array
// defini par: Array1: [ 'vars', a, c, "ELTTYPE" ]
//             Array2: [ 'vars', [a], [c], "ELTTYPE" ]
// ou ELTTYPE vaut: NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON
// avec ou sans star.
// Ne retourne que les champs et la varstring
// Retourne 1 si ok.
//=============================================================================
E_Int K_ARRAY::getFromArray2(PyObject* o,
                             char*& varString,
                             FldArrayF*& f)
{
  PyObject* tpl; PyObject* p;
  PyArrayObject* a;
  PyObject* ref; //PyObject* ref2;
  IMPORTNUMPY;

  // -- list --
  if (PyList_Check(o) == false)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
    return -1;
  }
  E_Int size = PyList_Size(o);
  if (size != 4 && size != 5)
  {
    PyErr_Warn(PyExc_Warning, 
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check number of elements in list.");
    return -2;
  }
 
  // -- varString --
  PyObject* l = PyList_GetItem(o,0);
  if (PyString_Check(l))
  {
    // pointeur sur la chaine python
    varString = PyString_AsString(l);
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(l))
  {
    varString = PyBytes_AsString(PyUnicode_AsUTF8String(l)); 
  }
#endif
  else
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. First element must be a string.");
    return -3;
  }
  E_Int nvar = getNumberOfVariables(varString);

  // -- field --
  tpl = PyList_GetItem(o, 1);
  //ref = tpl; Py_INCREF(ref); // trick
  Py_INCREF(tpl); ref = tpl;

  if (PyArray_Check(tpl) == true) // -- Array1 --
  {
    a = (PyArrayObject*)tpl;

    if (PyArray_NDIM(a) != 2)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: field must have two dimensions.");
      Py_DECREF(tpl);
      return -4;
    }
  
    E_Int s = PyArray_DIMS(a)[1];
    E_Int nfld = PyArray_DIMS(a)[0];

    if (nfld != nvar)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: number of variables different in varString and field.");
      Py_DECREF(tpl);
      return -4;
    }    
    f = new FldArrayF(s, nfld, (E_Float*)PyArray_DATA(a), true, true);
  }
  else if (PyList_Check(tpl) == true) // -- Array2 --
  {
    E_Int nfld = PyList_Size(tpl);
    E_Float** acu = new E_Float* [nfld];
    E_Int s = 0;
    for (E_Int i = 0; i < nfld; i++)
    {
      p = PyList_GetItem(tpl, i);
      if (PyArray_Check(p) == true)
      {
        a = (PyArrayObject*)p; //Py_INCREF(a);
        s = PyArray_SIZE(a); // must be constant
        acu[i] = (E_Float*)PyArray_DATA(a);
      }
      else
      {
        PyErr_Warn(PyExc_Warning,
                   "getFromArray: field must a list of numpys.");
        Py_DECREF(tpl);
        return -4;
      }
    }
    f = new FldArrayF(s, nfld, acu, true, true);
    delete [] acu;
  }
  else
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: second arg in array must be a numpy array.");
    return -4;
  }
  
  return 1;
}
// Extrait uniquement les champs
E_Int K_ARRAY::getFromArray2(PyObject* o,
                             FldArrayF*& f)
{
  char* varString;
  E_Int ret = getFromArray2(o, varString, f);
  return ret;
}
