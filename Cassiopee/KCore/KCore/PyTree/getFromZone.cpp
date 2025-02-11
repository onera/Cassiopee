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
#include "PyTree/PyTree.h"
#include "String/kstring.h"
using namespace K_FLD;
using namespace std; 

//==============================================================================
// Ajoute variable variable a la varString
// IN/OUT: varString
// IN/OUT: size: taille allouee de varString
// IN: variable: variable a ajouter 
// Si varString n'existe pas, il faut envoyer varString=NULL
//==============================================================================
void K_PYTREE::add2VarString(char*& varString, E_Int& size, 
                             const char* variable)
{
  E_Int l1, l2;
  char* s;
  l2 = strlen(variable);
  if (varString == NULL)
  {
    size = l2+512;
    s = new char [size+1];
    strcpy(s, variable);
    varString = s;
    return;
  }
  l1 = strlen(varString);
  if (l1+l2 >= size)
  {
    size += l2+512;
    s = new char [size+1];
    strcpy(s, varString);
    strcat(s, ",");
    strcat(s, variable);
    delete [] varString;
    varString = s;
    return;
  }
  else
  {
    strcat(varString, ",");
    strcat(varString, variable);
    return;
  }
}

//=============================================================================
// Extrait les donnees d'une zone d'un pyTree
// Il s'agit de donnees partagee avec python.
// Les donnees dont donc exactement les memes que celles du pyTree.
// IN: o: la zone python
// IN: xyz: 1: on prend les coords
// IN: loc: 0: en noeuds, 1: en centres, 2: les deux, -1: aucun
// IN: GridCoordinates: nom du container de coordonnees
// IN: FlowSolutionNodes: nom du container de champs en noeuds
// IN: FlowSolutionCenters: nom du container de champs en centres
// OUT: ni,nj,nk: les ni, nj,nk de la zone (structure) ou ni=npoints, 
// nj=nelements en non structure (y compris NGON)
// OUT: cn: connectivite 
// si zone Elements basiques: 1 tableau (comme en pyTree: cnfld * csize)
// si NGON: tableaux NGON [0], NFACE [1], PE [2] (opt)
// OUT: csize: 
// si zone Elements basiques ou NGON: ne (nbre d'elements)
// OUT: cnfld: 
// si zone Elements basiques: nbre de noeuds/element
// si zone NGON: nbre de faces
// OUT: eltType: type d'element (chaine a liberer par l'appelant)
// OUT: fields: un vecteur des pointeurs des champs (partage avec numpy)
// OUT: locs: localisation des champs locs[n] = 0 (en noeuds), = 1 (en centres)
// Retourne 1: zone structuree valide
// Retourne 2: zone non structuree valide
// Retourne 0: zone invalide
//=============================================================================
E_Int K_PYTREE::getFromZone(PyObject* o, E_Int xyz, E_Int loc,
                            char*& varString,
                            vector<E_Float*>& fields,
                            vector<E_Int>& locs,
                            E_Int& ni, E_Int& nj, E_Int& nk,
                            vector<E_Int*>& cn, 
                            E_Int& cnSize, E_Int& cnNfld,
                            char*& eltType,
                            vector<PyArrayObject*>& hook,
                            char* GridCoordinates,
                            char* FlowSolutionNodes,
                            char* FlowSolutionCenters)
{
  E_Int size=0; varString=NULL; eltType=NULL;
  PyObject* t; PyObject* t2;

  // get type
  t = getNodeFromName1(o, "ZoneType");
  if (t == NULL) return 0;
  E_Int s;
  char* type = getValueS(t, s, hook);
  E_Int ret = 2;
  if (K_STRING::cmp(type, s, "Structured") == 0) ret = 1;

  // get dims
  E_Int s0, s1;
  E_Int* d = getValueAI(o, s0, s1, hook);
  ni=0; nj=0; nk=0; cnSize=0; cnNfld=1;
  if (ret == 1) // structure
  {
    if (s0 == 1) { ni = d[0]; nj = 1; nk = 1; }
    else if (s0 == 2) { ni = d[0]; nj = d[1]; nk = 1; }
    else { ni = d[0]; nj = d[1]; nk = d[2]; }
  }
  else // non structure
  {
    ni = d[0]; nj = d[1]; // npoint, nelements
    eltType = new char [20];
    vector<PyObject*> ln;
    getNodesFromType1(o, "Elements_t", ln);
    E_Int* NGON=NULL; E_Int* NFACE=NULL; E_Int* BE=NULL;
    PyObject* NGONt=NULL;
    for (size_t i = 0; i < ln.size(); i++)
    {
      t = ln[i];
      E_Int* dtypep = getValueAI(t, hook);
      E_Int dtype = dtypep[0];
      E_Int dbnd = dtypep[1];
      if (dbnd == 0) // non boundary
      {
        t2 = getNodeFromName1(t, "ElementConnectivity");
        E_Int* v = getValueAI(t2, s0, s1, hook);
        switch (dtype)
        {
          case 2: // node
            BE = v;
            strcpy(eltType, "NODE");
            cnSize = s0; cnNfld = 1;
            break;
            
          case 3: // BAR
            BE = v;
            strcpy(eltType, "BAR"); 
            cnSize = s0; cnNfld = 1;
            break;
            
          case 5: // TRI
            BE = v;
            strcpy(eltType, "TRI"); 
            cnSize = s0/3; cnNfld = 3;
            break;
            
          case 7: // QUAD
            BE = v;
            strcpy(eltType, "QUAD"); 
            cnSize = s0/4; cnNfld = 4;
            break;
            
          case 10: // TETRA
            BE = v;
            strcpy(eltType, "TETRA"); 
            cnSize = s0/4; cnNfld = 4;
            break;
            
          case 12: // PYRA
            BE = v;
            strcpy(eltType, "PYRA"); 
            cnSize = s0/5; cnNfld = 5;
            break;
            
          case 14: // PENTA
            BE = v;
            strcpy(eltType, "PENTA"); 
            cnSize = s0/6; cnNfld = 6;
            break;
            
          case 17: // HEXA
            BE = v;
            strcpy(eltType, "HEXA"); 
            cnSize = s0/8; cnNfld = 8;
            break;
            
          case 20: // MIXED
            break;
            
          case 22: // NGON
            NGON = v;
            NGONt = t;
            cnNfld = s0;
            break;
            
          case 23: // NFACE
            NFACE = v;
            cnNfld = s0;
            break;
            
          default: break;
        }
      }
    }
    // Complete NGON if necessary
    if (NGON != NULL)
    {
      strcpy(eltType, "NGON");
      cn.push_back(NGON);
      cn.push_back(NFACE);
      t2 = getNodeFromName1(NGONt, "ParentElements");
      if (t2 != NULL)
      {
        E_Int* v = getValueAI(t2, s0, s1, hook);
        cn.push_back(v); // PE
      }
    }
    else cn.push_back(BE);
  }

  // get GridCoordinates
  if (xyz == 1)
  { 
    if (GridCoordinates == NULL) t = getNodeFromName1(o, "GridCoordinates");
    else t = getNodeFromName1(o, GridCoordinates);
    if (t != NULL)
    {
      t2 = getNodeFromName1(t, "CoordinateX");
      if (t2 != NULL)
      {
        E_Float* x = getValueAF(t2, hook);
        fields.push_back(x);
        locs.push_back(0);
        add2VarString(varString, size, "CoordinateX");
      }
      t2 = getNodeFromName1(t, "CoordinateY");
      if (t2 != NULL)
      {
        E_Float* y = getValueAF(t2, hook);
        fields.push_back(y);
        locs.push_back(0);
        add2VarString(varString, size, "CoordinateY");
      }
      t2 = getNodeFromName1(t, "CoordinateZ");
      if (t2 != NULL)
      {
        E_Float* z = getValueAF(t2, hook);
        fields.push_back(z);
        locs.push_back(0);
        add2VarString(varString, size, "CoordinateZ");
      }
    }
  }
  
  // Get node fields
  if (loc == 0 || loc == 2)
  {
    if (FlowSolutionNodes == NULL) t = getNodeFromName1(o, "FlowSolution");
    else t = getNodeFromName1(o, FlowSolutionNodes);
    if (t != NULL)
    {
      PyObject* childrens = PyList_GetItem(t, 2);
      E_Int n = PyList_Size(childrens);
      PyObject* l; PyObject* node; char* str;
      for (E_Int i = 0; i < n; i++)
      {
        l = PyList_GetItem(childrens, i);
        node = PyList_GetItem(l, 3);
        if (PyString_Check(node)) str = PyString_AsString(node); // type
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(node)) str = (char*)PyUnicode_AsUTF8(node); 
#endif
        else str = NULL;
        if (K_STRING::cmp(str, "DataArray_t") == 0)
        {
          node = PyList_GetItem(l, 0); // var name
          if (PyString_Check(node)) str = PyString_AsString(node);
#if PY_VERSION_HEX >= 0x03000000
          else if (PyUnicode_Check(node)) str = (char*)PyUnicode_AsUTF8(node); 
#endif
          E_Float* f = getValueAF(l, hook);
          fields.push_back(f);
          locs.push_back(0);
          add2VarString(varString, size, str);
        }
      }
    }
  }
  
  // Get centers field
  if (loc == 1 || loc == 2)
  {
    if (FlowSolutionCenters == NULL)
      t = getNodeFromName1(o, "FlowSolution#Centers");
    else t = getNodeFromName1(o, FlowSolutionCenters);
    if (t != NULL)
    {
      PyObject* childrens = PyList_GetItem(t, 2);
      E_Int n = PyList_Size(childrens);
      PyObject* l; PyObject* node; char* str;
      for (E_Int i = 0; i < n; i++)
      {
        l = PyList_GetItem(childrens, i);
        node = PyList_GetItem(l, 3);
        if (PyString_Check(node)) str = PyString_AsString(node); // type
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(node)) str = (char*)PyUnicode_AsUTF8(node);
#endif
        else str = NULL;
        if (K_STRING::cmp(str, "DataArray_t") == 0)
        {
          node = PyList_GetItem(l, 0); // var name
          if (PyString_Check(node)) str = PyString_AsString(node);
#if PY_VERSION_HEX >= 0x03000000
          else if (PyUnicode_Check(node)) str = (char*)PyUnicode_AsUTF8(node); 
#endif
          E_Float* f = getValueAF(l, hook);
          fields.push_back(f);
          locs.push_back(1);
          add2VarString(varString, size, str);
        }
      }
    }
  }
  
  //printf("varString %s\n", varString);
  return ret;
}
