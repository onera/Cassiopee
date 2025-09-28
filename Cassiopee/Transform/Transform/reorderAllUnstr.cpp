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
# include "stdio.h"
# include "transform.h"
# include "Nuga/include/GeomAlgo.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Redresseur de normales: etant donnee une liste de bloc non structures
   paroi (TRI, QUAD), oriente les blocs de telle sorte que les normales soient
   toutes orientees vers l'exterieur ou l'interieur au choix.
*/
//=============================================================================
PyObject* K_TRANSFORM::reorderAllUnstr(PyObject* self, PyObject* args)
{
  // Load block arrays
  PyObject* listBlks;
  E_Int outward=1; // direction of the normals vers l'exterieur par defaut
  if (!PYPARSETUPLE_(args, O_ I_,
                    &listBlks, &outward))
  {
      return NULL;
  }
  // Check dir
  if (outward != 1 && outward != -1)
  {
    printf("Warning: reorderAllUnstr: direction is invalid. outward direction is assumed.\n");
    outward = 1;
  }

  // Check every array in listBlks
  if (PyList_Check(listBlks) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "reorderAllUnstr: argument must be a list.");
    return NULL;
  }

  E_Int nzone = PyList_Size(listBlks);
  vector<FloatArray*> crds(nzone);
  vector<IntArray*> cnts(nzone); // the one that are retrieved and need to be deleted
  vector<IntArray*>  tri_cnts(nzone); // the one who are going to be processed by reversi
  E_Int posx, posy, posz;
  FloatArray* f; IntArray* cn;
  char* varString; char* eltType;
  vector<E_Bool> is_quad(nzone, false);

  // Extraction des infos pour chaque bloc
  for (E_Int i = 0; i < nzone; i++)
  {
    E_Int nil, njl, nkl;
    PyObject* tpl = PyList_GetItem(listBlks, i);

    E_Int res =
      K_ARRAY::getFromArray(tpl, varString, f, nil, njl, nkl, cn, eltType);

    if (res != 2)
    {
      PyErr_SetString(PyExc_TypeError,
                    "reorderAllUnstr: array is not unstructured.");
      return NULL;
    }

    if ((strcmp(eltType, "TRI") != 0) && (strcmp(eltType, "QUAD") != 0))
    {
      std::cout << "eltType " << eltType << std::endl;
      delete f; delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "reorderAllUnstr: currently only supported for TRI or QUAD arrays.");
      return NULL;
    }

    //check if coordinates exist
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f; delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "reorderAllUnstr: coordinates not found.");
      return NULL;
    }

    if (strcmp(eltType, "QUAD") == 0)
    {
      is_quad[i]=true;
      size_t sz = cn->cols();
      tri_cnts[i] = new IntArray;
      IntArray& tcnt = *tri_cnts[i];
      IntArray& qcnt=*cn;
      tcnt.reserve(3, 2*sz);
      E_Int T[3];
      for (size_t i=0; i < sz; ++i)
      {
        T[0]=qcnt(0,i); T[1]=qcnt(1,i);T[2]=qcnt(2,i);
        tcnt.pushBack(T, T+3);
        T[0]=qcnt(0,i); T[1]=qcnt(2,i);T[2]=qcnt(3,i);
        tcnt.pushBack(T, T+3);
      }
    }
    else tri_cnts[i]=cn;

    //std::cout << "zone type is quad : " << is_quad[i] << std::endl;

    crds[i]=f;
    cnts[i]=cn;
  }//parcours de toutes les zones

  bool otwd = (outward == 1);
  NUGA::GeomAlgo<K_MESH::Triangle>::reversi_chimera_skin(crds, tri_cnts, otwd);

  /*--------------*/
  /* build arrays */
  /*--------------*/
  PyObject* tpl;
  PyObject* l = PyList_New(0);

  for (E_Int i = 0; i < nzone; i++)
  {
    // apply reverse to quads
    if (is_quad[i])
    {
      IntArray& qcnt = *cnts[i];
      for (E_Int j=0; j < qcnt.cols(); ++j)
      {
        if ((*tri_cnts[i])(1,2*j) != qcnt(1,j)) // tri has been reversed
        {
          E_Int* ptr = qcnt.col(j);
          std::reverse(ptr, ptr+4);
        }
      }
    }

    tpl = K_ARRAY::buildArray(*crds[i], varString, *cnts[i], -1, is_quad[i] ? "QUAD" : "TRI");
    delete crds[i]; delete cnts[i];
    if (is_quad[i]) delete tri_cnts[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  return l;
}
