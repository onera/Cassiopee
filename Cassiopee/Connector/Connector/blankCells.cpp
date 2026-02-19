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

# include "connector.h"

using namespace std;
using namespace K_FLD;
extern "C"
{
  void k6searchblankedcellstrix_(
    const E_Int& npts, const E_Float* meshX, const E_Float* meshY,
    const E_Int& nelts, const E_Int* cn1, const E_Int* cn2, const E_Int* cn3,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellstrix2_(
    const E_Int& npts, const E_Float* meshX, const E_Float* meshY,
    const E_Int& nelts, const E_Int* cn1, const E_Int* cn2, const E_Int* cn3,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellstrixd_(
    const E_Int& npts, const E_Float* meshX, const E_Float* meshY,
    const E_Int& nelts, const E_Int* cn1, const E_Int* cn2, const E_Int* cn3,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta, const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellsquadx_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3, const E_Int* cn4,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellsquadx2_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3,const E_Int* cn4,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellsquadxd_(
    const E_Int& npts, const E_Float* meshX, const E_Float* meshY,
    const E_Int& nelts, const E_Int* cn1, const E_Int* cn2, const E_Int* cn3, const E_Int* cn4,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta, const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellstetrax_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3, const E_Int* cn4,
    const E_Float& xmin,  const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellstetrax2_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3, const E_Int* cn4,
    const E_Float& xmin,  const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellstetraxd_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3, const E_Int* cn4,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta, const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellspentax_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3,
    const E_Int* cn4, const E_Int* cn5, const E_Int* cn6,
    const E_Float& xmin,  const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellspentax2_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3,
    const E_Int* cn4, const E_Int* cn5, const E_Int* cn6,
    const E_Float& xmin,  const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellspentaxd_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3,
    const E_Int* cn4, const E_Int* cn5, const E_Int* cn6,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta, const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellshexax_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3, const E_Int* cn4,
    const E_Int* cn5, const E_Int* cn6, const E_Int* cn7, const E_Int* cn8,
    const E_Float& xmin,  const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellshexax2_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3, const E_Int* cn4,
    const E_Int* cn5, const E_Int* cn6, const E_Int* cn7, const E_Int* cn8,
    const E_Float& xmin,  const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellshexaxd_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Int& nelts,
    const E_Int* cn1, const E_Int* cn2, const E_Int* cn3, const E_Int* cn4,
    const E_Int* cn5, const E_Int* cn6, const E_Int* cn7, const E_Int* cn8,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta, const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellsx2d_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY,
    const E_Float& xmin, const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Int* indir, const E_Int& nz,
    const E_Float* z, const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellsx_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz,
    const E_Float* z, const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellsx12d_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked );

  void k6searchblankedcellsx1_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked );

  void k6searchblankednodesx2d_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int& isnot,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankednodesx_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int& isnot, E_Int* cellNatureField, E_Int& isMasked );

  void k6searchblankedcellsx22d_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6searchblankedcellsx2_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    E_Int* cellNatureField, E_Int& isMasked);

  void k6adjustcellnaturefield_(const E_Int& nbBlkCells,
                                const E_Int* blankedCells,
                                E_Int* cellNatFld );
  void k6searchblankedcellsxd_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta, const E_Int& isNot,
    E_Int* cellNatureField, E_Int& isMasked );

  void k6searchblankedcellsxd2d_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta, const E_Int& isNot,
    E_Int* cellNatureField,
    E_Int& isMasked );

  void k6searchblankednodesxd_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Int* diri, const E_Int* dirj,
    const E_Float& delta, const E_Int& isNot,
    E_Int* cellNatureField, E_Int& isMasked );

  void k6searchblankednodesxd2d_(
    const E_Int& npts,
    const E_Float* meshX, const E_Float* meshY,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Int* diri,
    const E_Float& delta, const E_Int& isNot,
    E_Int* cellNatureField,
    E_Int& isMasked );

  void k6searchblankedcellsxd2_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta,
    E_Int* cellNatureField, E_Int& isMasked );

  void k6searchblankedcellsxd1_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY, const E_Float* meshZ,
    const E_Float& xmin, const E_Float& ymin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray, const E_Float& hjray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta,
    E_Int* cellNatureField, E_Int& isMasked );

  void k6searchblankedcellsxd22d_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta,
    E_Int* cellNatureField, E_Int& isMasked );

  void k6searchblankedcellsxd12d_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* meshX, const E_Float* meshY,
    const E_Float& xmin,
    const E_Int& niray, const E_Int& njray,
    const E_Float& hiray,
    const E_Int* indir, const E_Int& nz, const E_Float* z,
    const E_Int* listOfInterpolatedPoints, const E_Int& np,
    const E_Float& delta,
    E_Int* cellNatureField, E_Int& isMasked );
}

# define RELEASEDATA1\
  for (size_t i = 0; i < rest1.size(); i++) \
    RELEASESHAREDB(rest1[i], vectOfObjs1[i], vectOfCoords[i], vectOfConnect1[i]);

# define RELEASEDATA2\
  for (size_t i = 0; i < rest2.size(); i++) \
    RELEASESHAREDB(rest2[i], vectOfObjs2[i], vectOfCellNs[i], vectOfConnect2[i]);

# define RELEASEDATA3\
  for (size_t i = 0; i < rest3.size(); i++) \
    RELEASESHAREDB(rest3[i], vectOfObjs3[i], vectOfBodies[i], vectOfConnect3[i]);

//============================================================================
/* Blank cells defined in arrays by a X-Ray mask
    version in place / getFromArray2 */
//============================================================================
PyObject* K_CONNECTOR::_blankCells(PyObject* self, PyObject* args)
{
  PyObject *coordArrays, *cellNArrays, *bodyArrays;
  E_Float delta; E_Float tol;
  E_Int isNot;
  E_Int dim;
  E_Int blankingType;
  E_Int xraydim1, xraydim2;
  char* cellNName;
  if (!PYPARSETUPLE_(args, OOO_ I_ R_ II_ R_ II_ S_,
                    &coordArrays, &cellNArrays,
                    &bodyArrays, &blankingType, &delta, &dim,
                    &isNot, &tol, &xraydim1, &xraydim2, &cellNName)) return NULL;
  if (PyList_Check(coordArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "_blankCells: first argument must be a list.");
    return NULL;
  }
  if (PyList_Check(cellNArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "_blankCells: second argument must be a list.");
    return NULL;
  }
  if (PyList_Check(bodyArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "_blankCells: third argument must be a list.");
    return NULL;
  }

  if (delta < 0.)
  {
    printf("Warning: _blankCells: delta must be a positive value. Set to default (1.e-10).\n");
    delta = 1.e-10;
  }
  if (isNot != 0 && isNot != 1)
  {
    printf("Warning: _blankCells: masknot is not valid. Set to default (0)\n");
    isNot = 0;
  }
  if (dim != 2 && dim != 3)
  {
    printf("Warning: _blankCells: dim is not valid. Set to default (3)\n");
    dim = 3;
  }
  if (blankingType < -2 || blankingType > 1)
  {
    printf("Warning: _blankCells: blankingType is invalid. Set to default (1)\n");
    blankingType = 1;
  }
  // verification de la coherence des arguments:
  // masknot + cell_intersect_opt est impossible
  if (isNot == 1 && blankingType < 0)
  {
    printf("Warning: _blankCells: cell_intersect_opt criterion and 'not' mask are not compatible.\n");
    printf(" cell_intersect criterion is activated.\n");
    blankingType = 1;
  }
  // Extract infos from coord arrays
  // seulement arrays structures avec coordonnees ici
  E_Int nzonesA = PyList_Size(coordArrays);
  vector<FldArrayF*> vectOfCoords;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<PyObject*> vectOfObjs1;
  vector<FldArrayI*> vectOfConnect1;
  vector<E_Int> rest1;
  for (E_Int i=0; i < nzonesA; i++)
  {
    E_Int nil, njl, nkl;
    FldArrayF* f; FldArrayI* cn;
    char* varString; char* eltType;
    PyObject* array = PyList_GetItem(coordArrays, i);
    E_Int ret = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl,
                                       cn, eltType);
    if (ret != 1)
    {
      RELEASEDATA1;
      PyErr_SetString(PyExc_TypeError,
                      "_blankCells: 1st arg must be a list of structured zones.");
      return NULL;
    }
    vectOfCoords.push_back(f);
    nit.push_back(nil); njt.push_back(njl); nkt.push_back(nkl);
    rest1.push_back(ret);
    vectOfConnect1.push_back(cn);
    vectOfObjs1.push_back(array);

    E_Int posxi = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posyi = K_ARRAY::isCoordinateYPresent(varString);
    E_Int poszi = K_ARRAY::isCoordinateZPresent(varString);
    if (posxi == -1 || posyi == -1 || poszi == -1)
    {
      RELEASEDATA1;
      PyErr_SetString(PyExc_TypeError,
                      "_blankCells: 1st arg must contain coordinates.");
      return NULL;
    }
    posxt.push_back(posxi+1); posyt.push_back(posyi+1); poszt.push_back(poszi+1);
  }

  // Extract infos from celln arrays:
  // si structure: localisation centre ou noeud ok
  // si non structure: localisation aux noeuds uniquement
  E_Int nzonesC = PyList_Size(cellNArrays);
  vector<FldArrayF*> vectOfCellNs;
  vector<E_Int> poscellNt;
  vector<E_Int> nitc; vector<E_Int> njtc; vector<E_Int> nktc;
  vector<PyObject*> vectOfObjs2;
  vector<FldArrayI*> vectOfConnect2;
  vector<E_Int> rest2;

  for (E_Int i = 0; i < nzonesC; i++)
  {
    E_Int nil, njl, nkl;
    FldArrayF* f; FldArrayI* cn;
    char* varString; char* eltType;
    PyObject* array = PyList_GetItem(cellNArrays, i);
    E_Int ret = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl,
                                       cn, eltType);
    if (ret != 1)
    {
      RELEASEDATA1; RELEASEDATA2;
      PyErr_SetString(PyExc_TypeError,
                      "_blankCells: 2nd arg must define a list of structured zones.");
      return NULL;
    }
    vectOfCellNs.push_back(f);
    nitc.push_back(nil); njtc.push_back(njl); nktc.push_back(nkl);
    rest2.push_back(ret);
    vectOfConnect2.push_back(cn);
    vectOfObjs2.push_back(array);

    E_Int posc = K_ARRAY::isNamePresent(cellNName,varString);
    if (posc == -1)
    {
      RELEASEDATA1; RELEASEDATA2;
      PyErr_SetString(PyExc_TypeError,
                      "_blankCells: 2nd arg must contain cellN variable.");
      return NULL;
    }
    poscellNt.push_back(posc+1);
  }

  if ( nzonesC != nzonesA)
  {
    RELEASEDATA1; RELEASEDATA2;
    PyErr_SetString(PyExc_TypeError,
                    "_blankCells: 1st and 2nd args must be of same size. ");
    return NULL;
  }

  // Extract infos from body arrays: non structures
  /* Extraction de la surface de masquage */
  E_Int nzonesB = PyList_Size(bodyArrays);
  E_Int elevationDir = 3;

  vector<PyObject*> vectOfObjs3;
  vector<FldArrayI*> vectOfConnect3;
  vector<E_Int> rest3;
  vector<E_Int> posxb; vector<E_Int> posyb; vector<E_Int> poszb;
  vector<FldArrayF*> vectOfBodies;

  for (E_Int i =0; i < nzonesB; i++)
  {
    E_Int nil, njl, nkl;
    FldArrayF* f; FldArrayI* cn;
    char* varString; char* eltType;
    PyObject* array = PyList_GetItem(bodyArrays, i);
    E_Int ret = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl,
                                       cn, eltType);
    if (ret != 2)
    {
      RELEASEDATA1; RELEASEDATA2; RELEASEDATA3;
      PyErr_SetString(PyExc_TypeError,
                      "_blankCells: 3rd arg must be a list of unstructured zones.");
      return NULL;
    }
    vectOfBodies.push_back(f);
    vectOfConnect3.push_back(cn);
    vectOfObjs3.push_back(array);
    rest3.push_back(ret);

    if (strcmp(eltType,"BAR") != 0 && strcmp(eltType,"TRI") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "_blankCells: body arrays must be all of TRI or BAR type.");
      RELEASEDATA1; RELEASEDATA2; RELEASEDATA3;
      return NULL;
    }

    E_Int posxi = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posyi = K_ARRAY::isCoordinateYPresent(varString);
    E_Int poszi = K_ARRAY::isCoordinateZPresent(varString);
    if (posxi == -1 || posyi == -1 || poszi == -1)
    {
      RELEASEDATA1; RELEASEDATA2; RELEASEDATA3;
      PyErr_SetString(PyExc_TypeError,
                      "_blankCells: 3rd arg must contain coordinates.");
      return NULL;
    }

    posxb.push_back(posxi+1);
    posyb.push_back(posyi+1);
    poszb.push_back(poszi+1);

    if (strcmp(eltType, "BAR") == 0 || dim == 2) elevationDir = 2;
  }

  // verification de la coherence des dimensions
  for (E_Int zone = 0; zone < nzonesA; zone++) // structured
  {
    E_Int ni = nit[zone]; E_Int nic = nitc[zone];
    E_Int nj = njt[zone]; E_Int njc = njtc[zone];
    E_Int nk = nkt[zone]; E_Int nkc = nktc[zone];
    if (blankingType != 0)
    {
      ni = K_FUNC::E_max(ni-1,1);
      nj = K_FUNC::E_max(nj-1,1);
      nk = K_FUNC::E_max(nk-1,1);
    }
    if (ni != nic || nj != njc || nk != nkc)
    {
      PyErr_SetString(PyExc_TypeError,
                      "_blankCells: dimensions of coord and celln arrays do not correspond.");
      RELEASEDATA1; RELEASEDATA2; RELEASEDATA3;
      return NULL;
    }
  }

  // Blanking...
  E_Int dim1 = E_Int(xraydim1); E_Int dim2 = E_Int(xraydim2);
  vector<FldArrayI*> vectOfCellNI;
  for (E_Int is = 0; is < nzonesA; is++)
  {
    E_Float* cellNz = vectOfCellNs[is]->begin();
    E_Int ncells = vectOfCellNs[is]->getSize();
    FldArrayI* cellnI = new FldArrayI(ncells);
    E_Int* cellnp = cellnI->begin();

    #pragma omp parallel for
    for (E_Int i = 0; i < ncells; i++) cellnp[i] = E_Int(cellNz[i]);
    vectOfCellNI.push_back(cellnI);
  }

  blankCellsStruct(elevationDir, isNot, blankingType, delta, tol, dim1, dim2,
                   posxt, posyt, poszt, nit, njt, nkt, vectOfCoords,
                   vectOfCellNI, posxb, posyb, poszb, vectOfBodies, vectOfConnect3);

  for (E_Int noc = 0; noc < nzonesA; noc++)
  {
    E_Int* cellNI = vectOfCellNI[noc]->begin();
    E_Int posc = poscellNt[noc];
    E_Int ncells = vectOfCellNI[noc]->getSize();
    E_Float* cellNp = vectOfCellNs[noc]->begin(posc);

    #pragma omp parallel for
    for (E_Int ind = 0; ind < ncells; ind++)
      cellNp[ind] = E_Float(cellNI[ind]);
    delete vectOfCellNI[noc];
  }

  RELEASEDATA1; RELEASEDATA2; RELEASEDATA3;
  Py_INCREF(Py_None);
  return Py_None;
}

//============================================================================
/* Blank cells defined in arrays by a X-Ray mask */
//============================================================================
PyObject* K_CONNECTOR::blankCells(PyObject* self, PyObject* args)
{
  PyObject* coordArrays; PyObject* cellnArrays;
  PyObject* bodyArrays;
  char* cellNName;
  E_Float delta; E_Float tol;
  E_Int isNot;
  E_Int dim;
  E_Int blankingType;
  E_Int xraydim1, xraydim2;

  if (!PYPARSETUPLE_(args, OOO_ I_ R_ II_ R_ II_ S_,
                    &coordArrays, &cellnArrays,
                    &bodyArrays, &blankingType, &delta, &dim,
                    &isNot, &tol, &xraydim1, &xraydim2, &cellNName))
  {
      return NULL;
  }

  if (PyList_Check(coordArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankCells: first argument must be a list.");
    return NULL;
  }
  if (PyList_Check(cellnArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankCells: second argument must be a list.");
    return NULL;
  }
  if (PyList_Check(bodyArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankCells: third argument must be a list.");
    return NULL;
  }

  if (delta < 0.)
  {
    printf("Warning: blankCells: delta must be a positive value. Set to default (1.e-10).\n");
    delta = 1.e-10;
  }
  if (isNot != 0 && isNot != 1)
  {
    printf("Warning: blankCells: masknot is not valid. Set to default (0)\n");
    isNot = 0;
  }
  if (dim != 2 && dim != 3)
  {
    printf("Warning: blankCells: dim is not valid. Set to default (3)\n");
    dim = 3;
  }
  if (blankingType < -2 || blankingType > 1)
  {
    printf("Warning: blankCells: blankingType is invalid. Set to default (1)\n");
    blankingType = 1;
  }
  // verification de la coherence des arguments:
  // masknot + cell_intersect_opt est impossible
  if (isNot == 1 && blankingType < 0)
  {
    printf("Warning: blankCells: cell_intersect_opt criterion and 'not' mask are not compatible.\n");
    printf(" cell_intersect criterion is activated.\n");
    blankingType = 1;
  }
  // Extract infos from coord arrays
  // seulement arrays structures avec coordonnees ici
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    coordArrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int ns = structF.size(); E_Int nu = unstrF.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankCells: 1st list of arrays is not valid.");
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objs[is], structF[is]);
    for (E_Int i = 0; i < nu; i++)
      RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
    return NULL;
  }

  // verification des coordonnees des domaines a masquer
  E_Int posxi, posyi, poszi;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu;
  vector<char*> varStrings;

  for (E_Int i = 0; i < ns; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarString[i]);
    posyi = K_ARRAY::isCoordinateYPresent(structVarString[i]);
    poszi = K_ARRAY::isCoordinateZPresent(structVarString[i]);
    if (posxi == -1 || posyi == -1 || poszi == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankCells: arrays must contain coordinates.");

      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objs[is], structF[is]);
      for (E_Int is = 0; is < nu; is++)
        RELEASESHAREDU(obju[is], unstrF[is], cnt[is]);
    }
    posxi++; posyi++; poszi++;
    posxs.push_back(posxi); posys.push_back(posyi); poszs.push_back(poszi);
  }
  for (E_Int i = 0; i < nu; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarString[i]);
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarString[i]);
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarString[i]);
    if (posxi == -1 || posyi == -1 || poszi == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankCells: arrays must contain coordinates.");

      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objs[is], structF[is]);
      for (E_Int is = 0; is < nu; is++)
        RELEASESHAREDU(obju[is], unstrF[is], cnt[is]);
    }
    posxi++; posyi++; poszi++;
    posxu.push_back(posxi); posyu.push_back(posyi); poszu.push_back(poszi);
  }

  // Extract infos from celln arrays:
  // si structure: localisation centre ou noeud ok
  // si non structure: localisation aux noeuds uniquement
  vector<E_Int> resc;
  vector<char*> structVarStringc; vector<char*> unstrVarStringc;
  vector<FldArrayF*> structFc; vector<FldArrayF*> unstrFc;
  vector<E_Int> nitc; vector<E_Int> njtc; vector<E_Int> nktc;
  vector<FldArrayI*> cntc; vector<char*> eltTypec;
  vector<PyObject*> objsc, objuc;
  skipNoCoord = false;
  skipStructured = false;
  skipUnstructured = false;
  skipDiffVars = true;
  E_Int isOkc = K_ARRAY::getFromArrays(
    cellnArrays, resc, structVarStringc, unstrVarStringc,
    structFc, unstrFc, nitc, njtc, nktc, cntc, eltTypec, objsc, objuc,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nuc = unstrFc.size(); E_Int nsc = structFc.size();
  if (nsc != ns || nuc != nu)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankCells: coordinates and cellN fields must be of same size.");
    for (E_Int i = 0; i < ns; i++)
      RELEASESHAREDS(objs[i], structF[i]);
    for(E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsc[is], structFc[is]);
    for (E_Int i = 0; i < nu; i++)
      RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
    for(E_Int is = 0; is < nuc; is++)
      RELEASESHAREDU(objuc[is], unstrFc[is], cntc[is]);
    return NULL;
  }
  if (isOkc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankCells: 2nd list of arrays is not valid.");
    for (E_Int i = 0; i < ns; i++)
    {
      RELEASESHAREDS(objs[i], structF[i]);
      RELEASESHAREDS(objsc[i], structFc[i]);
    }
    for (E_Int i = 0; i < nu; i++)
    {
      RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
      RELEASESHAREDU(objuc[i], unstrFc[i], cntc[i]);
    }
    return NULL;
  }

  // recherche de la variable celln
  char* varStringc = NULL;
  if (ns > 0) varStringc = structVarStringc[0];
  else if (nu > 0) varStringc = unstrVarStringc[0];

  E_Int posc = K_ARRAY::isNamePresent(cellNName,varStringc);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankCells: celln variable not found.");
    for (E_Int i = 0; i < ns; i++)
    {
      RELEASESHAREDS(objs[i], structF[i]);
      RELEASESHAREDS(objsc[i], structFc[i]);
    }
    for (E_Int i = 0; i < nu; i++)
    {
      RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
      RELEASESHAREDU(objuc[i], unstrFc[i], cntc[i]);
    }
    return NULL;
  }
  for (E_Int i = 0; i < ns; i++)
  {
    E_Int poscc = K_ARRAY::isNamePresent(cellNName,structVarStringc[i]);
    if (poscc != posc)
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankCells: celln variable not found for one structured array.");
      for (E_Int i = 0; i < ns; i++)
      {
        RELEASESHAREDS(objs[i], structF[i]);
        RELEASESHAREDS(objsc[i], structFc[i]);
      }
      for (E_Int i = 0; i < nu; i++)
      {
        RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
        RELEASESHAREDU(objuc[i], unstrFc[i], cntc[i]);
      }
      return NULL;
    }
  }
  for (E_Int i = 0; i < nu; i++)
  {
    E_Int poscc = K_ARRAY::isNamePresent(cellNName,unstrVarStringc[i]);
    if (poscc != posc)
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankCells: celln variable not found for one unstructured array.");
      for (E_Int i = 0; i < ns; i++)
      {
        RELEASESHAREDS(objs[i], structF[i]);
        RELEASESHAREDS(objsc[i], structFc[i]);
      }
      for (E_Int i = 0; i < nu; i++)
      {
        RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
        RELEASESHAREDU(objuc[i], unstrFc[i], cntc[i]);
      }
      return NULL;
    }
  }
  posc++;

  // Extract infos from body arrays: non structures
  /* Extraction de la surface de masquage */
  vector<E_Int> reslb;
  vector<char*> structVarStringb; vector<char*> unstrVarStringb;
  vector<FldArrayF*> structbF; vector<FldArrayF*> unstrbF;
  vector<E_Int> nib; vector<E_Int> njb; vector<E_Int> nkb;
  vector<FldArrayI*> cnb; vector<char*> eltTypeb;
  vector<PyObject*> objsb, objub;
  skipNoCoord = true;
  skipStructured = true;
  skipUnstructured = false;
  skipDiffVars = true;
  E_Int resb = K_ARRAY::getFromArrays(
    bodyArrays, reslb, structVarStringb, unstrVarStringb,
    structbF, unstrbF, nib, njb, nkb, cnb, eltTypeb, objsb, objub,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  if (resb != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankCells: body arrays must be unstructured.");
    for (E_Int i = 0; i < ns; i++)
    {
      RELEASESHAREDS(objs[i], structF[i]);
      RELEASESHAREDS(objsc[i], structFc[i]);
    }
    for (E_Int i = 0; i < nu; i++)
    {
      RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
      RELEASESHAREDU(objuc[i], unstrFc[i], cntc[i]);
    }
    E_Int nub = structbF.size();
    for (E_Int ib = 0; ib < nub; ib++)
      RELEASESHAREDU(objub[ib], unstrbF[ib], cnb[ib]);
    return NULL;
  }

  E_Int nzonesb = unstrbF.size();
  char* eltTypeb0 = NULL;

  // verification du type d'elements
  for (E_Int i = 0; i < nzonesb; i++)
  {
    if (eltTypeb0 == NULL)
    {
      if (strcmp(eltTypeb[i], "TRI") != 0 && strcmp(eltTypeb[i], "BAR") != 0)
      {
        PyErr_SetString(PyExc_TypeError,
                        "blankCells: body arrays must be all of TRI or BAR type.");
        for (E_Int i = 0; i < ns; i++)
        {
          RELEASESHAREDS(objs[i], structF[i]);
          RELEASESHAREDS(objsc[i], structFc[i]);
        }
        for (E_Int i = 0; i < nu; i++)
        {
          RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
          RELEASESHAREDU(objuc[i], unstrFc[i], cntc[i]);
        }
        for (E_Int ib = 0; ib < nzonesb; ib++)
          RELEASESHAREDU(objub[ib], unstrbF[ib], cnb[ib]);
        return NULL;
      }
      else eltTypeb0 = eltTypeb[i];
    }
    else
    {
      if (strcmp(eltTypeb0,eltTypeb[i] ) != 0)
      {
        PyErr_SetString(PyExc_TypeError,
                        "blankCells: body arrays must be all of TRI or BAR type.");
        for (E_Int i = 0; i < ns; i++)
        {
          RELEASESHAREDS(objs[i], structF[i]);
          RELEASESHAREDS(objsc[i], structFc[i]);
        }
        for (E_Int i = 0; i < nu; i++)
        {
          RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
          RELEASESHAREDU(objuc[i], unstrFc[i], cntc[i]);
        }
        for (E_Int ib = 0; ib < nzonesb; ib++)
          RELEASESHAREDU(objub[ib], unstrbF[ib], cnb[ib]);
        return NULL;
      }
    }
  }

  // verification des coordonnees
  vector<E_Int> posxb; vector<E_Int> posyb; vector<E_Int> poszb;
  E_Int nbodies = unstrbF.size();
  for (E_Int i = 0; i < nbodies; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarStringb[i]);
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarStringb[i]);
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarStringb[i]);
    if (posxi == -1 || posyi == -1 || poszi == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankCells: body arrays must contain coordinates.");
      for (E_Int i = 0; i < ns; i++)
      {
        RELEASESHAREDS(objs[i], structF[i]);
        RELEASESHAREDS(objsc[i], structFc[i]);
      }
      for (E_Int i = 0; i < nu; i++)
      {
        RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
        RELEASESHAREDU(objuc[i], unstrFc[i], cntc[i]);
      }
      for (E_Int ib = 0; ib < nzonesb; ib++)
        RELEASESHAREDU(objub[ib], unstrbF[ib], cnb[ib]);
      return NULL;
    }
    posxi++; posyi++; poszi++;
    posxb.push_back(posxi); posyb.push_back(posyi); poszb.push_back(poszi);
  }

  E_Int elevationDir = 3;
  if (strcmp(eltTypeb0, "BAR") == 0 || dim == 2) elevationDir = 2;
  // verification de la coherence des dimensions
  for (E_Int zone = 0; zone < ns; zone++) // structured
  {
    E_Int ni = nit[zone]; E_Int nic = nitc[zone];
    E_Int nj = njt[zone]; E_Int njc = njtc[zone];
    E_Int nk = nkt[zone]; E_Int nkc = nktc[zone];
    if (blankingType != 0)
    {
      ni = K_FUNC::E_max(ni-1,1);
      nj = K_FUNC::E_max(nj-1,1);
      nk = K_FUNC::E_max(nk-1,1);
    }
    if (ni != nic || nj != njc || nk != nkc)
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankCells: dimensions of coord and celln arrays do not correspond.");
      for (E_Int is = 0; is < ns; is++)
      {
        RELEASESHAREDS(objs[is], structF[is]);
        RELEASESHAREDS(objsc[is], structFc[is]);
      }
      for (E_Int iu = 0; iu < nu; iu++)
      {
        RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
        RELEASESHAREDU(objuc[iu], unstrFc[iu], cntc[iu]);
      }
      for (E_Int ib = 0; ib < nzonesb; ib++)
        RELEASESHAREDU(objub[ib], unstrbF[ib], cnb[ib]);
      return NULL;
    }
  }
  for (E_Int zone = 0; zone < nu; zone++) // unstructured
  {
    if (blankingType == 0) // masquage en noeuds
    {
      if (unstrFc[zone]->getSize() != unstrF[zone]->getSize())
      {
        PyErr_SetString(PyExc_TypeError,
                        "blankCells: dimensions of coord and celln arrays do not correspond.");
        for (E_Int is = 0; is < ns; is++)
        {
          RELEASESHAREDS(objs[is], structF[is]);
          RELEASESHAREDS(objsc[is], structFc[is]);
        }
        for (E_Int iu = 0; iu < nu; iu++)
        {
          RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
          RELEASESHAREDU(objuc[iu], unstrFc[iu], cntc[iu]);
        }
        for (E_Int ib = 0; ib < nzonesb; ib++)
          RELEASESHAREDU(objub[ib], unstrbF[ib], cnb[ib]);
        return NULL;
      }
    }
    else // Masquage en centres
    {
        E_Int iserr=0;
        if  (strcmp(eltType[zone],"NGON")==0)
        {
            iserr=1;
            PyErr_SetString(PyExc_TypeError,
                            "blankCells: blanking with cell intersect option is not implemented for NGONs.");
        }
        else
        {
            if (unstrFc[zone]->getSize() != cnt[zone]->getSize())
            {
                iserr=1;
                PyErr_SetString(PyExc_TypeError,
                                "blankCells: dimensions of coord and celln arrays do not correspond.");
            }
        }
        if (iserr==1)
        {
            for (E_Int is = 0; is < ns; is++)
            {
              RELEASESHAREDS(objs[is], structF[is]);
              RELEASESHAREDS(objsc[is], structFc[is]);
            }
            for (E_Int iu = 0; iu < nu; iu++)
            {
              RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
              RELEASESHAREDU(objuc[iu], unstrFc[iu], cntc[iu]);
            }
            for (E_Int ib = 0; ib < nzonesb; ib++)
              RELEASESHAREDU(objub[ib], unstrbF[ib], cnb[ib]);
            return NULL;
        }
    }
  } // fin test coherence dimensions
  // Build arrays
  PyObject* l = PyList_New(0);
  PyObject* tpl;

  // Blanking...
  E_Int dim1 = E_Int(xraydim1); E_Int dim2 = E_Int(xraydim2);
  if (ns > 0) // structure
  {
    vector<FldArrayI*> cellns; // recup du celln structure
    for (E_Int is = 0; is < ns; is++)
    {
      E_Int ncells = structFc[is]->getSize();
      FldArrayI* celln = new FldArrayI(ncells);
      E_Float* fp = structFc[is]->begin(posc);
      E_Int* cellnp = celln->begin();
      #pragma omp parallel for
      for (E_Int i = 0; i < ncells; i++) cellnp[i] = E_Int(fp[i]);
      cellns.push_back(celln);
    }
    blankCellsStruct(elevationDir, isNot, blankingType, delta, tol, dim1, dim2,
                     posxs, posys, poszs, nit, njt, nkt, structF,
                     cellns, posxb, posyb, poszb, unstrbF, cnb);
    for (E_Int is = 0; is < ns; is++)
    {
      E_Int api = 1; // TODO
      E_Int ncells = cellns[is]->getSize();
      E_Int* fp = cellns[is]->begin();
      FldArrayF* cellnout = new FldArrayF(ncells);
      E_Float* cellnp = cellnout->begin();
      #pragma omp parallel for
      for (E_Int i = 0; i < ncells; i++) cellnp[i] = E_Float(fp[i]);
      tpl = K_ARRAY::buildArray3(*cellnout, cellNName,
                                 nitc[is], njtc[is], nktc[is], api);
      delete cellns[is];
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
      delete cellnout;
      RELEASESHAREDS(objs[is], structF[is]);
      RELEASESHAREDS(objsc[is], structFc[is]);
    }
  }

  if (nu > 0)
  {
    vector<FldArrayI*> cellnu;// recup du celln non structure
    for (E_Int iu = 0; iu < nu; iu++)
    {
      E_Int ncells = unstrFc[iu]->getSize();
      FldArrayI* celln = new FldArrayI(ncells);
      E_Float* fp = unstrFc[iu]->begin(posc);
      E_Int* cellnp = celln->begin();
      #pragma omp parallel for
      for (E_Int i = 0; i < ncells; i++) cellnp[i] = E_Int(fp[i]);
      cellnu.push_back(celln);
    }
    blankCellsUnstr(elevationDir, isNot, blankingType, delta, tol, dim1, dim2,
                    posxu, posyu, poszu, unstrF, cnt,
                    cellnu, cntc, posxb, posyb, poszb, unstrbF, cnb);
    for (E_Int iu = 0; iu < nu; iu++)
    {
      E_Int ncells = cellnu[iu]->getSize();
      E_Int* fp = cellnu[iu]->begin();
      FldArrayF* cellnout = new FldArrayF(ncells);
      E_Float* cellnp = cellnout->begin();
      E_Int api = cellnout->getApi();
      #pragma omp parallel for
      for (E_Int i = 0; i < ncells; i++) cellnp[i] = E_Float(fp[i]);
      FldArrayI* cnout = new K_FLD::FldArrayI(*cntc[iu]);
      tpl = K_ARRAY::buildArray3(*cellnout, cellNName, *cnout, eltTypec[0], api);
      delete cellnu[iu];
      delete cnout;
      PyList_Append(l, tpl); Py_DECREF(tpl);
      delete cellnout;
      RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
      RELEASESHAREDU(objuc[iu], unstrFc[iu], cntc[iu]);
    }
  }

  for (E_Int ib = 0; ib < nzonesb; ib++)
    RELEASESHAREDU(objub[ib], unstrbF[ib], cnb[ib]);
  return l;
}
//============================================================================
/* blankCells (c function) for unstructured arrays: NODE type is only valid */
//============================================================================
void K_CONNECTOR::blankCellsUnstr(
  E_Int elevationDir, E_Int isNot,  E_Int blankingType,
  E_Float delta, E_Float tol, E_Int dim1, E_Int dim2,
  vector<E_Int>& posxt, vector<E_Int>& posyt,
  vector<E_Int>& poszt,
  vector<FldArrayF*>& blankedCoords,  vector<FldArrayI*>& cnt,
  vector<FldArrayI*>& cellns, vector<FldArrayI*>& cntc,
  vector<E_Int>& posxb, vector<E_Int>& posyb,
  vector<E_Int>& poszb,
  vector<FldArrayF*>& fieldsb,
  vector<FldArrayI*>& cnb)
{
  list<XRayPlane*> planes;
  // creation du masque
  E_Float xmin, ymin, zmin, xmax, ymax, zmax;
  E_Float xminz, yminz, zminz, xmaxz, ymaxz, zmaxz;

  compCharacteristics(isNot, elevationDir, dim1, dim2, tol, delta,
                      posxb, posyb, poszb, fieldsb, cnb, planes,
                      xmin, ymin, zmin, xmax, ymax, zmax);
  xmin = xmin - delta; xmax = xmax + delta;
  ymin = ymin - delta; ymax = ymax + delta;
  zmin = zmin - delta; zmax = zmax + delta;

  // masquage des domaines a masquer
  E_Int nzones = blankedCoords.size();
  E_Int npts;
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    // intersection des bbox ?
    npts = blankedCoords[zone]->getSize();
    K_COMPGEOM::boundingBoxUnstruct(npts,
                            blankedCoords[zone]->begin(posxt[zone]),
                            blankedCoords[zone]->begin(posyt[zone]),
                            blankedCoords[zone]->begin(poszt[zone]),
                            xminz, yminz, zminz, xmaxz, ymaxz, zmaxz);
    E_Int intersect =
      K_COMPGEOM::compBoundingBoxIntersection(
        xmin, xmax, ymin, ymax, zmin, zmax,
        xminz, xmaxz, yminz, ymaxz, zminz, zmaxz, 1.e-6);
    if (intersect != 0)
    {
      FldArrayI& cellN0 = *cellns[zone]; // 0: M / 2: I / 1: N
      E_Int ncells = cellN0.getSize();
      FldArrayI cellN(ncells); cellN.setAllValuesAt(1);
      E_Int* cellN0p = cellN0.begin();
      E_Int* cellNp = cellN.begin();
      #pragma omp parallel for
      for (E_Int ind = 0; ind < ncells; ind++)
      {
        if (cellN0p[ind] == 0) cellNp[ind] = -1;
        else if (cellN0p[ind] == 2) cellNp[ind] = 0;
      }
      FldArrayI blankedCells(ncells);// -1: masque / 0: interpole / 1: normal
      blankedCells.setAllValuesAt(1);
      E_Int isMasked = searchForBlankedCellsUnstr(
        elevationDir, blankingType, isNot, delta, planes,
        posxt[zone], posyt[zone], poszt[zone],
        *blankedCoords[zone], *cnt[zone],
        blankedCells);

      if (isMasked == 1)
      {
        k6adjustcellnaturefield_(blankedCells.getSize(),
                                 blankedCells.begin(),
                                 cellN.begin());
        E_Int* cellN0p = cellN0.begin(); // 1/0/2
        E_Int* cellNp = cellN.begin();
        #pragma omp parallel for
        for (E_Int ind = 0; ind < cellN.getSize(); ind++)
        {
          if (cellNp[ind] == -1) cellN0p[ind] = 0;
          else if (cellNp[ind] == 0) cellN0p[ind] = 2;
        }
      }
    }
  }
  // nettoyage
  for (list<XRayPlane*>::iterator itr = planes.begin();
       itr != planes.end(); itr++)
  {delete [] (*itr)->tempZ; delete *itr; }
  return;
}
//============================================================================
/* blankCells (c function) */
//============================================================================
void K_CONNECTOR::blankCellsStruct(
  E_Int elevationDir, E_Int isNot, E_Int blankingType,
  E_Float delta, E_Float tol, E_Int dim1, E_Int dim2,
  vector<E_Int>& posxt, vector<E_Int>& posyt,
  vector<E_Int>& poszt,
  vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt,
  vector<FldArrayF*>& blankedCoords,
  vector<FldArrayI*>& cellns,
  vector<E_Int>& posxb, vector<E_Int>& posyb,
  vector<E_Int>& poszb,
  vector<FldArrayF*>& fieldsb,
  vector<FldArrayI*>& cnb)
{
  list<XRayPlane*> planes;
  // creation du masque
  E_Float xmin, ymin, zmin, xmax, ymax, zmax;
  E_Float xminz, yminz, zminz, xmaxz, ymaxz, zmaxz;

  compCharacteristics(isNot, elevationDir, dim1, dim2, tol, delta,
                      posxb, posyb, poszb, fieldsb, cnb, planes,
                      xmin, ymin, zmin, xmax, ymax, zmax);

  xmin = xmin - delta; xmax = xmax + delta;
  ymin = ymin - delta; ymax = ymax + delta;
  zmin = zmin - delta; zmax = zmax + delta;

  // masquage des domaines a masquer
  E_Int nzones = blankedCoords.size();
  E_Int intersect = 1;
  E_Int npts;
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    intersect = 1;
    if (isNot == 0) // test intersection des bbox que pour les cas de masques classiques
    {
      npts = blankedCoords[zone]->getSize();
      // intersection des bbox ?
      K_COMPGEOM::boundingBoxUnstruct(npts,
                              blankedCoords[zone]->begin(posxt[zone]),
                              blankedCoords[zone]->begin(posyt[zone]),
                              blankedCoords[zone]->begin(poszt[zone]),
                              xminz, yminz, zminz, xmaxz, ymaxz, zmaxz);
      
      intersect = K_COMPGEOM::compBoundingBoxIntersection(
        xmin, xmax, ymin, ymax, zmin, zmax,
        xminz, xmaxz, yminz, ymaxz, zminz, zmaxz, 1.e-6);
    }
    if (intersect != 0)
    {
      FldArrayI& cellN0 = *cellns[zone]; // 0: M / 2: I / 1: N
      E_Int ncells = cellN0.getSize();
      FldArrayI cellN(ncells); cellN.setAllValuesAt(1);

      #pragma omp parallel for
      for (E_Int ind = 0; ind < ncells; ind++)
      {
        if (cellN0[ind] == 0) cellN[ind] = -1;
        else if (cellN0[ind] == 2) cellN[ind] = 0;
      }
      FldArrayI blankedCells(ncells);// -1: masque / 0: interpole / 1: normal
      blankedCells.setAllValuesAt(1);

      E_Int isMasked =
      searchForBlankedCellsStruct(elevationDir, blankingType,
        isNot, delta, planes,
        nit[zone], njt[zone], nkt[zone],
        posxt[zone], posyt[zone], poszt[zone],
        *blankedCoords[zone],
        blankedCells);

      if (isMasked == 1)
      {
        k6adjustcellnaturefield_(blankedCells.getSize(),
                                    blankedCells.begin(),
                                    cellN.begin());

        E_Int* cellN0p = cellN0.begin(); // 1/0/2
        E_Int* cellNp = cellN.begin();

        #pragma omp parallel for
        for (E_Int ind = 0; ind < cellN.getSize(); ind++)
        {
          if (cellNp[ind] == -1) cellN0p[ind] = 0;
          else if (cellNp[ind] == 0) cellN0p[ind] = 2;
        }
      }
    }
  }
  // nettoyage
  for (list<XRayPlane*>::iterator itr = planes.begin();
       itr != planes.end(); itr++)
  {delete [] (*itr)->tempZ; delete *itr; }
}

//=============================================================================
/* Recherche des cellules masquees pour 1 domaine donne
   retourne 1: masque
            0: pas masque
           -1: echec erreur de critere
*/
//=============================================================================
E_Int K_CONNECTOR::searchForBlankedCellsStruct(
  E_Int elevationDir, E_Int blankingType,
  E_Int isNot, E_Float delta,
  list<XRayPlane*>& planes,
  E_Int ni, E_Int nj, E_Int nk,
  E_Int posx, E_Int posy, E_Int posz,
  FldArrayF& field,
  FldArrayI& blankedCells)
{
  E_Int masked = 0;

  list<XRayPlane*>::iterator itr;
  for (itr = planes.begin(); itr != planes.end(); itr++)
  {
    XRayPlane* p = *itr;
    E_Int maskedl = 0;
    switch (blankingType)
    {
      case 1: // cell_intersect
        if (elevationDir == 2)
        {
          k6searchblankedcellsx2d_(
            ni, nj, nk, field.begin(posx), field.begin(posy),
            p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
            p->Z.getSize(), p->Z.begin(), isNot,
            blankedCells.begin(), maskedl);
        }
        else
        {
          k6searchblankedcellsx_(ni, nj, nk, field.begin(posx),
                                  field.begin(posy),field.begin(posz),
                                  p->xmin, p->ymin, p->ni, p->nj,
                                  p->hi, p->hj, p->indir.begin(),
                                  p->Z.getSize(), p->Z.begin(), isNot,
                                  blankedCells.begin(), maskedl);

        }
        break;

      case 0: // node_in
        if (elevationDir == 2)
        {
          k6searchblankednodesx2d_(
            field.getSize(), field.begin(posx), field.begin(posy),
            p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
            p->Z.getSize(), p->Z.begin(), isNot,
            blankedCells.begin(), maskedl);
        }
        else // dim 3
        {
          k6searchblankednodesx_(field.getSize(), field.begin(posx),
                                 field.begin(posy),field.begin(posz),
                                 p->xmin, p->ymin, p->ni, p->nj,
                                 p->hi, p->hj, p->indir.begin(),
                                 p->Z.getSize(), p->Z.begin(), isNot,
                                 blankedCells.begin(), maskedl);
        }
        break;
      case -1: // cell_intersect_opt + depth=1
        if (elevationDir == 2)
        {
          k6searchblankedcellsx12d_(
            ni, nj, nk, field.begin(posx), field.begin(posy),
            p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
            p->Z.getSize(), p->Z.begin(),
            blankedCells.begin(), maskedl);
        }
        else
        {
          k6searchblankedcellsx1_(ni, nj, nk, field.begin(posx),
                                  field.begin(posy),field.begin(posz),
                                  p->xmin, p->ymin, p->ni, p->nj,
                                  p->hi, p->hj, p->indir.begin(),
                                  p->Z.getSize(), p->Z.begin(),
                                  blankedCells.begin(), maskedl);
        }
        break;

      case -2: // cell_intersect_opt + depth=2
        if (elevationDir == 2)
        {
          k6searchblankedcellsx22d_(
            ni, nj, nk, field.begin(posx), field.begin(posy),
            p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
            p->Z.getSize(), p->Z.begin(),
            blankedCells.begin(), maskedl);
        }
        else
        {
          k6searchblankedcellsx2_(ni, nj, nk, field.begin(posx),
                                  field.begin(posy),field.begin(posz),
                                  p->xmin, p->ymin, p->ni, p->nj,
                                  p->hi, p->hj, p->indir.begin(),
                                  p->Z.getSize(), p->Z.begin(),
                                  blankedCells.begin(), maskedl);
        }
        break;

      default:
        printf("Warning: searchForBlankedCellsStruct: blankingType value is invalid.\n");
        return -1;
    }
    masked = K_FUNC::E_max(maskedl, masked);
  }
  if (delta > 0.)
  {
    return holeExpansionStruct(elevationDir, blankingType, isNot,
                               delta, planes, ni, nj, nk, posx, posy, posz,
                               field, blankedCells);

  }
  return (masked > 0);
}

//=============================================================================
/* Recherche des cellules masquees pour 1 domaine donne
   retourne 1: masque
            0: pas masque
           -1: echec erreur de critere
*/
//=============================================================================
E_Int K_CONNECTOR::searchForBlankedCellsUnstr(
  E_Int elevationDir, E_Int blankingType,
  E_Int isNot, E_Float delta,
  list<XRayPlane*>& planes,
  E_Int posx, E_Int posy, E_Int posz,
  FldArrayF& field, FldArrayI& cn,
  FldArrayI& blankedCells)
{
  E_Int nelts = cn.getNfld();
  E_Int masked = 0;
  list<XRayPlane*>::iterator itr;
  for (itr = planes.begin(); itr != planes.end(); itr++)
  {
    XRayPlane* p = *itr;
    E_Int maskedl = 0;
    switch (blankingType)
    {
      case 1: //cell_intersect
        if (elevationDir == 2)
        {
          if (nelts == 3) // TRI
            k6searchblankedcellstrix_(
              field.getSize(), field.begin(posx), field.begin(posy),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),
              p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(), isNot,
              blankedCells.begin(), maskedl);

          else if (nelts == 4) // QUAD
            k6searchblankedcellsquadx_(
              field.getSize(), field.begin(posx), field.begin(posy),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),cn.begin(4),
              p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(), isNot,
              blankedCells.begin(), maskedl);
          else {printf("Warning: searchForBlankedCellsUnstr: element type not valid.\n"); return -1;}
        }
        else
        {
          if (nelts == 4) // TETRA
            k6searchblankedcellstetrax_(
              field.getSize(),
              field.begin(posx), field.begin(posy), field.begin(posz),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3), cn.begin(4),
              p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(), isNot,
              blankedCells.begin(), maskedl);
          else if (nelts == 6) // PENTA
            k6searchblankedcellspentax_(
              field.getSize(),
              field.begin(posx), field.begin(posy), field.begin(posz),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),
              cn.begin(4), cn.begin(5), cn.begin(6),
              p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(), isNot,
              blankedCells.begin(), maskedl);
          else if (nelts == 8) // HEXA
            k6searchblankedcellshexax_(
              field.getSize(),
              field.begin(posx), field.begin(posy), field.begin(posz),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3), cn.begin(4),
              cn.begin(5), cn.begin(6), cn.begin(7), cn.begin(8),
              p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(), isNot,
              blankedCells.begin(), maskedl);

          else {printf("Warning: searchForBlankedCellsUnstr: element type not valid.\n"); return -1;}
        }
        break;

      case 0: // node_in
        if (elevationDir == 2)
        {
          k6searchblankednodesx2d_(
            field.getSize(), field.begin(posx), field.begin(posy),
            p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
            p->Z.getSize(), p->Z.begin(), isNot,
            blankedCells.begin(), maskedl);
        }
        else // dim 3
        {
          k6searchblankednodesx_(field.getSize(), field.begin(posx),
                                 field.begin(posy),field.begin(posz),
                                 p->xmin, p->ymin, p->ni, p->nj,
                                 p->hi, p->hj, p->indir.begin(),
                                 p->Z.getSize(), p->Z.begin(), isNot,
                                 blankedCells.begin(), maskedl);
        }
        break;

      case -1: // cell_intersect_opt + depth=1
      case -2: // cell_intersect_opt + depth=2
        if (elevationDir == 2)
        {
          if (nelts == 3) // TRI
            k6searchblankedcellstrix2_(
              field.getSize(), field.begin(posx), field.begin(posy),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),
              p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(),
              blankedCells.begin(), maskedl);
          else if (nelts == 4)
            k6searchblankedcellsquadx2_(
              field.getSize(), field.begin(posx), field.begin(posy),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),cn.begin(4),
              p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(),
              blankedCells.begin(), maskedl);

          else {printf("Warning: searchForBlankedCellsUnstr: element type is invalid.\n"); return -1;}
        }
        else //3D
        {
          if (nelts == 4) // TETRA
            k6searchblankedcellstetrax2_(
              field.getSize(),
              field.begin(posx), field.begin(posy), field.begin(posz),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3), cn.begin(4),
              p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(),
              blankedCells.begin(), maskedl);
          else if (nelts == 6) // PENTA
            k6searchblankedcellspentax2_(
              field.getSize(),
              field.begin(posx), field.begin(posy), field.begin(posz),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),
              cn.begin(4), cn.begin(5), cn.begin(6),
              p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(),
              blankedCells.begin(), maskedl);

          else if (nelts == 8) // HEXA
            k6searchblankedcellshexax2_(
              field.getSize(),
              field.begin(posx), field.begin(posy), field.begin(posz),
              cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),cn.begin(4),
              cn.begin(5), cn.begin(6), cn.begin(7), cn.begin(8),
              p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(),
              blankedCells.begin(), maskedl);

          else {printf("Warning: searchForBlankedCellsUnstr: element type is invalid.\n"); return -1;}
        }
        break;
      default:
        printf("Warning: searchForBlankedCellUnstr: blankingCriteria value is invalid.\n");
        return -1;
    }
    masked = K_FUNC::E_max(maskedl, masked);
  }

//   if (delta > 0.)
//   {
//     return holeExpansionUnstr(elevationDir, blankingType, isNot, delta, planes,
//                               posx, posy, posz, field, cn, blankedCells);

//   }
  return (masked > 0);
}
//=============================================================================
// Expand hole of a distance _delta
//=============================================================================
E_Int K_CONNECTOR::holeExpansionUnstr(
  E_Int elevationDir, E_Int blankingType,
  E_Int isNot, E_Float delta,
  list<XRayPlane*>& planes,
  E_Int posx, E_Int posy, E_Int posz,
  FldArrayF& field, FldArrayI& cn,
  FldArrayI& cellNatFld)
{
  FldArrayI listOfInterpolatedPoints;
  FldArrayI diri; FldArrayI dirj;
  E_Int isMasked = 1;
  E_Int nelts = cn.getNfld();

  list<XRayPlane*>::iterator itr;
  //E_Int depth = 2;
  //if (blankingType == 1 || blankingType == -1 ) depth = 1;
  if (blankingType < 0) blankingType = 1;// cell_intersect_opt + delta trop long -> cell_intersect + delta

  // debug
  isMasked = 0;
  while (isMasked == 1)
  {
    // Calcul la liste des points interpoles
//     compListOfInterpolatedPointsUnstr(depth, cn,
//                                       cellNatFld, listOfInterpolatedPoints);
    isMasked = 0;

    for (itr = planes.begin(); itr != planes.end(); itr++)
    {
      XRayPlane* p = *itr;
      switch (blankingType)
      {
        case 1: // cell_intersect
          if (elevationDir == 2)
          {
            if (nelts == 3) // TRI
              k6searchblankedcellstrixd_(
                field.getSize(), field.begin(posx), field.begin(posy),
                cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),
                p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
                p->Z.getSize(), p->Z.begin(),
                listOfInterpolatedPoints.begin(),
                listOfInterpolatedPoints.getSize(),
                delta, isNot,
                cellNatFld.begin(), isMasked);
            else if (nelts == 4) // QUAD
              k6searchblankedcellsquadxd_(
                field.getSize(), field.begin(posx), field.begin(posy),
                cn.getSize(), cn.begin(1), cn.begin(2),
                cn.begin(3), cn.begin(4),
                p->xmin, p->ni, p->nj, p->hi, p->indir.begin(),
                p->Z.getSize(), p->Z.begin(),
                listOfInterpolatedPoints.begin(),
                listOfInterpolatedPoints.getSize(),
                delta, isNot, cellNatFld.begin(), isMasked);
            else {printf("Warning: holeExpansionUnstr: element type not valid.\n"); return -1;}
          }
          else
          {
            if (nelts == 4) // TETRA
              k6searchblankedcellstetraxd_(
                field.getSize(),
                field.begin(posx), field.begin(posy), field.begin(posz),
                cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),
                cn.begin(4),
                p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
                p->Z.getSize(), p->Z.begin(),
                listOfInterpolatedPoints.begin(),
                listOfInterpolatedPoints.getSize(),
                delta, isNot, cellNatFld.begin(), isMasked);
            else if (nelts == 6) // PENTA
              k6searchblankedcellspentaxd_(
                field.getSize(),
                field.begin(posx), field.begin(posy), field.begin(posz),
                cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),
                cn.begin(4), cn.begin(5), cn.begin(6),
                p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
                p->Z.getSize(), p->Z.begin(),
                listOfInterpolatedPoints.begin(),
                listOfInterpolatedPoints.getSize(),
                delta, isNot, cellNatFld.begin(), isMasked);
            else if (nelts == 8) // HEXA
              k6searchblankedcellshexaxd_(
                field.getSize(),
                field.begin(posx), field.begin(posy), field.begin(posz),
                cn.getSize(), cn.begin(1), cn.begin(2), cn.begin(3),
                cn.begin(4), cn.begin(5), cn.begin(6), cn.begin(7),
                cn.begin(8),
                p->xmin, p->ymin, p->ni, p->nj, p->hi, p->hj, p->indir.begin(),
                p->Z.getSize(), p->Z.begin(),
                listOfInterpolatedPoints.begin(),
                listOfInterpolatedPoints.getSize(),
                delta, isNot, cellNatFld.begin(), isMasked);
            else {printf("Warning: holeExpansionUnstr: element type not valid.\n"); return -1;}
          }
          break;

        case 0: // node_in
          if (elevationDir == 2)
          {
            k6searchblankednodesxd2d_(
              field.getSize(), field.begin(posx), field.begin(posy),
              p->xmin, p->ni, p->nj, p->hi,
              p->indir.begin(), p->Z.getSize(), p->Z.begin(),
              listOfInterpolatedPoints.begin(),
              listOfInterpolatedPoints.getSize(),
              diri.begin(),
              delta, isNot,
              cellNatFld.begin(), isMasked );
          }
          else
          {
            k6searchblankednodesxd_(
              field.getSize(), field.begin(posx),
              field.begin(posy),field.begin(posz),
              p->xmin, p->ymin, p->ni, p->nj,
              p->hi, p->hj, p->indir.begin(),
              p->Z.getSize(), p->Z.begin(),
              listOfInterpolatedPoints.begin(),
              listOfInterpolatedPoints.getSize(),
              diri.begin(), dirj.begin(),
              delta, isNot,
              cellNatFld.begin(), isMasked );
          }
        break;

        default:
        printf("Warning: holeExpansionUnstr: blankingType value is invalid.\n");
        return -1;
      }
    }
  }
  return 1;
}
//=============================================================================
// Expand hole of a distance _delta
//=============================================================================
E_Int K_CONNECTOR::holeExpansionStruct(E_Int elevationDir, E_Int blankingType,
                                       E_Int isNot, E_Float delta,
                                       list<XRayPlane*>& planes,
                                       E_Int ni, E_Int nj, E_Int nk,
                                       E_Int posx, E_Int posy, E_Int posz,
                                       FldArrayF& field,
                                       FldArrayI& cellNatFld)
{
  FldArrayI listOfInterpolatedPoints;
  FldArrayI diri; FldArrayI dirj;
  list<XRayPlane*>::iterator itr;
  E_Int imc = ni-1;
  E_Int jmc = nj-1; if (nj == 1) jmc = 1;
  E_Int kmc = nk-1; if (nk == 1) kmc = 1;
  E_Int depth = 2;
  if (blankingType == 1 || blankingType == -1) depth = 1;
  else if (blankingType == 0)
  {
    imc = ni; jmc = nj; kmc = nk; depth = 1;
  }
  E_Int isMasked = 1;
  while (isMasked == 1)
  {
    // Calcul la liste des points interpoles
    compListOfInterpolatedPoints(elevationDir, blankingType, depth, imc,jmc,kmc,
                                 cellNatFld, listOfInterpolatedPoints, diri, dirj);
    isMasked = 0;
    for (itr = planes.begin(); itr != planes.end(); itr++)
    {
      XRayPlane* p = *itr;
      switch (blankingType)
      {
        case 1: // cell_intersect
          if (elevationDir == 2)
          {
            k6searchblankedcellsxd2d_(
              ni, nj, nk, field.begin(posx),
              field.begin(posy),
              p->xmin, p->ni, p->nj, p->hi,
              p->indir.begin(),
              p->Z.getSize(), p->Z.begin(),
              listOfInterpolatedPoints.begin(),
              listOfInterpolatedPoints.getSize(),
              delta, isNot,
              cellNatFld.begin(), isMasked );
          }
          else
          {
            k6searchblankedcellsxd_(ni, nj, nk, field.begin(posx),
                                    field.begin(posy),field.begin(posz),
                                    p->xmin, p->ymin, p->ni, p->nj,
                                    p->hi, p->hj, p->indir.begin(),
                                    p->Z.getSize(), p->Z.begin(),
                                    listOfInterpolatedPoints.begin(),
                                    listOfInterpolatedPoints.getSize(),
                                    delta, isNot,
                                    cellNatFld.begin(), isMasked );
          }
          break;

        case 0: // node_in
          if (elevationDir == 2)
          {
            k6searchblankednodesxd2d_(field.getSize(), field.begin(posx),
                                      field.begin(posy),
                                      p->xmin, p->ni, p->nj, p->hi,
                                      p->indir.begin(),
                                      p->Z.getSize(), p->Z.begin(),
                                      listOfInterpolatedPoints.begin(),
                                      listOfInterpolatedPoints.getSize(),
                                      diri.begin(),
                                      delta, isNot,
                                      cellNatFld.begin(), isMasked );
          }
          else
          {
            k6searchblankednodesxd_(field.getSize(), field.begin(posx),
                                    field.begin(posy),field.begin(posz),
                                    p->xmin, p->ymin, p->ni, p->nj,
                                    p->hi, p->hj, p->indir.begin(),
                                    p->Z.getSize(), p->Z.begin(),
                                    listOfInterpolatedPoints.begin(),
                                    listOfInterpolatedPoints.getSize(),
                                    diri.begin(), dirj.begin(),
                                    delta, isNot,
                                    cellNatFld.begin(), isMasked );
          }
          break;

        case -2: // cell_intersect_opt + depth=2
          if (elevationDir == 2)
          {
            k6searchblankedcellsxd22d_( ni, nj, nk, field.begin(posx),
                                        field.begin(posy),
                                        p->xmin, p->ni, p->nj, p->hi,
                                        p->indir.begin(),
                                        p->Z.getSize(), p->Z.begin(),
                                        listOfInterpolatedPoints.begin(),
                                        listOfInterpolatedPoints.getSize(),
                                        delta,
                                        cellNatFld.begin(), isMasked );
          }
          else
          {
            k6searchblankedcellsxd2_(ni, nj, nk, field.begin(posx),
                                     field.begin(posy),field.begin(posz),
                                     p->xmin, p->ymin, p->ni, p->nj,
                                     p->hi, p->hj, p->indir.begin(),
                                     p->Z.getSize(), p->Z.begin(),
                                     listOfInterpolatedPoints.begin(),
                                     listOfInterpolatedPoints.getSize(),
                                     delta, cellNatFld.begin(), isMasked);
          }
          break;

        case -1: // cell_intersect_opt + depth=1
          if (elevationDir == 2)
            k6searchblankedcellsxd12d_(ni, nj, nk, field.begin(posx),
                                       field.begin(posy),
                                       p->xmin, p->ni, p->nj, p->hi,
                                       p->indir.begin(),
                                       p->Z.getSize(), p->Z.begin(),
                                       listOfInterpolatedPoints.begin(),
                                       listOfInterpolatedPoints.getSize(),
                                       delta, cellNatFld.begin(), isMasked);
          else
            k6searchblankedcellsxd1_(ni, nj, nk, field.begin(posx),
                                     field.begin(posy),field.begin(posz),
                                     p->xmin, p->ymin, p->ni, p->nj,
                                     p->hi, p->hj, p->indir.begin(),
                                     p->Z.getSize(), p->Z.begin(),
                                     listOfInterpolatedPoints.begin(),
                                     listOfInterpolatedPoints.getSize(),
                                     delta, cellNatFld.begin(), isMasked);
          break;

        default:
        printf("Warning: holeExpansionStruct: blankingType value is invalid.\n");
        return -1;
      }
    }

  }
  return 1;
}

//=============================================================================
// Compute the list of interpolated points
//=============================================================================
void K_CONNECTOR::compListOfInterpolatedPoints(
  E_Int elevationDir, E_Int type, E_Int depth,
  E_Int imc, E_Int jmc, E_Int kmc,
  FldArrayI& cellNatFld,
  FldArrayI& listOfInterpolatedPoints,
  FldArrayI& diri, FldArrayI& dirj)

{
  E_Int ncells = imc*jmc*kmc;

  FldArrayI cellNTemp = cellNatFld;
  E_Int dir = 0;
  searchMaskInterpolatedCellsStruct(imc, jmc, kmc, depth, dir, cellNatFld, cellNTemp);
  listOfInterpolatedPoints.malloc(ncells);
  E_Int* cellNt = cellNTemp.begin();
  E_Int* iP = listOfInterpolatedPoints.begin();

  E_Int n = 0;
  for (E_Int i = 0; i < ncells; i++)
  {
    if (cellNt[i] == 0) { iP[n] = i; n++; }
  }
  listOfInterpolatedPoints.reAlloc(n);
  // Search for sign of delta expansion in i and j directions
  iP = listOfInterpolatedPoints.begin();
  if (type == 0)
  {
    if (elevationDir == 2 ) // 2D
    {
      diri.malloc(n); diri.setAllValuesAtNull();
      for (E_Int i = 0; i < n; i++)
      {
        E_Int indN = iP[i];
        if (indN == 0)
        {
          if (cellNt[indN+1] == -1) diri[i] = 1;
        }
        else if ( indN == imc-1)
        {
          if (cellNt[indN-1] == -1) diri[i] = -1;
        }
        else
        {
          if (cellNt[indN-1] == -1) diri[i] = -1;
          else if (cellNt[indN+1] == -1) diri[i] = 1;
        }
      }
    }
    else
    {
      diri.malloc(n); diri.setAllValuesAtNull();
      dirj.malloc(n); dirj.setAllValuesAtNull();
      for (E_Int i = 0; i < n; i++)
      {
        E_Int indN = iP[i];
        if (indN == 0)
        {
          if (cellNt[indN+1] == -1 ) diri[i] = 1;
        }
        else if (indN == imc-1)
        {
          if (cellNt[indN-1] == -1 ) diri[i] = -1;
        }
        else
        {
          if (cellNt[indN-1] == -1 ) diri[i] = -1;
          else if (cellNt[indN+1] == -1 ) diri[i] = 1;
        }

        if (indN-imc < 0)
        {
          if (cellNt[indN+imc] == -1 ) dirj[i] = 1;
        }
        else if (indN+imc > ncells)
        {
          if (cellNt[indN-imc] == -1 ) dirj[i] = -1;
        }
        else
        {
          if (cellNt[indN-imc] == -1 ) dirj[i] = -1;
          else if (cellNt[indN+imc] == -1 ) dirj[i] = 1;
        }
      }
    }
  }
}
