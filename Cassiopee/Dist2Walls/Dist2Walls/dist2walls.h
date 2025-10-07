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

# ifndef _DIST2WALLS_DIST2WALLS_H_
# define _DIST2WALLS_DIST2WALLS_H_

# include "kcore.h"

namespace K_DIST2WALLS
{
  void computeMininterf(
    std::vector<E_Int>& ncellst,
    E_Int posx, E_Int posy, E_Int posz,
    std::vector<K_FLD::FldArrayF*>& fields,
    std::vector<E_Int>& posxv, std::vector<E_Int>& posyv, std::vector<E_Int>& poszv,
    std::vector<E_Int>& poscv, std::vector<K_FLD::FldArrayF*>& fieldsw,
    std::vector<K_FLD::FldArrayF*>& distances,
    std::vector<K_FLD::FldArrayI*>& cntw,E_Int isminortho);
  void computeMininterfSigned(
    std::vector<E_Int>& ncellst,
    E_Int posx, E_Int posy, E_Int posz,
    std::vector<K_FLD::FldArrayF*>& fields,
    std::vector<E_Int>& posxv, std::vector<E_Int>& posyv, std::vector<E_Int>& poszv,
    std::vector<E_Int>& poscv, std::vector<K_FLD::FldArrayF*>& fieldsw,
    E_Int possx, E_Int possy, E_Int possz,
    std::vector<K_FLD::FldArrayF*>& distances);
  void computeOrthoDist(
    std::vector<E_Int>& ncellst,
    E_Int posx, E_Int posy, E_Int posz, std::vector<E_Int>& posflag,
    std::vector<K_FLD::FldArrayF*>& fields,
    std::vector<E_Int>& posxv, std::vector<E_Int>& posyv, std::vector<E_Int>& poszv,
    std::vector<E_Int>& poscv,
    std::vector<K_FLD::FldArrayF*>& fieldsw,
    std::vector<K_FLD::FldArrayI*>& cntw,
    std::vector<K_FLD::FldArrayF*>& distances,E_Int isminortho,E_Int isIBM_F1,E_Float dTarget);
  void computeSignedOrthoDist(
    std::vector<E_Int>& ncellst,
    E_Int posx, E_Int posy, E_Int posz,
    std::vector<K_FLD::FldArrayF*>& fields,
    std::vector<E_Int>& posxv, std::vector<E_Int>& posyv, std::vector<E_Int>& poszv,
    std::vector<E_Int>& poscv,
    std::vector<K_FLD::FldArrayF*>& fieldsw, std::vector<K_FLD::FldArrayI*>& cntw,
    E_Int possx, E_Int possy, E_Int possz,
    std::vector<K_FLD::FldArrayF*>& distances);

  PyObject* distance2Walls(PyObject* self, PyObject* args);
  PyObject* distance2WallsSigned(PyObject* self, PyObject* args);
  PyObject* distance2WallsOrtho(PyObject* self, PyObject* args);
  PyObject* distance2WallsOrthoSigned(PyObject* self, PyObject* args);
  PyObject* eikonal(PyObject* self, PyObject* args);
}
#endif
