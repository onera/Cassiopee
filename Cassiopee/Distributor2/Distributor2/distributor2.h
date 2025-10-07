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

# ifndef _DISTRIBUTOR2_DISTRIBUTOR2_H_
# define _DISTRIBUTOR2_DISTRIBUTOR2_H_

# include "kcore.h"
# include <vector>

namespace K_DISTRIBUTOR2
{
E_Float eval(E_Int nb, E_Int NProc, E_Float meanPtsPerProc,
             std::vector<E_Float>& solver, std::vector<E_Float>& latence,
             std::vector<E_Float>& comSpeed, E_Int* com, E_Int* comd, E_Int sizeComd,
             K_FLD::FldArrayF& nbPtsPerProcs, std::vector<E_Float>& nbPts,
             E_Int* nicePeople );
void genetic(std::vector<E_Float>& nbPts, std::vector<E_Int>& setBlocks,
             E_Int NProc, E_Int* com, E_Int* comd, E_Int sizeComd,
             std::vector<E_Float>& solver,
             std::vector<E_Float>& latence, std::vector<E_Float>& comSpeed,
             E_Int param,
             std::vector<E_Int>& out, E_Float& meanPtsPerProc, E_Float& varMin,
             E_Float& varMax, E_Float& varRMS, E_Int& nptsCom,
             E_Float& volRatio, E_Float& bestAdapt);
void graph(std::vector<E_Float>& nbPts, std::vector<E_Int>& setBlocks,
           E_Int NProc, E_Int* com, E_Int* comd, E_Int sizeComd,
           std::vector<E_Float>& solver,
           std::vector<E_Float>& latence, std::vector<E_Float>& comSpeed,
           E_Int param,
           std::vector<E_Int>& out, E_Float& meanPtsPerProc, E_Float& varMin,
           E_Float& varMax, E_Float& varRMS, E_Int& nptsCom,
           E_Float& volRatio, E_Float& bestAdapt);
void gradient(std::vector<E_Float>& nbPts, std::vector<E_Int>& setBlocks,
              E_Int NProc, E_Int* com, E_Int* comd, E_Int sizeComd,
              std::vector<E_Float>& solver,
              std::vector<E_Float>& latence, std::vector<E_Float>& comSpeed,
              E_Int param, std::vector<E_Int>& out,
              E_Float& meanPtsPerProc, E_Float& varMin,
              E_Float& varMax, E_Float& varRMS, E_Int& nptsCom,
              E_Float& volRatio, E_Float& bestAdapt);
void stats(std::vector<E_Float>& nbPts, E_Int NProc,
           E_Int* com, E_Int* comd, E_Int sizeComd,
           std::vector<E_Int>& out,
           E_Int& empty, E_Float& varMin, E_Float& varMax, E_Float& varRMS,
           E_Float& volRatio);

PyObject* distribute(PyObject* self, PyObject* args);
}
#endif
