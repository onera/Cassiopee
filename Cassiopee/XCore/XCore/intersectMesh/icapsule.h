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
#pragma once

#include "mesh.h"

struct ICapsule {
    IMesh M;
    std::vector<IMesh> Ss;

    E_Float NEAR_EDGE_TOL = 1e-3;
    E_Float NEAR_VERTEX_TOL = 1e-3;

    ICapsule(){}

    ICapsule(const Karray &marray, const std::vector<Karray> &sarrays,
        const std::vector<E_Float *> &ptags);
    
    static std::vector<PointLoc> refine(Smesh &Mf,
        std::set<E_Int> &mfids, Smesh &Sf);
};
