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

#ifndef _CPLOT_UNSTRUCTZONE_H_
#define _CPLOT_UNSTRUCTZONE_H_

#include "Zone.h"
#include <vector>

/* Define a unstructured zone */
class UnstructZone : public Zone
{
  public:
    UnstructZone(CPlotState* states, ZoneImpl* impl);
    virtual ~UnstructZone();
    void compNorm();

  public:
    /* Types of basic elements */
    enum enumType { NODE=0, BAR, TRI, QUAD, TETRA, PENTA, PYRA, HEXA, EMPTY1, EMPTY2, NGON, ENDTYPE };

    E_Int np;                       // number of points
    E_Int ne;                       // total number of elements
    /* An unstruct zone can have multiple connectivities but they must be of the same dimension 
       For NGONs, zone has only one connectivity. */
    std::vector<E_Int> nec;         // number of elements for each connectivity
    std::vector<int> eltType;     // element type for each connectivity (see enum)
                                  // 0: NODE, 1: BAR, 2: TRI, 3: QUAD, 4: TETRA
                                  // 5: PENTA, 6: PYRA, 7: HEXA
                                  // 10: NGON
    std::vector<E_Int> eltSize;   // nbre of nodes for each element
    E_Int _is_high_order;         // Dit si on est en presence d'elements high order ou non (all connects)
                                  // sauf pour NGONS (=1) 
    std::vector<E_Int*> connect;  // connectivities (ne*eltSize sauf pour NGONS)
                                  // Seulement pour NGONS:
    E_Int* posFaces;              // position de la face no i dans connect
    E_Int nelts1D;                // nombre d'elements 1D
    E_Int nelts2D;                // nombre d'elements 2D
    E_Int* posElts1D;             // position des elements 1D dans connect
    E_Int* posElts2D;             // position des elements 2D dans connect
};
#endif
