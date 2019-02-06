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

#ifndef _CPLOT_UNSTRUCTZONE_H_
#define _CPLOT_UNSTRUCTZONE_H_

//#define __SHADERS__ // A ENLEVER

#include "Zone.h"
//#if defined(__SHADERS__)
//#  include "CPlot/Shaders/ShaderBufferObject.hpp"
//#endif

/* Define a unstructured zone */
class UnstructZone : public Zone
{
  public:
  UnstructZone( CPlotState* states, ZoneImpl* impl );
    virtual ~UnstructZone();

    void compNorm();

  public:
    int np;                       // number of points
    int ne;                       // number of elements
    enum enumType {
      NODE = 0, BAR, TRI, QUAD, TETRA, PENTA, PYRA, HEXA, NGON, ENDTYPE
    };
    int eltType;                  // 0:   NODE
                                  // 1: BAR
                                  // 2: TRI
                                  // 3: QUAD
                                  // 4: TETRA
                                  // 7: HEXA
                                  // 10: NGON
    int eltSize;                  // nbre de noeuds par elt
                                  // sauf pour NGONS (=1) 
    int* connect;                 // connectivity (ne*eltSize sauf pour NGONS)
                                  // Seulement pour NGONS
    int* posFaces;                // position de la face no i dans connect
    int nelts1D;                  // nombre d'elements 1D
    int nelts2D;                  // nombre d'elements 2D
    int* posElts1D;               // position des elements 1D dans connect
    int* posElts2D;               // position des elements 2D dans connect
    int  _is_high_order;          // Dit si on est en presence d'elements high order ou non
};
#endif
