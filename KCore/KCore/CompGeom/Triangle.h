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

//=============================================================================
/* Class defining a triangle  */
//=============================================================================
class Triangle
{
  public:    
    ///+ 1- Constructors / Destructor
    /** Constructor :
     IN: indA: index of point A
     IN: indB: for point B
     IN: indC: for point C */
    Triangle(E_Float coefa, E_Float coefb, 
             E_Float coefc, E_Float coefd,
             E_Int indA, E_Int indB, E_Int indC, K_FLD::FldArrayF& coord);
    
    /** Destructor */
    ~Triangle();

    /* Get the 3 vertices indices in the coord array*/
    void getVertexIndices(E_Int& indA, E_Int& indB, E_Int& indC);

    /* Return OC^2 + r^2, r radius of circum circle */
    E_Float getDistMaxOfCC();

    /* Return true if point of index ind in coord array is in circum circle */
    E_Boolean isPointInCircumCircle(E_Int ind, K_FLD::FldArrayF& coord,
                                    E_Float tol=1.e-14);

  private:
    /* Compute the circum circle radius and center*/
    void compCircumCircle(E_Float coefa, E_Float coefb, 
                          E_Float coefc, E_Float coefd,
                          K_FLD::FldArrayF& coord);
  
  private:
    E_Int _indA; // indices of the 3 vertices in coord array
    E_Int _indB;
    E_Int _indC;

    E_Float _distmax; //distance max for completed status criterium
    E_Float _radius2; // square radius of the circumcircle of the triangle ABC
    E_Float _coordCC[3]; // coordinates of the circumcircle center of ABC 
};
//-----------------------------------------------------------------------------
// inline functions
//-----------------------------------------------------------------------------
inline void K_COMPGEOM::Triangle::getVertexIndices(
  E_Int& indA, E_Int& indB, E_Int& indC)
{
  indA = _indA;
  indB = _indB;
  indC = _indC;
}
//-----------------------------------------------------------------------------
inline E_Float K_COMPGEOM::Triangle::getDistMaxOfCC()
{
  return _distmax;
}


//========================= KCore/CompGeom/Triangle.h =========================
