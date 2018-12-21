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
#ifndef _POST_ZIPPER_ZIPLIB_H
#define _POST_ZIPPER_ZIPLIB_H

# include "SingleSegment.h"
# include "Pocket.h"
# include <vector>
# define FldArrayF K_FLD::FldArrayF

/** Test if edge1 (field1) and edge2 (field2) coming from different
    triangles are matching */
E_Boolean testIfEdgesAreMatching(E_Float matchTol,
                                 FldArrayF& field1, 
                                 FldArrayF& field2);

/** Test if pt1 (field1(i1,:)) and one extremity of field2 coming 
    from different singleSegments are matching. If true return i2. 
    if not found return -1  */
E_Int testIfPointsAreMatching(E_Float matchTol, E_Int i1, 
                              FldArrayF& field1, 
                              FldArrayF& field2);

/** Write the single segments to be pocketted */
void writeSingleSegments(std::vector<SingleSegment*>& singleSegments);

/** Build pockets
    In: matchTol : matching Tolerance
        singleSegments : vector of segments not set in a pair 
    Out : vector of pockets built by consecutive singleSegments*/
void buildPocket(E_Float matchTol,
                 std::vector<SingleSegment*>& singleSegments,
                 std::vector<Pocket*>& pockets);

/** Add segment SORTED in list of pockets */
void addSegmentInGlobalTab(E_Int istart, FldArrayF& field,
                           FldArrayF& globalTab, E_Int& cnt);

/* Given an array of size sizeIni*nfld erase mutiply defined pts
   and resize the array*/
void eraseDoublePts(E_Int sizeIni, E_Int nfld, FldArrayF& globalTab);


/* project point P on plane ABC. Return the coordinates of the projection H
   If projection is not possible return False
   IN : (x,y,z) : coordinates of point to be projected
   IN : (x0,y0,z0),(x1,y1,z1) (x2,y2,z2) coordinates of 3 points of the plane
   OUT : xp, yp, zp coordinates of the projection of (x,y,z) on plane (ABC)
   OUT : distProj : distance between (x,y,z) and (xp,yp,zp)
*/
E_Boolean projectPointOnPlane(E_Float x, E_Float y, E_Float z,
                              E_Float x0, E_Float y0, E_Float z0,
                              E_Float x1, E_Float y1, E_Float z1,
                              E_Float x2, E_Float y2, E_Float z2,
                              E_Float& xp, E_Float& yp, E_Float& zp,
                              E_Float& distProj);
# undef FldArrayF
#endif
