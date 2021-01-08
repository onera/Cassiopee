/*    
    Copyright 2013-2020 Onera.

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

#include "occ.h"
#include <TopoDS_Face.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRep_Tool.hxx>

// Project coords on CAD face
void projectOnFace__(E_Int npts, E_Float* px, E_Float* py, E_Float* pz, TopoDS_Face& F)
{
  gp_Pnt Point;
  Handle(Geom_Surface) face = BRep_Tool::Surface(F);
  
  for (E_Int i=0; i < npts; i++)
  {
    Point.SetCoord(px[i], py[i], pz[i]);
    GeomAPI_ProjectPointOnSurf o(Point, face, Extrema_ExtAlgo_Tree);
    gp_Pnt Pj = o.NearestPoint();
    //printf("projection %f %f %f -> %f %f %f\n",x,y,z,Pj.X(),Pj.Y(),Pj.Z());
    px[i] = Pj.X(); py[i] = Pj.Y(); pz[i] = Pj.Z();
  }
}
