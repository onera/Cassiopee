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
//Authors : Sam Landier (sam.landier@onera.fr)


#include "Nuga/include/Tetrahedron.h"

namespace K_MESH
{

///
bool Tetrahedron::is_inside
(const E_Float* Ni, const E_Float* Nj, const E_Float* Nk, const E_Float* Nl, const E_Float* pt)
{
  if (NUGA::zzdet4(Ni, Nj, Nl, pt) > EPSILON)
    return false;
  if (NUGA::zzdet4(Nj, Nk, Nl, pt) > EPSILON)
    return false;
  if (NUGA::zzdet4(Nk, Ni, Nl, pt) > EPSILON)
    return false;
  if (NUGA::zzdet4(Ni, Nk, Nj, pt) > EPSILON)
    return false;
  return true;
}

void Tetrahedron::triangulate(E_Int* target)
{  
  target[0]=_nodes[0]; target[1]=_nodes[1]; target[2]=_nodes[3];
  target[3]=_nodes[0]; target[4]=_nodes[3]; target[5]=_nodes[2];
  target[6]=_nodes[0]; target[7]=_nodes[2]; target[8]=_nodes[1];
  target[9]=_nodes[1]; target[10]=_nodes[2]; target[11]=_nodes[3];
}

E_Float Tetrahedron::quality(const K_FLD::FloatArray& crd, E_Float* Vol)
{
  E_Float S(.0), V;
  E_Int i0, i1, i2;
  const E_Float *p0, *p1, *p2, *p3;
  E_Float Q=12.*sqrt(1.5);
    
  for (int i=0; i<4; i++)
  {
    i0= node(i%4);
    i1= node((i+1)%4);
    i2= node((i+2)%4);
    p0= crd.col(i0);
    p1= crd.col(i1);
    p2= crd.col(i2);
    
    S += K_MESH::Triangle::surface<3>(p0,p1,p2); 
  }
  
  E_Float Lmin2,Lmax2;
  edge_length_extrema(crd, Lmin2, Lmax2);
    
  p0= crd.col(node(0));    
  p1= crd.col(node(1));
  p2= crd.col(node(2));
  p3= crd.col(node(3));
  
  *Vol = V =::fabs(volume(p0, p1, p2, p3));
  
  return Q*V/(sqrt(Lmax2)*S);
}

///
void Tetrahedron::edge_length_extrema(const K_FLD::FloatArray& crd, E_Float& Lmin2, E_Float& Lmax2)
{
  Lmin2 = NUGA::FLOAT_MAX;
  Lmax2 = -1.;
  
  for (size_t i=0; i < 4; ++i)
  {
    const E_Float* Pi = crd.col(node(i));
    
    for (size_t j=i+1; j < 4; ++j)
    {
      const E_Float* Pj = crd.col(node(j));
      
      E_Float L2 = NUGA::sqrDistance(Pi,Pj, 3);
      
      Lmin2 = MIN(Lmin2, L2);
      Lmax2 = MAX(Lmax2, L2);
    }
  }
}

}


