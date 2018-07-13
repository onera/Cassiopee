/*    
    Copyright 2013-2018 Onera.

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

#include "Hexahedron.h"

namespace K_MESH
{
  
const E_Int K_MESH::Hexahedron::NB_NODES=8;
const E_Int K_MESH::Hexahedron::NB_TRIS=12;

void Hexahedron::triangulate(E_Int* target)
{
  // WARNING: connectT3 is Apended (not cleared upon entry)
  
  K_MESH::Quadrangle q4;
  
  for (size_t f = 0; f < 6; ++f)
  {
    getBoundary(f, q4);
    K_MESH::Quadrangle::triangulate(q4.nodes(), target + f*6);
  }
}

void Hexahedron::reorder_pgs(ngon_type& ng, const K_FLD::IntArray& F2E, E_Int i)
{
    std::map<E_Int,E_Int> glmap; // crd1 to 0-26 indexes
    E_Int nb_faces = ng.PHs.stride(i); 
    E_Int* faces = ng.PHs.get_facets_ptr(i);
    E_Int PGi = faces[0] - 1;
    E_Int* pN = ng.PGs.get_facets_ptr(PGi);

    glmap[*pN] = 0; // PHi(0,0) -> 0   

    if (F2E(1,PGi) == i) // for BOT, PH is the right element = well oriented
    { 
        for (int k = 1; k < 4; ++k){
            glmap[*(pN+k)] = k;
        }
    }
    else // wrong orientation : swap of 1 and 3
    { 
        glmap[*(pN+3)] = 1;
        glmap[*(pN+2)] = 2;
        glmap[*(pN+1)] = 3;
    }
    E_Int TopId(E_IDX_NONE),LeftId(E_IDX_NONE),RightId(E_IDX_NONE),FrontId(E_IDX_NONE),BackId(E_IDX_NONE);

    for (int k = 1; k < 6; ++k)
    {
        int count = 0;
        Vector_t<E_Boolean> commonNodes(4,false);
        E_Int testedPG = faces[k]-1;
        E_Int* pNode = ng.PGs.get_facets_ptr(testedPG);

        for (int j = 0; j < 4; ++j)
        {

            auto it = glmap.find(pNode[j]);
            if (it != glmap.end())
            {
                // found
                count++;
                commonNodes[it->second] = true;
            }
        }
        if (count == 0) // no common point, the ith PG is the TOP
            TopId = k;
        else if (commonNodes[0] && commonNodes[1])
            FrontId = k;
        else if (commonNodes[1] && commonNodes[2])
            RightId = k;
        else if (commonNodes[2] && commonNodes[3])
            BackId = k;
        else if (commonNodes[0] && commonNodes[3])
            LeftId = k;
    }
    
    E_Int mol[6];

    mol[0] = faces[0];
    mol[1] = faces[TopId];
    mol[2] = faces[LeftId];
    mol[3] = faces[RightId];
    mol[4] = faces[FrontId];
    mol[5] = faces[BackId];

    for (int i = 0; i < nb_faces; ++i)
    {
        faces[i] = mol[i];
    }  
}

void Hexahedron::get_internal(E_Int* nodes, E_Int* p)
{ 
    p[0] = nodes[9];  p[1] = nodes[11];  p[2] = nodes[16];  p[3] = nodes[14];  p[4] = nodes[12];  p[5] = nodes[25];  p[6] = nodes[17];  p[7] = nodes[23];  p[8] = nodes[26];
    p[9] = nodes[18];  p[10] = nodes[19];  p[11] = nodes[20];  p[12] = nodes[21];  p[13] = nodes[22];  p[14] = nodes[23];  p[15] = nodes[24];  p[16] = nodes[25];  p[17] = nodes[26];
    p[18] = nodes[8];  p[19] = nodes[10];  p[20] = nodes[15];  p[21] = nodes[13];  p[22] = nodes[12];  p[23] = nodes[24];  p[24] = nodes[17];  p[25] = nodes[22];  p[26] = nodes[26];    
}

void Hexahedron::get_orient(const ngon_type& ng, const K_FLD::IntArray& F2E, E_Int PHi, E_Int* PHi_orient)
{
    const E_Int* p = ng.PHs.get_facets_ptr(PHi);
    
    for (int i = 0; i < 6; ++i)
      PHi_orient[i] = (F2E(1,p[i]-1) == PHi) ? -1 : 1;
}

bool Hexahedron::pt_is_inside(const ngon_type& ng, const K_FLD::FloatArray& crd1, E_Int PHi, const E_Int* PHi_orient, const E_Float* pt)
{
    const E_Int* p = ng.PHs.get_facets_ptr(PHi);
    
    for (int i = 0; i < 6; i++)
    {
        const E_Int* pN = ng.PGs.get_facets_ptr(p[i]-1);
        E_Float det = K_FUNC::zzdet4(crd1.col(pN[0]-1), crd1.col(pN[1]-1), crd1.col(pN[2]-1), pt);

        if ((PHi_orient[i] == 1) && (det > E_EPSILON) ) return false;
        else if ((PHi_orient[i] == -1) && (det < E_EPSILON) ) return false;
    }
    return true;
}


}
