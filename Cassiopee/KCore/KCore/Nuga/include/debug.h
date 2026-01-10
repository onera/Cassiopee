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

#ifndef __DEBUG_H__
#define	__DEBUG_H__

#include <fstream>

class Gen_debug
{
public:
  static void add_box_to_extract_file(const K_FLD::FloatArray& crd, std::ofstream& extract_file, std::vector<int>* indices=0)
  {
    K_SEARCH::BBox3D box;
    if (!indices)
      box.compute(K_FLD::ArrayAccessor<K_FLD::FloatArray >(crd));
    else
      box.compute(K_FLD::ArrayAccessor<K_FLD::FloatArray >(crd), *indices);

    /*extract_file << "boxes.push_back(K_SEARCH::BBox3D());" << std::endl;
    extract_file << "boxes[boxes.size()-1].minB[0]=" << box.minB[0] << ";"<< std::endl;
    extract_file << "boxes[boxes.size()-1].minB[1]=" << box.minB[1] << ";"<< std::endl;
    extract_file << "boxes[boxes.size()-1].minB[2]=" << box.minB[2] << ";"<< std::endl;
    extract_file << "boxes[boxes.size()-1].maxB[0]=" << box.maxB[0] << ";"<< std::endl;
    extract_file << "boxes[boxes.size()-1].maxB[1]=" << box.maxB[1] << ";"<< std::endl;
    extract_file << "boxes[boxes.size()-1].maxB[2]=" << box.maxB[2] << ";"<< std::endl << std::endl;*/
    
    extract_file << box.minB[0] << ";" << box.minB[1] << ";" << box.minB[2] << ";" << box.maxB[0] << ";" << box.maxB[1] << ";" << box.maxB[2] << ";" << std::endl;
    
  }
  
  static void get_PHT3_points
  (E_Int PHT3i, const K_FLD::IntArray& connectT3, 
   const std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s, 
   const K_FLD::FloatArray& icrd, K_FLD::FloatArray& ocrd)
{
  std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >::const_iterator it = PH_to_PGT3s.find(PHT3i);
  if (it == PH_to_PGT3s.end())
    return;
  
  Vector_t<E_Int> keep(icrd.cols(), false);
  
  const std::map<E_Int, Vector_t<E_Int> >& Aggs = it->second;
  std::map<E_Int, Vector_t<E_Int> >::const_iterator itPG = Aggs.begin();
  for (; itPG != Aggs.end(); ++itPG)
  {
    const E_Int& PGi = itPG->first;
    const Vector_t<E_Int> & T3s = itPG->second;
    for (size_t i=0; i < T3s.size(); ++i)
    {
      keep[connectT3(0, T3s[i])]=keep[connectT3(1, T3s[i])]=keep[connectT3(2, T3s[i])]=true;
    }
  }
  
  ocrd=icrd;
  Vector_t<E_Int> nids;
  K_FLD::FloatArray::compact(ocrd, keep, nids);
}
};

#endif
