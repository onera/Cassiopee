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
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef __MASKGEN_CPP__
#define __MASKGEN_CPP__

#include "maskGen.h"

#include <algorithm>

#include "CompGeom/compGeom.h"

#ifdef DEBUG_MASK
#include <fstream>
bool K_CONNECTOR::maskGen::dbg_switch = false;
#include "IO/io.h"
#endif

namespace K_CONNECTOR
{  
  ///
  maskGen::~maskGen()
  {
    //std::cout << "clearing the tree" << std::endl;
    for (size_t i = 0; i < _boxes.size(); ++i)
      delete _boxes[i];
    delete _tree;
    if (_reoriented)
    {
      typedef AConnec_t::array_type acon_t;
      acon_t* ptr = const_cast<acon_t*>(&_connT4->array());
      delete ptr;
    }
    delete _connT4;
    
    if (_kdtree) delete _kdtree; 
  }
 
  ///
  E_Int maskGen::blank
  (const K_FLD::FldArrayF& coord, E_Int px, E_Int py, E_Int pz, K_FLD::FldArrayI& isBlanked, E_Int cellnval, bool overwrite)
  {
#ifdef DEBUG_MASK
  std::cout << " INPUT BLANKING " << std::endl;
  std::cout << "coord cols/rows : " << coord.getSize() << "/" << coord.getNfld() << std::endl;
  std::cout << "px py pz " << px << " " << py << " " << pz << std::endl;
  std::cout << "isBlanked cols/rows : " <<  isBlanked.getSize() << "/" << isBlanked.getNfld() << std::endl;
#endif
  
    K_FLD::ArrayAccessor<K_FLD::FldArrayF> cA(coord, px, py, pz);
    
    Vector_t<E_Int> b;
    if (_ELType == TH4)
      __blank<K_MESH::Tetrahedron>(cA, b);
    else if (_ELType == TRI)
      __blank<K_MESH::Triangle>(cA, b);

   if (overwrite) isBlanked.resize(0, 1);
     
    size_t bsz(b.size());
    if ((size_t)isBlanked.getSize() < bsz)
      isBlanked.resize(1, b.size(), VISIBLE);
  
    for (size_t i = 0; i < bsz; ++i)
      if (b[i] == MASKED) 
        isBlanked[i]=cellnval;
    
#ifdef DEBUG_MASK
    E_Int count_blk(0);
    for (size_t i = 0; i < bsz; ++i)
    {
      if (isBlanked[i]==cellnval)
        ++count_blk;
    }
    std::cout << "nb of blanked : " << count_blk << " over " << bsz << std::endl;
#endif
  
    return 0;
  }
   
  ///
  template <typename T>
  void maskGen::__blank
  (const ACoord_t& coord, Vector_t<E_Int>& isBlanked)
  {
    isBlanked.clear();
    isBlanked.resize(coord.size(), VISIBLE);
    
    E_Float P[3];
    size_t sz(coord.size());
    Vector_t<E_Int> caught_boxes;
    
#ifdef DEBUG_MASK
    E_Int nb_box_av(0), nb_box_max(0), nb_box_min(_connT4->size());
#endif

#ifndef DEBUG_MASK
#pragma omp parallel for private(caught_boxes, P)
#endif
    for (size_t i = 0; i < sz; ++i)
    {
      coord.getEntry(i, P);
      
#ifdef DEBUG_MASK
  //if (i == 10+(99*11)+(99*99*12))
    if (i == 101111)
    {
      maskGen::dbg_switch=true;
      //bool b = is_inside<T>(P, caught_boxes);
      //maskGen::dbg_switch=false;
      //std::cout << "is [" << P[0] << "," << P[1] << "," << P[2] << "] in ?? : " << b << std::endl;
    }
#endif
      
      //if (i == 111011)
        if (is_inside<T>(P, caught_boxes))
          isBlanked[i] = MASKED;
      
#ifdef DEBUG_MASK
      nb_box_av += caught_boxes.size();
      if (caught_boxes.size() < nb_box_min) nb_box_min = caught_boxes.size();
      if (caught_boxes.size() > nb_box_max) nb_box_max = caught_boxes.size();
#endif
    }
    
#ifdef DEBUG_MASK
    nb_box_av /= sz;
    std::cout << " nb of min/av/max boxes : " << nb_box_min << "/" << nb_box_av << "/" << nb_box_max << std::endl;
#endif
    
    //std::cout << " STATS : direct/almost/heavy : " << count_direct << "/" << count_almost_direct << "/" << count_heavy << std::endl;
    
  }
  
  ///
  void maskGen::__init()
  {
    //std::cout << "init : box creation" << std::endl;
     // compute T4 boxes
    __create_boxes(_boxes);
    //std::cout << "init : tree creation" << std::endl;
    // build the tree.
    _tree = new K_SEARCH::BbTree3D(_boxes, _tolerance);
    
    //std::cout << "init : reorient" << std::endl;
    /*{
    K_FLD::FloatArray pos(_coordT4.array());
    K_FLD::IntArray connect(_connT4->array(), _connT4->shift());
    meshIO::write("/home/slandier/tmp/incorrect.mesh", pos, connect);
    }*/
     // reorient
    if (_ELType == TH4)
      _reoriented = __reorient<K_MESH::Tetrahedron>();
    //else if (_ELType == TRI)
    //  _reoriented = __reorient<K_MESH::Triangle>();
    
    //K_FLD::FloatArray pos(_coordT4.array());
    //K_FLD::IntArray connect(_connT4->array(), _connT4->shift());
    //meshIO::write("/home/slandier/tmp/correct.mesh", pos, connect);
    
    if (_ELType == TRI)
    {
      _kdtree= new K_SEARCH::KdTree<K_FLD::FldArrayF>(_coordT4, _tolerance);
      
      _ancestors.resize(_coordT4.size(), E_IDX_NONE);
      for (E_Int i=0; i < _connT4->size(); ++i)
      {
        _ancestors[_connT4->getVal(i,0)]=i;
        _ancestors[_connT4->getVal(i,1)]=i;
        _ancestors[_connT4->getVal(i,2)]=i;
      }
      
      E_Int T3[3], sz;
      sz = _connT4->size();
      _isoG.resize(3, sz);
      for (size_t i=0;i< _ancestors.size(); ++i)
      {
        if (_ancestors[i] == E_IDX_NONE)
          continue;
        _connT4->getEntry(_ancestors[i], T3);
        K_MESH::Triangle::isoG(_coordT4, T3, _isoG.col(_ancestors[i]));
      }
      
      // Pre-compue normals
      _normals.resize(3, sz);
      E_Float Q0[3], Q1[3], Q2[3];
      E_Int Ti[3];
      for (E_Int i = 0; i < sz; ++i)
      {
        _connT4->getEntry(i, Ti);
        _coordT4.getEntry(Ti[0], Q0);
        _coordT4.getEntry(Ti[1], Q1);
        _coordT4.getEntry(Ti[2], Q2);
      
        K_MESH::Triangle::normal(Q0,Q1,Q2, _normals.col(i));
        K_FUNC::normalize<3>(_normals.col(i));
      }
      
    }
    
#ifdef DEBUG_MASK
    /*std::ofstream f("/home/slandier/tmp/MASK/tree.txt");
    std::cout << " ze tree " << std::endl;
    f << _tree->_tree << std::endl;
    f.close();*/
#endif
  }
  
 /* ///
  void maskGen::__create_boxes(Vector_t<K_SEARCH::BBox3D*>& boxes)
  {
    size_t nb_elts(_connT4->size());
    const E_Int NB_NODES = _connT4->stride();
    
#ifdef DEBUG_MASK
    std::cout << "__create_boxes : nb of elements : " << nb_elts <<std::endl;
    std::cout << "__create_boxes : NB_NODES : " << NB_NODES <<std::endl;
    
    std::ofstream fc("/home/slandier/tmp/MASK/coordT4.txt");
    fc << _coordT4.array() << std::endl;
    fc.close();
    
    std::ofstream ft("/home/slandier/tmp/MASK/connectT4.txt");
    ft << _connT4->array() << std::endl;
    ft.close();
    
#endif
    
    E_Int* elt = new E_Int[NB_NODES];
    
    for (size_t i = 0; i < nb_elts; ++i)
    {
      _connT4->getEntry(i, elt);
      
      K_SEARCH::BBox3D* box = new K_SEARCH::BBox3D(_coordT4, elt, NB_NODES);
      boxes.push_back(box);
    }
    
    delete [] elt;
    
#ifdef DEBUG_MASK
    std::ofstream f("/home/slandier/tmp/MASK/boxes.txt");
    for (size_t i = 0; i < nb_elts; ++i)
    {
      f <<boxes[i]->minB[0] << " " << boxes[i]->minB[1] << " " << boxes[i]->minB[2] << std::endl;
      f <<boxes[i]->maxB[0] << " " << boxes[i]->maxB[1] << " " << boxes[i]->maxB[2] << std::endl;
      f << std::endl;
    }
    f.close();
#endif
  }*/
  
  ///
  inline void maskGen::__create_boxes(Vector_t<K_SEARCH::BBox3D*>& boxes)
  {
    size_t nb_elts(_connT4->size());
    
#ifdef DEBUG_MASK
    const E_Int NB_NODES = _connT4->stride();
    std::cout << "__create_boxes : nb of elements : " << nb_elts <<std::endl;
    std::cout << "__create_boxes : NB_NODES : " << NB_NODES <<std::endl;
    
    std::ofstream fc("/home/slandier/tmp/MASK/coordT4.txt");
    fc << _coordT4.array() << std::endl;
    fc.close();
    
    std::ofstream ft("/home/slandier/tmp/MASK/connectT4.txt");
    ft << _connT4->array() << std::endl;
    ft.close();
    
#endif
    
    boxes.resize(nb_elts);// liste des bbox de ts les elements de a2
    K_FLD::FldArrayF bbox(nb_elts,6);// xmin, ymin, zmin, xmax, ymax, zmax
    K_FLD::FldArrayF& arr = const_cast<K_FLD::FldArrayF&>(_coordT4.array());
    K_FLD::FldArrayI& con = const_cast<K_FLD::FldArrayI&>(_connT4->array());
    K_COMPGEOM::boundingBoxOfUnstrCells(
      con, arr.begin(_coordT4.posX(0)),
      arr.begin(_coordT4.posX(1)), arr.begin(_coordT4.posX(2)),  
      bbox);
    E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
    E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
    E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
    E_Float minB[3];  E_Float maxB[3];
    for (size_t et = 0; et < nb_elts; et++)
    {
      minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
      maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
      boxes[et] = new K_SEARCH::BBox3D(minB, maxB);
    }
    
#ifdef DEBUG_MASK
    /*std::ofstream f("/home/slandier/tmp/MASK/boxes.txt");
    for (size_t i = 0; i < nb_elts; ++i)
    {
      f <<boxes[i]->minB[0] << " " << boxes[i]->minB[1] << " " << boxes[i]->minB[2] << std::endl;
      f <<boxes[i]->maxB[0] << " " << boxes[i]->maxB[1] << " " << boxes[i]->maxB[2] << std::endl;
      f << std::endl;
    }
    f.close();*/
#endif
  }
}
#endif
