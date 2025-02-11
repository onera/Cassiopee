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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef __COLLIDER_H__
#define	__COLLIDER_H__

#include <vector>
//#include <algorithm>
#define Vector_t std::vector

#include "ArrayAccessor.h"
#include "ArrayWriter.h"
#include "BbTree.h"
#include "Predicates.h"

#ifdef COLLIDE_DBG
//#include "Search/chrono.h"
#include "IO/DynArrayIO.h"
#endif

// For specialisations
#include "Nuga/include/Hexahedron.h"


#define TEMPLATE_COORD_CONNECT_DIM_ELT_PRED template <typename Coordinate_t, typename Connectivity_t, short DIM, typename Element_t, template<class, class> class Predicate_t >
#define COLLIDER Collider<Coordinate_t, Connectivity_t, DIM, Element_t, Predicate_t>

namespace K_CONNECT
{  
TEMPLATE_COORD_CONNECT_DIM_ELT_PRED
class Collider
{
public:
  
  typedef K_FLD::ArrayAccessor<Coordinate_t>    ACoordinate_t;
  typedef K_FLD::ArrayAccessor<Connectivity_t>  AConnectivity_t;
  typedef K_SEARCH::BoundingBox<DIM>            BBox_t;
  typedef Vector_t<BBox_t*>                     BoxVector_t;
  typedef K_SEARCH::BbTree<DIM>                 Tree_t;
  typedef Predicate_t<Coordinate_t, Connectivity_t> Pred_t; 
    
  /// For FldArrays
  Collider(const K_FLD::FldArrayF& coords,  E_Int px, E_Int py, E_Int pz,
             const K_FLD::FldArrayI& connect, E_Float tolerance = EPSILON)
            :_aCoords(coords, px, py, pz), _aConnect(0),_connect(0),
             _predicate(new Pred_t),  _tree(0), _tolerance(tolerance), _own_tree(true){__init(connect, -1);}
  
  Collider(const K_FLD::ArrayAccessor<K_FLD::FldArrayF>& aCoord, K_FLD::ArrayAccessor<K_FLD::FldArrayI>& aConnect,
           Tree_t* tree, E_Float tolerance = EPSILON)
            :_aCoords(aCoord), _aConnect(&aConnect), _connect(0),
             _predicate(new Pred_t),  _tree(tree), _tolerance(tolerance), _own_tree(false)
  {
    _predicate->setLeft(&_aCoords, _aConnect);  
  }
      
  /// For DynArrays
  Collider(const K_FLD::FloatArray& coords, const K_FLD::IntArray& connect, E_Float tolerance =EPSILON)
            :_aCoords(coords), _aConnect(0),_connect(0),
             _predicate(new Pred_t), _tree(0), _tolerance(tolerance){__init(connect);}
   
  ///
  ~Collider(){if (_own_tree){__destroy_tree(); delete _connect; delete _aConnect;} delete _predicate;}
  
  ///
  inline const ACoordinate_t& getCoordinates() const { return _aCoords;}
  
  /// Get the colliding ids (FldArrays) with specified predicate.
  E_Int compute (const K_FLD::FldArrayF& coords,  E_Int px, E_Int py, E_Int pz,
                const K_FLD::FldArrayI& connect, Vector_t<E_Int>& colliding_ids);
  
  /// Get the colliding ids (DynArrays) with specified predicate.
  E_Int compute (const K_FLD::FloatArray& coords, const K_FLD::IntArray& connect,
                Vector_t<E_Int>& colliding_ids);
  
private:
  ///
  void __compute (const ACoordinate_t& coords, const AConnectivity_t& connect,
                  Vector_t<E_Int>& colliding_ids);
  ///
  void __init(const Connectivity_t& connect, E_Int shift = 0);
  
  ///
  void __create_tree(const ACoordinate_t& coords, const AConnectivity_t& connect, BoxVector_t& boxes, Tree_t*& tree);
  ///
  void __destroy_tree();
  
  ///
  bool __reorient(K_FLD::ArrayAccessor<Connectivity_t>*& connT4);
  
protected:
  ACoordinate_t    _aCoords;
  AConnectivity_t* _aConnect;
  Connectivity_t*  _connect;
  Vector_t<E_Int>  _colors;
  Pred_t*          _predicate;
  BoxVector_t      _boxes;
  Tree_t          *_tree;
  E_Float          _tolerance;
  bool             _own_tree;
  
};

template <typename Element_t>
class FldT3T3Collider3D : public Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, T3T3_XPredicate >
{
  public:
    typedef Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, T3T3_XPredicate > parent_t;
    ///
    FldT3T3Collider3D(const K_FLD::FldArrayF& coords,  E_Int px, E_Int py, E_Int pz,
                      const K_FLD::FldArrayI& connect, E_Float tolerance = EPSILON)
            :Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, T3T3_XPredicate >(coords, px, py, pz, connect, tolerance)
             {}
    ///
    virtual ~FldT3T3Collider3D(){}
    
};

template <typename Element_t>
class FldTH4T3Collider3D : public Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, TH4T3_XPredicate >
{
  public:
    typedef Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, TH4T3_XPredicate > parent_t;
    ///
    FldTH4T3Collider3D(const K_FLD::FldArrayF& coords,  E_Int px, E_Int py, E_Int pz,
                      const K_FLD::FldArrayI& connect, E_Float tolerance = EPSILON, typename parent_t::Tree_t* tree = 0)
            :Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, TH4T3_XPredicate >(coords, px, py, pz, connect, tolerance, tree)
             {}
    ///
    virtual ~FldTH4T3Collider3D(){}
    
};

template <typename Element_t>
class FldTH4HX6Collider3D : public Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, TH4HX6_XPredicate >
{
  public:
    typedef Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, TH4HX6_XPredicate > parent_t;
    ///
    FldTH4HX6Collider3D(const K_FLD::FldArrayF& coords,  E_Int px, E_Int py, E_Int pz,
                      const K_FLD::FldArrayI& connect, E_Float tolerance = EPSILON)
            :Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, TH4HX6_XPredicate >(coords, px, py, pz, connect, tolerance)
             {}
    FldTH4HX6Collider3D(const K_FLD::ArrayAccessor<K_FLD::FldArrayF>& coords,
                      K_FLD::ArrayAccessor<K_FLD::FldArrayI>& connect, typename parent_t::Tree_t* tree, E_Float tolerance = EPSILON)
            :Collider<K_FLD::FldArrayF, K_FLD::FldArrayI, 3, Element_t, TH4HX6_XPredicate >(coords, connect, tree, tolerance)
             {}
    ///
    virtual ~FldTH4HX6Collider3D(){}
    
};

template <typename Element_t>
class DynT3T3Collider3D : public Collider<K_FLD::FloatArray, K_FLD::IntArray, 3, Element_t, T3T3_XPredicate >
{
  public:
    typedef Collider<K_FLD::FloatArray, K_FLD::IntArray, 3, Element_t, T3T3_XPredicate > parent_t;
    ///
    DynT3T3Collider3D(const K_FLD::FloatArray& coords, const K_FLD::IntArray& connect, E_Float tolerance =EPSILON)
            : Collider<K_FLD::FloatArray, K_FLD::IntArray, 3, Element_t, T3T3_XPredicate >(coords, connect, tolerance)
            {}
    ///
    virtual ~DynT3T3Collider3D(){}
};


TEMPLATE_COORD_CONNECT_DIM_ELT_PRED
E_Int COLLIDER::compute
(const K_FLD::FldArrayF& coords,  E_Int px, E_Int py, E_Int pz,
 const K_FLD::FldArrayI& connect, Vector_t<E_Int>& is_colliding)
{
  ACoordinate_t aco(coords, px, py, pz);
  AConnectivity_t acv(connect, -1);
  _predicate->setRight(&aco, &acv);
  
  __compute(aco, acv, is_colliding);
  
#ifdef COLLIDE_DBG
  E_Int counter(0);
  for (size_t i=0; i <is_colliding.size(); ++i)
  {
    if (is_colliding[i])++counter;
  } 
  std::cout << "number of masked : " << counter << std::endl << std::endl;
#endif

  return 0;
}

TEMPLATE_COORD_CONNECT_DIM_ELT_PRED
E_Int COLLIDER::compute
(const K_FLD::FloatArray& coords, const K_FLD::IntArray& connect,
 Vector_t<E_Int>& is_colliding)
{
  is_colliding.resize(connect.cols(), 0);
  
  ACoordinate_t aco(coords);
  AConnectivity_t acv(connect);
  _predicate->setRight(&aco, &acv);
  
  __compute(aco, acv, is_colliding);

  return 0;
}

TEMPLATE_COORD_CONNECT_DIM_ELT_PRED
void COLLIDER::__compute
(const ACoordinate_t& coords, const AConnectivity_t& connect,
 Vector_t<E_Int>& is_colliding)
{
  K_FLD::IntArray e;
  E_Int s;
  size_t sz((size_t)connect.size());
  
  is_colliding.resize(sz, 0);

  if (_predicate)
  {
    Vector_t<E_Int> tmp;
    size_t i,j;
    Pred_t pred/*(*_predicate)*/;

#pragma omp parallel for private(i, j,e,s,tmp, pred) 
    for (i = 0; i < sz; ++i)
    {
      pred.set(*_predicate); //hack until correct optimization : _tmpT3s etc.. predicate attributes are guilty, somehow shared instead of being privates
      
      s = connect.stride(i);
      e.reserve(1, s);
      connect.getEntry(i, e.begin());
      K_SEARCH::BoundingBox<DIM> bb(coords, e.begin(), s);
    
      tmp.clear();
      _tree->getOverlappingBoxes(bb.minB, bb.maxB, tmp);
      
      /*if (i == 3574)
      {
        std::cout << " element : " << i << std::endl;
        std::cout << " stride : " << s << std::endl;
        std::cout << " box min : " << bb.minB[0] << " " <<bb.minB[1]<< " " << bb.minB[2]<< std::endl;
        std::cout << " box max : " << bb.maxB[0] << " " <<bb.maxB[1]<< " " << bb.maxB[2]<< std::endl;
        std::cout << "nb of caught boxes : " << tmp.size() << std::endl;
        
        //(pred)(6861, i) ;// N1:tmp[j], N2:i 
        
        K_FLD::IntArray out, candidates;
        K_FLD::FloatArray toto(coords.array());
        out.reserve(3, 12+tmp.size());
        out.resize(3,12);
        K_MESH::Hexahedron E;
        E.set(connect, i);
        E.triangulate(out.begin());
        Vector_t<E_Int> nids;
        //DELAUNAY::MeshTool::compact_to_mesh(toto, out, nids);
        
        E_Int T[3];
        for (size_t u=0; u< tmp.size(); ++u)
        {
          _aConnect->getEntry(tmp[u], T);
          //T[0] +=1; T[1]+=1; T[2]+=1;
          candidates.pushBack(T, T+3);
        }
        
        K_FLD::FloatArray tutu(_aCoords.array());
        meshIO::write("/home/slandier/tmp/cands.mesh", tutu, candidates);
           
        candidates.shift(toto.cols());
        toto.pushBack(tutu);
        out.pushBack(candidates);
        
        toto.pushBack(bb.minB, bb.minB+3);
        toto.pushBack(bb.maxB, bb.maxB+3);
        
        T[0]=T[1]=toto.cols()-2;T[2]=toto.cols()-1;
        out.pushBack(T,T+3);
        
        meshIO::write("/home/slandier/tmp/H6candiadtes.mesh", toto, out);
      }*/
    
      for (j=0; j<tmp.size(); ++j)
      {
        if ((pred)(tmp[j], i)) // N1:tmp[j], N2:i 
        {
          is_colliding[i]=1;
          /*if (i == 3574)
            std::cout << "colliding triangle : " << tmp[j] << std::endl;*/
          break;
        }
      }
    }
  }
  else
  {
    Vector_t<E_Int> tmp;
    size_t sz = connect.size();

#pragma omp parallel for private(e,s,tmp)
    for (size_t i= 0; i < sz; ++i)
    {
      s = connect.stride(i);
      e.reserve(1, s);
      connect.getEntry(i, e.begin());
      K_SEARCH::BoundingBox<DIM> bb(coords, e.begin(), s);

      if (_tree->hasAnOverlappingBox(bb.minB, bb.maxB))
        is_colliding[i]=1;
    }
  }
}

TEMPLATE_COORD_CONNECT_DIM_ELT_PRED
void COLLIDER::__init(const Connectivity_t& connect, E_Int shift)
{
  _aConnect = new AConnectivity_t(connect, shift);
  __reorient(_aConnect);
  __create_tree(_aCoords, *_aConnect, _boxes, _tree);
  _predicate->setLeft(&_aCoords, _aConnect);  
}

///
TEMPLATE_COORD_CONNECT_DIM_ELT_PRED
void COLLIDER::__create_tree(const ACoordinate_t& coords, const AConnectivity_t& connect, BoxVector_t& boxes, Tree_t*& tree)
{

#ifdef COLLIDE_DBG
  DELAUNAY::chrono c;
  c.start();
#endif
  
  K_FLD::IntArray e;
  E_Int s, sz(connect.size());
 
#ifdef COLLIDE_DBG 
  std::cout << "create " << sz << " boxes.." <<  std::endl;
#endif

  for (E_Int i = 0; i < sz; ++i)
  {
    s = connect.stride(i);
    e.reserve(1, s);
    connect.getEntry(i, e.begin());
    
    K_SEARCH::BoundingBox<DIM>* bb = new K_SEARCH::BoundingBox<DIM>(coords, e.begin(), s);
    boxes.push_back(bb);
  }
  
#ifdef COLLIDE_DBG
  std::cout << "boxes done : " << c.elapsed() << std::endl;
  
  std::cout << "create tree.." <<  std::endl;
  c.start();
#endif
  
  tree = new K_SEARCH::BbTree3D(boxes, _tolerance);
  
#ifdef COLLIDE_DBG
  std::cout << "tree done : " << c.elapsed() << std::endl;
#endif
}

///
TEMPLATE_COORD_CONNECT_DIM_ELT_PRED
void COLLIDER::__destroy_tree()
{
  for (size_t i = 0; i < _boxes.size();++i)delete _boxes[i];
  _boxes.clear();
  
  delete _tree;
  _tree=0;
}

///
  TEMPLATE_COORD_CONNECT_DIM_ELT_PRED
  bool COLLIDER::__reorient(K_FLD::ArrayAccessor<Connectivity_t>*& connT4)
  {
    Vector_t<size_t> elts_to_reorient;
    size_t nb_elts(connT4->size());
    const E_Int NB_NODES = connT4->stride();
    E_Float Ni[3], Nj[3], Nk[3], Nl[3];
    E_Int t4[NB_NODES];
    
    for (size_t i = 0; i < nb_elts; ++i)
    {
      connT4->getEntry(i, t4);
      
      _aCoords.getEntry(t4[0], Ni);
      _aCoords.getEntry(t4[1], Nj);
      _aCoords.getEntry(t4[2], Nk);
      _aCoords.getEntry(t4[3], Nl);
      
      if (NUGA::zzdet4(Ni, Nj, Nk, Nl) < -EPSILON) //must be reoriented
        elts_to_reorient.push_back(i);
    }
    
    size_t s = elts_to_reorient.size();
    
    if (s == 0)
      return false;
    //need to create a copy of the input connectivity
    Connectivity_t* new_connect = new Connectivity_t(connT4->array());
    K_FLD::ArrayWriter<Connectivity_t> writer(*connT4, *new_connect);
    //permut the first and second nodes
    for (size_t i = 0; i < s; ++i)
      std::swap(writer.getVal(elts_to_reorient[i], 0), writer.getVal(elts_to_reorient[i], 1));
    
    //reassign the new array to the accessor
    E_Int shift = connT4->shift();
    delete connT4;
    connT4 = new K_FLD::ArrayAccessor<Connectivity_t>(*new_connect, shift);
    
    return true;
  }

}

#endif
