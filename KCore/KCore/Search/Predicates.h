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

#ifndef __PREDICATES_H__
#define	__PREDICATES_H__

#include "Fld/ArrayAccessor.h"
#include "MeshElement/Tetrahedron.h"
#include "MeshElement/Hexahedron.h"

#ifdef DEBUG_COLLIDE_PRED
#include "meshIO.h"
inline void fooX(const E_Float* P0, const E_Float* P1, const E_Float* P2,
                 const E_Float* Q0, const E_Float* Q1, const E_Float* Q2)
{
  K_FLD::FloatArray pos;
  pos.pushBack(P0,P0+3);pos.pushBack(P1,P1+3);pos.pushBack(P2,P2+3);
  pos.pushBack(Q0,Q0+3);pos.pushBack(Q1,Q1+3);pos.pushBack(Q2,Q2+3);
  
  K_FLD::IntArray crd(3,2);
  crd(0,0)=0; crd(1,0)=1;crd(2,0)=2;
  crd(0,1)=3;crd(1,1)=4;crd(2,1)=5;
  
  meshIO::write("/home/slandier/tmp/x.mesh", pos, crd);
}
#endif

#define TEMPLATE_COORD_CONNECT template <typename Coordinate_t, typename Connectivity_t>

TEMPLATE_COORD_CONNECT
struct TruePredicate : public std::binary_function <E_Int, E_Int, bool>
{inline bool operator() (E_Int i, E_Int j) const {return true;}};


TEMPLATE_COORD_CONNECT
struct T3T3_XPredicate
{
  typedef K_FLD::ArrayAccessor<Coordinate_t>   ACoordinate_t;
  typedef K_FLD::ArrayAccessor<Connectivity_t> AConnectivity_t;
  
  //typedef K_MESH::Triangle Left_t;
  //typedef K_MESH::Triangle Right_t;
  
  T3T3_XPredicate(){}
  
  void setLeft(const ACoordinate_t* coords1, const AConnectivity_t* connect1){_coords1 = coords1; _connect1 = connect1;}
  void setRight(const ACoordinate_t* coords2, const AConnectivity_t* connect2){_coords2=coords2; _connect2=connect2;}
  
  void set(const T3T3_XPredicate& rhs){_coords1 = rhs._coords1; _connect1 = rhs._connect1; _coords2 = rhs._coords2; _connect2 = rhs._connect2;}
  
  inline bool operator() (E_Int N1, E_Int N2) const  
  {
    E_Int t1[3], t2[3];
    const E_Int DIM = 3;
    E_Float P0[DIM], P1[DIM], P2[DIM], Q0[DIM], Q1[DIM], Q2[DIM];
    
    _connect1->getEntry(N1, t1);    
    _coords1->getEntry(t1[0], P0);
    _coords1->getEntry(t1[1], P1);
    _coords1->getEntry(t1[2], P2);
    
    _connect2->getEntry(N2, t2);
    _coords2->getEntry(t2[0], Q0);
    _coords2->getEntry(t2[1], Q1);
    _coords2->getEntry(t2[2], Q2);
    
    return K_MESH::Triangle::fast_intersectT3<DIM>(P0, P1, P2, Q0, Q1, Q2, -1.);
  }
          
  const ACoordinate_t* _coords1;
  const ACoordinate_t*_coords2;
  const AConnectivity_t* _connect1;
  const AConnectivity_t* _connect2;
  
};

TEMPLATE_COORD_CONNECT
struct TH4T3_XPredicate
{
  typedef K_FLD::ArrayAccessor<Coordinate_t>   ACoordinate_t;
  typedef K_FLD::ArrayAccessor<Connectivity_t> AConnectivity_t;
  
  //typedef K_MESH::Triangle Left_t;
  //typedef K_MESH::Triangle Right_t;
  
  TH4T3_XPredicate(){_tmpT3s.resize(3,4);}
  
  TH4T3_XPredicate(const TH4T3_XPredicate& rhs):_coords1(rhs._coords1), _coords2(rhs._coords2), _connect1(rhs._connect1), _connect2(rhs._connect2)
  {
  }
  
  TH4T3_XPredicate& operator=(const TH4T3_XPredicate& rhs)
  {
    _coords1=rhs._coords1;
    _coords2=rhs._coords2;
    _connect1=rhs._connect1;
    _connect2=rhs._connect2;
    return *this;
  }
  
  
  void setLeft(const ACoordinate_t* coords1, const AConnectivity_t* connect1){_coords1 = coords1; _connect1 = connect1;}
  void setRight(const ACoordinate_t* coords2, const AConnectivity_t* connect2){_coords2=coords2; _connect2=connect2;}
  
  
  inline bool operator() (E_Int N1, E_Int N2) const  
  {
    E_Int t1[4], t2[3];
    const E_Int DIM = 3;
    E_Float P0[DIM], P1[DIM], P2[DIM], P3[DIM], Q0[DIM], Q1[DIM], Q2[DIM];
    
    _connect1->getEntry(N1, t1);
    _connect2->getEntry(N2, t2);
    
    _coords1->getEntry(t1[0], P0);
    _coords1->getEntry(t1[1], P1);
    _coords1->getEntry(t1[2], P2);
    _coords1->getEntry(t1[3], P3);
    
    //fast check one of the triangle nodes is inside
    
    _coords2->getEntry(t2[0], Q0);
    if (K_MESH::Tetrahedron::is_inside(P0, P1, P2, P3, Q0))
      return true;
    
    _coords2->getEntry(t2[1], Q1);
    if (K_MESH::Tetrahedron::is_inside(P0, P1, P2, P3, Q1))
      return true;
    
    _coords2->getEntry(t2[2], Q2);
    if (K_MESH::Tetrahedron::is_inside(P0, P1, P2, P3, Q2))
      return true;
    
    _TH4i.setNodes(&t1[0]);
    _TH4i.triangulate(_tmpT3s.begin());
    
    K_FLD::IntArray::const_iterator pK;
    for (size_t i=0; i < 4; ++i)
    {
      pK = _tmpT3s.col(i);
      _coords1->getEntry(*pK, P0);
      _coords1->getEntry(*(pK+1), P1);
      _coords1->getEntry(*(pK+2), P2);
      
      if (K_MESH::Triangle::fast_intersectT3<DIM>(P0, P1, P2, Q0, Q1, Q2, -1.))
        return true;
    }
    return false;
  }
          
  const ACoordinate_t* _coords1;
  const ACoordinate_t*_coords2;
  const AConnectivity_t* _connect1;
  const AConnectivity_t* _connect2;
  
  mutable K_FLD::IntArray _tmpT3s;
  mutable K_MESH::Tetrahedron _TH4i;
  
};

TEMPLATE_COORD_CONNECT
struct TH4HX6_XPredicate
{
  typedef K_FLD::ArrayAccessor<Coordinate_t>   ACoordinate_t;
  typedef K_FLD::ArrayAccessor<Connectivity_t> AConnectivity_t;
  
  //typedef K_MESH::Triangle Left_t;
  //typedef K_MESH::Triangle Right_t;
  
  TH4HX6_XPredicate(){_tmpT3s.resize(3,4);}
  
  TH4HX6_XPredicate(const TH4HX6_XPredicate& rhs):_coords1(rhs._coords1), _coords2(rhs._coords2), _connect1(rhs._connect1), _connect2(rhs._connect2)
  {
  }
  
  TH4HX6_XPredicate& operator=(const TH4HX6_XPredicate& rhs)
  {
    _coords1=rhs._coords1;
    _coords2=rhs._coords2;
    _connect1=rhs._connect1;
    _connect2=rhs._connect2;
    return *this;
  }
  
  void setLeft(const ACoordinate_t* coords1, const AConnectivity_t* connect1){_coords1 = coords1; _connect1 = connect1;}
  void setRight(const ACoordinate_t* coords2, const AConnectivity_t* connect2){_coords2=coords2; _connect2=connect2;}
  
  void set(const TH4HX6_XPredicate& rhs){_coords1 = rhs._coords1; _connect1 = rhs._connect1; _coords2 = rhs._coords2; _connect2 = rhs._connect2;}
  
  
  inline bool operator() (E_Int N1, E_Int N2) const  
  {
    E_Int t1[4], t2[8];
    const E_Int DIM = 3;
    E_Float P0[DIM], P1[DIM], P2[DIM], P3[DIM], Q0[DIM], Q1[DIM], Q2[DIM];
    
    _connect1->getEntry(N1, t1);
    _connect2->getEntry(N2, t2);
    
    _coords1->getEntry(t1[0], P0);
    _coords1->getEntry(t1[1], P1);
    _coords1->getEntry(t1[2], P2);
    _coords1->getEntry(t1[3], P3);
    
    //fast check one of the hexa nodes is inside
    for (size_t i=0; i < 8; ++i)
    {
      _coords2->getEntry(t2[i], Q0);
      if (K_MESH::Tetrahedron::is_inside(P0, P1, P2, P3, Q0))
       return true;
    }
    
    _TH4i.setNodes(&t1[0]);
    _TH4i.triangulate(_tmpT3s.begin());
    
    _HX6i.setNodes(t2);
    
    K_FLD::IntArray::const_iterator pK;
    E_Int q4[4];
    for (size_t f=0; f <6; ++f)
    {
      _HX6i.getBoundary(f, q4);
      
      _coords2->getEntry(q4[0], Q0);
      _coords2->getEntry(q4[1], Q1);
      _coords2->getEntry(q4[2], Q2);
       
      for (size_t i=0; i < 4; ++i)
      {
        pK = _tmpT3s.col(i);
       _coords1->getEntry(*pK, P0);
       _coords1->getEntry(*(pK+1), P1);
       _coords1->getEntry(*(pK+2), P2);
      
       if (K_MESH::Triangle::fast_intersectT3<DIM>(P0, P1, P2, Q0, Q1, Q2, -1.))
         return true;
       
       _coords2->getEntry(q4[3], Q1);
       
       if (K_MESH::Triangle::fast_intersectT3<DIM>(P0, P1, P2, Q0, Q2, Q1, -1.))
         return true;
      }
    }
    
    return false;
  }
          
  const ACoordinate_t* _coords1;
  const ACoordinate_t*_coords2;
  const AConnectivity_t* _connect1;
  const AConnectivity_t* _connect2;
  
  mutable K_FLD::IntArray _tmpT3s;
  mutable K_MESH::Tetrahedron _TH4i;
  mutable  K_MESH::Hexahedron _HX6i;
  
};

TEMPLATE_COORD_CONNECT
struct T3HX6_XPredicate
{
  typedef K_FLD::ArrayAccessor<Coordinate_t>   ACoordinate_t;
  typedef K_FLD::ArrayAccessor<Connectivity_t> AConnectivity_t;
  
  //typedef K_MESH::Triangle Left_t;
  //typedef K_MESH::Triangle Right_t;
  
  T3HX6_XPredicate(){;}
  
  T3HX6_XPredicate(const T3HX6_XPredicate& rhs):_coords1(rhs._coords1), _coords2(rhs._coords2), _connect1(rhs._connect1), _connect2(rhs._connect2)
  {
  }
  
  T3HX6_XPredicate& operator=(const T3HX6_XPredicate& rhs)
  {
    _coords1=rhs._coords1;
    _coords2=rhs._coords2;
    _connect1=rhs._connect1;
    _connect2=rhs._connect2;
    return *this;
  }
  
  void setLeft(const ACoordinate_t* coords1, const AConnectivity_t* connect1){_coords1 = coords1; _connect1 = connect1;}
  void setRight(const ACoordinate_t* coords2, const AConnectivity_t* connect2){_coords2=coords2; _connect2=connect2;}
  
  void set(const T3HX6_XPredicate& rhs){_coords1 = rhs._coords1; _connect1 = rhs._connect1; _coords2 = rhs._coords2; _connect2 = rhs._connect2;}
  
  
  inline bool operator() (E_Int N1, E_Int N2) const  
  {
    E_Int t1[3], t2[8];
    const E_Int DIM = 3;
    E_Float P0[DIM], P1[DIM], P2[DIM], Q0[DIM], Q1[DIM], Q2[DIM];
    
 #ifdef DEBUG_COLLIDE_PRED
    bool enable=(N1==6861 && N2==3574);
#endif
    
    _connect1->getEntry(N1, t1);
    _connect2->getEntry(N2, t2);
    
    _coords1->getEntry(t1[0], P0);
    _coords1->getEntry(t1[1], P1);
    _coords1->getEntry(t1[2], P2);    
        
    _HX6i.setNodes(t2);
    
    E_Int q4[4];
    for (size_t f=0; f <6; ++f)
    {
      _HX6i.getBoundary(f, q4);
      
      _coords2->getEntry(q4[0], Q0);
      _coords2->getEntry(q4[1], Q1);
      _coords2->getEntry(q4[2], Q2);
            
      if (K_MESH::Triangle::fast_intersectT3<DIM>(P0, P1, P2, Q0, Q1, Q2, -1.))
      {
#ifdef DEBUG_COLLIDE_PRED
        if( enable)
        {
          std::cout << "first : " << f << std::endl;
          fooX(P0, P1, P2, Q0, Q1, Q2);
        }
#endif
         return true;
      }
      
      _coords2->getEntry(q4[3], Q1);
      
      if (K_MESH::Triangle::fast_intersectT3<DIM>(P0, P1, P2, Q0, Q2, Q1, -1.))
      {
#ifdef DEBUG_COLLIDE_PRED
        if( enable)
        {
          std::cout << "second : " << f << std::endl;
          fooX(P0, P1, P2, Q0, Q1, Q2);
        }
#endif
         return true;
      }
    }
    
    return false;
  }
          
  const ACoordinate_t* _coords1;
  const ACoordinate_t*_coords2;
  const AConnectivity_t* _connect1;
  const AConnectivity_t* _connect2;
  
  mutable  K_MESH::Hexahedron _HX6i;
  
};


#endif	/* PREDICATES_H */
