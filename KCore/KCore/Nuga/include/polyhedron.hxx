/*
 
 
 
              NUGA 
 
 
 
 */
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_POLYHEDRON_HXX
#define NUGA_POLYHEDRON_HXX

#include "MeshElement/Polyhedron.h"

namespace NUGA
{
///autonomous (i.e. self contained) polyhedron
template <int TopoShape>
struct aPolyhedron : public K_MESH::Polyhedron<TopoShape>
{
  using parent_type = K_MESH::Polyhedron<TopoShape>;
  
  ngon_unit          m_pgs;
  std::vector<E_Int> m_faces;
  K_FLD::FloatArray  m_crd;
  
  aPolyhedron(){plug(); parent_type::_triangles = nullptr;}
  aPolyhedron(const ngon_unit* pgs, const E_Int* faces, E_Int nb_faces, const K_FLD::FloatArray& crd); // from "mesh" to autonmous
  aPolyhedron(const parent_type& ph, const K_FLD::FloatArray& crd); // from "mesh" to autonmous
  aPolyhedron(const ngon_unit& lpgs, const K_FLD::FloatArray& crd);
 
  aPolyhedron(const aPolyhedron& r):m_pgs(r.pgs), m_faces(r.m_faces), m_crd(r.m_crd){plug();}
  
    
  aPolyhedron& operator=(const parent_type& rhs) = delete;
  aPolyhedron& operator=(const aPolyhedron& rhs) = delete;
          
  void set(const parent_type& ph, const K_FLD::FloatArray& crd);
  
  void reorient();
  
  void clear() {m_pgs.clear(); m_faces.clear(); m_crd.clear();}
  
  template<typename TriangulatorType> E_Int volume(E_Float& v, bool need_reorient);
  
  void plug(){parent_type::_pgs = &m_pgs; parent_type::_faces = &m_faces[0]; parent_type::_nb_faces = m_faces.size();}
  
};

///
template <int TopoShape>
aPolyhedron<TopoShape>::aPolyhedron(const ngon_unit* pgs, const E_Int* faces, E_Int nb_faces, const K_FLD::FloatArray& crd)
{
  parent_type ph(pgs, faces, nb_faces);
  set(ph, crd);
  
}

///
template <int TopoShape>
aPolyhedron<TopoShape>::aPolyhedron(const parent_type& ph, const K_FLD::FloatArray& crd)
{
  set(ph, crd);
}

///
template <int TopoShape>
aPolyhedron<TopoShape>::aPolyhedron(const ngon_unit& lpgs, const K_FLD::FloatArray& lcrd)
{
  m_pgs = lpgs;
  m_crd = lcrd;
  
  m_faces.clear();
  K_CONNECT::IdTool::init_inc(m_faces, m_pgs.size(), 1);
  
  plug();
}

///
template <int TopoShape>
void aPolyhedron<TopoShape>::set(const parent_type& ph, const K_FLD::FloatArray& crd)
{
  ph.compact(crd, m_pgs, m_crd);
  m_faces.clear();
  K_CONNECT::IdTool::init_inc(m_faces, m_pgs.size(), 1);
  
  plug();
}

///
template <int TopoShape>
void aPolyhedron<TopoShape>::reorient()
{
  ngon_unit neighbors;
  K_MESH::Polygon::build_pg_neighborhood(m_pgs, neighbors);
  std::vector<E_Int> oids, orient(m_faces.size(), 1);
  K_CONNECT::EltAlgo<K_MESH::Polygon>::reversi_connex(m_pgs, neighbors, 0/*ref*/, orient);

  //Apply new orientation
  for (E_Int i = 0; i < m_pgs.size(); ++i)
  {
    if (orient[i] == -1)
    {
      E_Int s = m_pgs.stride(i);
      E_Int* p = m_pgs.get_facets_ptr(i);
      std::reverse(p, p + s);
    }
  }
}

template <int TopoShape>
template<typename TriangulatorType>
E_Int aPolyhedron<TopoShape>::volume(E_Float& v, bool need_reorient)
{
  return parent_type::template volume<TriangulatorType>(m_crd, v, need_reorient);
}


// autonomous with history
template <int TopoShape>
struct haPolyhedron : public aPolyhedron<TopoShape>
{
  using parent_type = NUGA::aPolyhedron<TopoShape>;
  using base_type = K_MESH::Polyhedron<TopoShape>;
  
  std::vector<E_Int> poids;

  haPolyhedron():parent_type(){};//=default;
  
  haPolyhedron& operator=(const base_type& rhs) = delete; // crd must be also specified so disable this operation (permanent)
  haPolyhedron& operator=(const parent_type& rhs) = delete;// {set(rhs, rhs.m_crd);} delet until required
  haPolyhedron& operator=(const haPolyhedron& rhs);
          
  void set(const base_type& ph, const K_FLD::FloatArray& crd);
  void set(ngon_type& ng, E_Int PHi, const K_FLD::FloatArray& crd);
  
};

///
template <int TopoShape>
haPolyhedron<TopoShape>& haPolyhedron<TopoShape>::operator=(const haPolyhedron<TopoShape>& rhs)
{
  poids = rhs.poids;
  parent_type::m_pgs = rhs.m_pgs;
  parent_type::m_faces = rhs.m_faces;
  parent_type::m_crd = rhs.m_crd;
  
  parent_type::plug();
}

///
template <int TopoShape>
void haPolyhedron<TopoShape>::set(const base_type& ph, const K_FLD::FloatArray& crd)
{
  ph.compact(crd, parent_type::m_pgs, parent_type::m_crd, &poids);
  parent_type::m_faces.clear();
  K_CONNECT::IdTool::init_inc(parent_type::m_faces, parent_type::m_pgs.size(), 1);
  
  parent_type::plug();
}

///
template <int TopoShape>
void haPolyhedron<TopoShape>::set(ngon_type& ng, E_Int PHi, const K_FLD::FloatArray& crd)
{
  base_type e(ng,PHi);
  set(e, crd);
}

}

#endif
