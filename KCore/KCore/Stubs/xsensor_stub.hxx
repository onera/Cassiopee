/*
 
 
 
              NUGA 
 
 
 
 */
// Authors: Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#ifndef NUGA_XSENSOR_HXX
#define NUGA_XSENSOR_HXX

#include "Nuga/include/geom_sensor.hxx"
#include "Nuga/include/hierarchical_mesh.hxx"
#include "MeshElement/Hexahedron.h"


namespace NUGA
{

/// X geometric sensor
template <typename ELT_t, typename mesh_t, typename crd_t = K_FLD::FloatArray> //ngu for surfacic (PGs) or ngon_t for volumic
class xsensor : public geom_sensor<mesh_t, crd_t>

{  
  public:
    using elt_type = ELT_t;
    
    using data_type = crd_t; //point cloud
    
    xsensor(mesh_t& mesh, const K_FLD::IntArray& cntS, E_Int itermax = -1):geom_sensor<mesh_t,crd_t>(mesh, 1/*max_pts_per_cell*/, itermax), _cntS(cntS){}

    E_Int init(data_type& data) override;
    
    void add_x_points(data_type& data, K_SEARCH::BbTree3D& tree);
    
  private:
    const K_FLD::IntArray& _cntS;
            
};

///
template <>
void xsensor<K_MESH::Hexahedron, NUGA::hierarchical_mesh<K_MESH::Hexahedron>, K_FLD::FloatArray>::add_x_points(data_type& data, K_SEARCH::BbTree3D& tree)
{
  Vector_t<E_Int> ids;
  std::map<K_MESH::NO_Edge,E_Int> unique_edges;
  // add the x points (x <=> intersection)
  // cntS for hexa : 8 lines for nodes and one column per PHi
  E_Int nb_nodes = 8;
  for (int i = 0; i < _cntS.cols(); ++i)// for each elts of the source file, check if there is a x situation
  {
    E_Int nodes[8];
    // fill nodes with the 8 nodes of PHi
    for (int j = 0; j < nb_nodes; ++j)
        nodes[j] = _cntS.col(i)[j];
      
    // get the 12 edges
    Vector_t<K_MESH::NO_Edge> edges(12);
    edges.clear();
    K_MESH::Hexahedron::get_edges(nodes,edges);
    
    Vector_t<E_Float> lambdas;
    
    // are there edges crossing the NGON mesh ?
    for (int j = 0; j < 12; ++j) // for each edge
    {
      K_MESH::NO_Edge noE = edges[j];
      auto it = unique_edges.find(noE);
      if (it != unique_edges.end()) continue;
      unique_edges[noE] = E_IDX_NONE;
      E_Int pt1 = noE.node(0);
      E_Int pt2 = noE.node(1);

      E_Float* P0 = data.col(pt1);
      E_Float* P1 = data.col(pt2);

      ids.clear();
      tree.getOverlappingBoxes(P0,P1,ids);
      lambdas.clear();
      
      if (ids.empty()) continue;
      
      E_Float lambda0;
      E_Float lambda1;
      E_Float tolerance = 1.e-8;
      
      for (size_t k = 0; k < ids.size(); ++k) // for each element of hmesh that may intersect with noE
      {
        E_Int PHi = ids[k];
        E_Int* face = _hmesh._ng.PHs.get_facets_ptr(PHi);
        E_Int nb_faces = _hmesh._ng.PHs.stride(PHi);

        // check if noE intersect with PHi
        bool x = K_MESH::Hexahedron::cross(_hmesh._ng, _hmesh._crd, face, nb_faces, data, P0, P1, lambda0, lambda1, tolerance);
        
        if (!x) continue; // no intersection
        
        lambdas.push_back(lambda0);
        lambdas.push_back(lambda1);
      }
      
      if (lambdas.empty()) continue;
      // sort the lambdas
      std::sort(lambdas.begin(),lambdas.end());
      // ne pas ajouter deux points dont les lambdas ont un écart plus petit que tolérance
      // add the x points in data
      E_Float coord0[3];
      E_Float coord1[3];
      for (int p = 0; p < 3; ++p)
      {
        coord0[p] = P0[p];
        coord1[p] = P1[p];
      }
      for (size_t l = 0; l < lambdas.size()-1; ++l) // we ignore the 2 parts from a pt " edge to a lambda " (only add a point between 2 lambdas)
      {
        // O -/-/-/- x --N-- x --N-- x -/-/-/- O : N are the points created, x are the lambda points and 0 the 2 points which define the edge
        E_Float x1 = lambdas[l];
        E_Float x2 = lambdas[l+1];

        if (fabs(x1 - x2) < tolerance) continue;
        
        //add the points between the 2 lambdas in data
        E_Float x3 = (x1 + x2) /2;
        E_Float xpoint[3];
        for (int m = 0; m < 3; ++m)
          xpoint[m] = coord0[m] + x3*(coord1[m] - coord0[m]); // the new point
        data.pushBack(xpoint, xpoint+3);
      }
    }
  }
}

///
template <>
E_Int xsensor<K_MESH::Hexahedron, NUGA::hierarchical_mesh<K_MESH::Hexahedron>, K_FLD::FloatArray>::init(data_type& data)
{
  // fill in _points_to_cell with the initial mesh (giving for each source point the highest level cell containing it)
  _points_to_cell.clear();
  
  K_SEARCH::BbTree3D tree(_hmesh._crd, _hmesh._ng);

  // add the intersection points
  add_x_points(data, tree);
  
  // localize points : fill _points_to_cell
  geom_sensor<NUGA::hierarchical_mesh<K_MESH::Hexahedron>, K_FLD::FloatArray>::locate_points(tree, data);
  
  return 0;
}


}


#endif

