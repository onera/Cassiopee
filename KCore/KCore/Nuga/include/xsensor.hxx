/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_XSENSOR_HXX
#define NUGA_XSENSOR_HXX

#include "Nuga/include/geom_sensor.hxx"
#include "Nuga/include/Hexahedron.h"


namespace NUGA
{

struct mesh_s
{
  K_FLD::FloatArray crd;
  K_FLD::IntArray  cnt;

  E_Int getSize() { return crd.cols(); }
  E_Float* get(E_Int i) { return crd.col(i); }
};

/// X geometric sensor
template <typename ELT_t, typename mesh_t> //ngu for surfacic (PGs) or ngon_t for volumic
class xsensor : public geom_sensor<mesh_t, mesh_s>
{  
  public:
    using elt_type = ELT_t;
    using parent_t = geom_sensor<mesh_t, mesh_s>;
    using sensor_input_t = mesh_s;
    
    xsensor(mesh_t& mesh, eSmoother smootyp, E_Int max_pts_per_cell, E_Int itermax):parent_t(mesh, smootyp, max_pts_per_cell, itermax){}

    E_Int assign_data(const sensor_input_t& data) ; //not an override as args are different

 private:
    void add_x_points(sensor_input_t& data);            
};
///
template <typename ELT_t, typename mesh_t>
void xsensor<ELT_t, mesh_t>::add_x_points(sensor_input_t& data)
{
  Vector_t<E_Int> ids;
  std::map<K_MESH::NO_Edge,E_Int> unique_edges;
  //
  for (int i = 0; i < data.cnt.cols(); ++i)// for each elts of the source file, check if there is a x situation
  {      
    // get the edges
    Vector_t<K_MESH::NO_Edge> edges(ELT_t::NB_EDGES);
    ELT_t::get_edges(data.cnt.col(i),edges);
    
    Vector_t<E_Float> lambdas;
    
    // are there edges crossing the NGON mesh ?
    for (int j = 0; j < 12; ++j) // for each edge
    {
      K_MESH::NO_Edge noE = edges[j];
      auto it = unique_edges.find(noE);
      if (it != unique_edges.end()) continue;
      unique_edges[noE] = IDX_NONE;
      E_Int pt1 = noE.node(0);
      E_Int pt2 = noE.node(1);

      E_Float* P0 = data.crd.col(pt1);
      E_Float* P1 = data.crd.col(pt2);

      ids.clear();
      parent_t::_bbtree->getOverlappingBoxes(P0,P1,ids);
      lambdas.clear();
      
      if (ids.empty()) continue;
      
      E_Float lambda0;
      E_Float lambda1;
      E_Float tolerance = 1.e-8;
      
      for (size_t k = 0; k < ids.size(); ++k) // for each element of hmesh that may intersect with noE
      {
        E_Int PHi = ids[k];
        const E_Int* face = parent_t::_hmesh._ng.PHs.get_facets_ptr(PHi);
        E_Int nb_faces = parent_t::_hmesh._ng.PHs.stride(PHi);

        // check if noE intersect with PHi
        bool x = ELT_t::cross(parent_t::_hmesh._ng, parent_t::_hmesh._crd, face, nb_faces, parent_t::_data.crd, P0, P1, lambda0, lambda1, tolerance);
        
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
        data.crd.pushBack(xpoint, xpoint+3);
      }
    }
  }
}

///
template <typename ELT_t, typename mesh_t>
E_Int xsensor<ELT_t, mesh_t>::assign_data(const sensor_input_t& data)
{
  parent_t::_data = data;
  add_x_points(parent_t::_data);

  parent_t::_iter = 0; //reset
             // fill in _points_to_cell with the initial mesh (giving for each source point the highest level cell containing it)
  parent_t::_points_to_cell.clear();

  // Count the number of data per cell

  // Check if hmesh has been initialized
  if (not parent_t::_hmesh.is_initialised())
  {
    std::cout << "Error in geom_sensor.hxx / function assign_data: hmesh has not been initialised." << std::endl;
    return 1;
  }

  parent_t::locate_points();
  
  return 0;
}


}


#endif

