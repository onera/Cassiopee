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

#ifndef NUGA_NODAL_SENSOR_HXX
#define NUGA_NODAL_SENSOR_HXX

#include "Nuga/include/sensor.hxx"


namespace NUGA
{

template <typename mesh_t>
class nodal_sensor : public sensor<mesh_t, std::vector<typename mesh_t::cell_incr_t>> // Vector_t might be templatized to have ISO mode with a metric field
{
  public:
  
    using cell_incr_t = typename mesh_t::cell_incr_t;
    using sensor_input_t = std::vector<cell_incr_t>;
    using parent_t = sensor<mesh_t, sensor_input_t>;
    using output_t = typename mesh_t::output_t; //fixme: static assert to add : must be ISO => IntVec

    nodal_sensor(mesh_t& mesh): parent_t(mesh, new V1_smoother<mesh_t>()){}

    virtual E_Int assign_data(const sensor_input_t& data) override;
    
    bool fill_adap_incr(output_t& adap_incr, bool do_agglo) override;
    bool update() override;
};

/// 
template <typename mesh_t>
E_Int nodal_sensor<mesh_t>::assign_data(const sensor_input_t& data)
{
  E_Int ncrd = parent_t::_hmesh._crd.cols();

  if (data.size() <= parent_t::_hmesh.pthids0.size() && !parent_t::_hmesh.pthids0.empty())
  {
    //converts input data to hmesh data
    sensor_input_t hmdat;
    hmdat.resize(ncrd, cell_incr_t(0));

    for (size_t k=0; k < data.size(); ++k){
      //std::cout << "id : " << k << " --> hmid : " << parent_t::_hmesh.pthids0[k] << " over " << ncrd << "points" << std::endl;
      hmdat[parent_t::_hmesh.pthids0[k]] = data[k];
    }

    parent_t::assign_data(hmdat);

  }
  else
  {
    parent_t::assign_data(data);
    parent_t::_data.resize(ncrd, cell_incr_t(0)); // resize anyway to ensure same size as crd
  }

  return 0;
}


///
template <typename mesh_t>
bool nodal_sensor<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  //
  bool filled{ false };
  sensor_input_t& Ln = parent_t::_data;
  
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  E_Int nb_elt= parent_t::_hmesh._ng.PHs.size();
  bool flag(false);

  adap_incr.cell_adap_incr.clear();
  adap_incr.face_adap_incr.clear();

  using cell_incr_t = typename output_t::cell_incr_t;
  using face_incr_t = typename output_t::face_incr_t;

  adap_incr.cell_adap_incr.resize(nb_elt, cell_incr_t(0));
  adap_incr.face_adap_incr.resize(nb_faces, face_incr_t(0));

  if (Ln.empty()) return false;

  for (int i=0; i< nb_elt; i++){
    if (parent_t::_hmesh._PHtree.is_enabled(i)){
      const E_Int* faces= parent_t::_hmesh._ng.PHs.get_facets_ptr(i);
      E_Int n_faces= parent_t::_hmesh._ng.PHs.stride(i);
      for (int k=0; k< n_faces; k++){
        E_Int PGk= *(faces+k)-1;
        const E_Int* pN= parent_t::_hmesh._ng.PGs.get_facets_ptr(PGk);
        E_Int n_face_nodes = parent_t::_hmesh._ng.PGs.stride(PGk);
        
        for (int l=0; l< n_face_nodes; l++){
          E_Int nodes_faces= *(pN+l)-1;
          if (Ln[nodes_faces]>0){
            adap_incr.cell_adap_incr[i]= 1;
            flag=true;
            filled = true;
            break;
          }
        }
        if (flag==true) break;
      }
    }
    flag=false;
  }
  return filled;
}

///
template <typename mesh_t>
bool nodal_sensor<mesh_t>::update()
{
  sensor_input_t& Ln = parent_t::_data;
  E_Int nb_pts = parent_t::_hmesh._crd.size();
  E_Int nb_new_pts = nb_pts - Ln.size();
  Ln.resize(nb_pts, cell_incr_t(IDX_NONE));

  bool updated{ false };

  //std::cout << "updating..." << std::endl;
  
  for (int i=0; i< nb_pts; i++)
  {
    if (Ln[i] > 0 && Ln[i] != IDX_NONE)
    {
      --Ln[i];
      updated = true;
    }
  }

  // interpolate new points 
  if (updated)
  {
    std::vector<E_Int> weight;
    E_Int nb_iter = nb_new_pts; //to avoid eventual infinite loop
    //E_Int count{0};
    while (nb_new_pts != 0 && nb_iter-- > 0)
    {
      weight.clear();
      weight.resize(nb_pts, 0);

      E_Int nb_new_pts0 = nb_new_pts;
      const ngon_unit& pgs = parent_t::_hmesh._ng.PGs;
      for (E_Int i = 0; i < pgs.size(); ++i)
      {
        E_Int nnodes = pgs.stride(i);
        const E_Int* pnodes = pgs.get_facets_ptr(i);

        for (E_Int j = 0; j < nnodes; ++j)
        {
          E_Int Nm1 = pnodes[j] - 1;
          E_Int N = pnodes[(j + 1) % nnodes] - 1;
          E_Int Np1 = pnodes[(j + 2) % nnodes] - 1;

          if (Ln[N] == IDX_NONE || weight[N] != 0) continue; // interpolated value (at this iter) does not contribute
          
          if (Ln[Nm1] == IDX_NONE || weight[Nm1] != 0) // to interpolate
          {
            if (Ln[Nm1] == IDX_NONE)
            {
              --nb_new_pts; //first update => reduce counter
              Ln[Nm1] = 0;
            }
            Ln[Nm1] += Ln[N];
            ++weight[Nm1];
          }

          if (Ln[Np1] == IDX_NONE || weight[Np1] != 0) // to interpolate
          {
            if (Ln[Np1] == IDX_NONE)
            {
              --nb_new_pts; //first update => reduce counter
              Ln[Np1] = 0;
            }
            Ln[Np1] += Ln[N];
            ++weight[Np1];
          }
        }
      }

      for (size_t i = 0; i < Ln.size(); ++i)
      {
        if (weight[i] != 0) Ln[i] /= weight[i];
      }

      if (nb_new_pts0 == nb_new_pts) break;

      //std::cout << "iter : " << count++ << " nb of processed : " << nb_new_pts0 - nb_new_pts<< std::endl;
    }

    //reset eventual non processed values
    for (size_t i = 0; i < Ln.size(); ++i)
    {
      if (Ln[i] == IDX_NONE) Ln[i] = 0;
    }

    //E_Int minv = *std::min_element(ALL(Ln));
    //E_Int maxv = *std::max_element(ALL(Ln));

    //std::cout << "VALS are bewteen " << minv << " and " << maxv << std::endl;
  }

  //std::cout << "updating done" << std::endl;

  return updated;
}

}
#endif
