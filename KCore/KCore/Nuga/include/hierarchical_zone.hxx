/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr), Alexis Rouil (alexis.rouil@onera.fr)

#ifndef NUGA_HIERACHICAL_ZONE_HXX
#define NUGA_HIERACHICAL_ZONE_HXX

#include "Nuga/include/hierarchical_mesh.hxx"
#include "Nuga/include/join_sensor.hxx"
#include "Nuga/include/communicator.hxx"
#include "Nuga/include/join_t.hxx"

namespace NUGA
{

  using crd_t = K_FLD::FloatArray;

  ///
  template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t = ngon_type>
  struct hierarchical_zone : public hierarchical_mesh<ELT_t, STYPE, ngo_t>
  {
    using hmesh_t = hierarchical_mesh<ELT_t, STYPE, ngo_t>;
    using jsensor_t = join_sensor<hmesh_t>;
    using join_data_t = std::map<E_Int, std::vector<E_Int>>;
    using jcom_t = jsensor_com_agent<hmesh_t, typename jsensor_t::input_t>;
    using communicator_t = NUGA::communicator<jcom_t>;

    hierarchical_zone(E_Int id, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, const join_data_t& jdata, E_Int idx_start, communicator_t& com);

    ~hierarchical_zone()
    { 
      if (join            != nullptr) delete join;
      if (jsensor         != nullptr) delete jsensor;
      if (COM.agents[zid] != nullptr) delete COM.agents[zid];
    }

    E_Int            zid;
    join_t<hmesh_t>* join;
    jsensor_t*       jsensor;
    communicator_t&  COM;

  };
  
  ///
  template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t = ngon_type>
  hierarchical_zone<ELT_t, STYPE, ngo_t>::hierarchical_zone
  (E_Int id, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, const join_data_t& jdata, E_Int idx_start, communicator_t& com) :
   hmesh_t(crd, cnt), zid(id), join(nullptr), jsensor(nullptr), COM(com)
  {
    join_data_t jmp = jdata;
    if (!jmp.empty()) // join is specified
    {
      jsensor = new jsensor_t(*this);
      join    = new join_t<hmesh_t>(id, *this, idx_start);

      assert(COM.agents.size() > zid); //COM agents must be resized by the caller before this call

      com.agents[zid] = new jcom_t(zid, join, jsensor);

      for (auto j : jmp)
      {
        E_Int joinedZid = j.first;
        std::vector<E_Int>& ptlist = j.second;

        join->link(joinedZid, ptlist);
        COM.agents[zid]->link(joinedZid);
      }
    }
  }
}

#endif
