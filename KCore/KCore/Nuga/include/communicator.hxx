/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_COMMUNICATOR_HXX
#define NUGA_COMMUNICATOR_HXX

#include <vector>
#include "Nuga/include/join_sensor.hxx"
#include "Nuga/include/join_t.hxx"
#include "Nuga/include/adaptor.hxx"

namespace NUGA
{
 
  ///
  template <typename ExData> // IntArray (or ngon_unit, but not impl)
  struct com_slots
  {
    ExData inslot;
    ExData outslot;
  };


  struct com_agent
  {
    E_Int id; // agent id

    com_agent(E_Int i) :id(i) {};
    virtual ~com_agent(){} //for msys2/gcc : complain even if no attributes to delete

    virtual void link(E_Int agent_id/*id to link with*/) = 0;
    virtual bool pack()   = 0;  // for sending
    virtual bool unpack() = 0;  // at reception
    virtual bool process() = 0; // do something with the received data
  };

  ///
  template <typename agent_t>
  struct communicator
  {
    std::vector<agent_t*> agents;

    void isend(E_Int sender_id);

  };

  ///
  template <typename agent_t>
  void communicator<agent_t>::isend(E_Int sender_id)
  {
    const agent_t& Sag = *agents[sender_id];
    if (Sag.jsensor == nullptr) return;

    //const auto& Smesh = Sag.jsensor->get_hmesh();

    for (auto& j : Sag.jid_to_boxes)
    {
      E_Int Rid = j.first;
      auto& outslot = j.second.outslot;           // departure slot

      agent_t& Rag = *agents[Rid];
      if (Rag.jsensor == nullptr) return;

      auto it1 = Rag.jid_to_boxes.find(sender_id);
      auto& inslot = it1->second.inslot;         // arrival slot

      auto it2 = Rag.join->jid_to_joinlist.find(sender_id);
      auto& Rptlist = it2->second;

      for (auto it : outslot)
      {
        E_Int PGi = Rptlist[it.first] - Rag.join->idx_start;
        inslot[PGi] = it.second;
      }
    }
  }

  ///
  template <typename mesh_t, typename ExData>
  struct jsensor_com_agent : public com_agent
  {
    using data_t = ExData;
    using join_sensor_t = join_sensor<mesh_t>;

    jsensor_com_agent(E_Int id, join_t<mesh_t>* j, join_sensor_t* js) :com_agent(id), join(j), jsensor(js) {}

    join_t<mesh_t>*                     join;
    join_sensor_t*                      jsensor;
    std::map<E_Int, com_slots<ExData>>  jid_to_boxes;

    void link(E_Int jzid) override;

    bool pack()   override; // for sending
    bool unpack() override; // at reception
    bool process() override; // do something with the received data

  };

  ///
  template <typename mesh_t, typename ExData>
  void jsensor_com_agent<mesh_t, ExData>::link(E_Int jzid)
  {
    jid_to_boxes[jzid]; //create empty slots
  }


  /// prepare packs to send to all joined zones, put them on each output slot
  template <typename mesh_t, typename ExData>
  bool jsensor_com_agent<mesh_t, ExData>::pack()
  {
    using pg_arr_t = typename mesh_t::pg_arr_t;
    //std::cout << "pack : begin" << std::endl;
    if (jsensor == nullptr) return false;
    //std::cout << "pack : get_hmesh" << std::endl;
    const mesh_t& mesh = jsensor->get_hmesh();
    bool has_packs{ false };
    //std::cout << "pack : loop on joins : " << join << std::endl;
    for (auto& j : join->jid_to_joinlist)
    {
      E_Int jid = j.first;
      auto& ptlist = j.second;

      //std::cout << "jid/ptl sz : " << jid << "/" << ptlist.size() << std::endl;

      jid_to_boxes[jid].outslot.clear();

      for (size_t i = 0; i < ptlist.size(); ++i)
      {
        pg_arr_t p; //only implemented for IntArray
        E_Int PGi = ptlist[i] - join->idx_start;
        //std::cout << "i/PGi : " << i << "/" << PGi << std::endl;
        mesh.extract_plan(PGi, true/*reverse*/, 0/*because previous sync*/, p);
        //std::cout << "after extract_plan" << std::endl;
        if (p.getSize())
        {
          jid_to_boxes[jid].outslot[i] = p; //assume same sorting of joinlists
          has_packs = true;
        }
      }
    }
    //std::cout << "pack : end" << std::endl;
    return has_packs;
  }

  ///
  template <typename mesh_t, typename ExData>
  bool jsensor_com_agent<mesh_t, ExData>::unpack()
  {
    if (jsensor == nullptr) return false;
    
    // gather inputs
    data_t data;
    for (auto& j : jid_to_boxes)
      data.insert(ALL(j.second.inslot));
 
    jsensor->assign_data(data); // pass it to sensor

    return !data.empty();
  }

  template <typename mesh_t, typename ExData>
  bool jsensor_com_agent<mesh_t, ExData>::process()
  {
    //std::cout << "process: jsensor : " << jsensor << std::endl;
    if (jsensor == nullptr) return false;
    //std::cout << "process: before cast" << std::endl;
    mesh_t& mesh = const_cast<mesh_t&>(jsensor->get_hmesh());
    //std::cout << "process: before sizes" << std::endl;
    E_Int npgs0 = mesh._ng.PGs.size();
    E_Int nphs0 = mesh._ng.PHs.size();
    //std::cout << "process: sizes:" << npgs0 << "/" << nphs0 << std::endl;
    //std::cout << "process: run adaptor" << std::endl;
    NUGA::adaptor<mesh_t, join_sensor_t>::run(mesh, *jsensor, false/*do_agglo*/); //fixme : agglo
    //std::cout << "after run" << std::endl;
    if (mesh._ng.PGs.size() != npgs0) return true;
    if (mesh._ng.PHs.size() != nphs0) return true;
    return false;
  }

}

#endif
