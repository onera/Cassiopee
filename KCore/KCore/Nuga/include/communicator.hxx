/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_COMMUNICATOR_HXX
#define NUGA_COMMUNICATOR_HXX

#include <vector>
#include "Nuga/include/join_sensor.hxx"
#include "Nuga/include/join_t.hxx"

namespace NUGA
{
 
  ///
  template <typename ExData>
  struct com_slots
  {
    ExData inslot;
    ExData outslot;
  };


  struct com_agent
  {
    E_Int id; // agent id

    com_agent(E_Int i) :id(i) {};

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

    const mesh_t& Smesh = Sag.jsensor->get_hmesh();

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
    if (jsensor == nullptr) return false;
    const mesh_t& mesh = jsensor->get_hmesh();
    bool has_packs{ false };

    for (auto& j : join->jid_to_joinlist)
    {
      E_Int jid = j.first;
      auto& ptlist = j.second;

      jid_to_boxes[jid].outslot.clear();

      for (size_t i = 0; i < ptlist.size(); ++i)
      {
        K_FLD::IntArray p;
        E_Int PGi = ptlist[i] - join->idx_start;
        mesh.extract_plan(PGi, true/*reverse*/, 0/*because previous sync*/, p);
        if (p.cols())
        {
          jid_to_boxes[jid].outslot[i] = p; //assume same sorting of joinlists
          has_packs = true;
        }
      }
    }
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
    if (jsensor == nullptr) return false;
    mesh_t& mesh = const_cast<mesh_t&>(jsensor->get_hmesh());

    E_Int npgs0 = mesh._ng.PGs.size();
    E_Int nphs0 = mesh._ng.PHs.size();

    NUGA::adaptor<mesh_t, join_sensor_t>::run(mesh, *jsensor, false/*do_agglo*/); //fixme : agglo

    if (mesh._ng.PGs.size() != npgs0) return true;
    if (mesh._ng.PHs.size() != nphs0) return true;
    return false;
  }

}

#endif
