/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __T3_MESHER_H__
#define __T3_MESHER_H__

#include "MeshData.h"
#include "MesherMode.h"
#include "Mesher.h"

namespace DELAUNAY
{
template <typename T>
class T3Mesher : public Mesher<T, VarMetric<T>>
{
public:
  using parent_t = Mesher<T, VarMetric<T>>;
public:
  T3Mesher()
#ifdef DEBUG_MESHER
   :dbg_flag(false)
#endif
  {}
  T3Mesher(MesherMode& mode):parent_t()
#ifdef DEBUG_MESHER
   ,dbg_flag(false)
#endif
  {
    parent_t::mode = mode;
  }
  ~T3Mesher(void){}
  E_Int run (MeshData& data);
  
#if defined(DEBUG_MESHER)
    public:
      bool dbg_flag;
#endif
  public:
    std::vector<edge_error_t>  _edge_errors;
};

template <typename T>
E_Int
T3Mesher<T>::run (MeshData& data)
{
  typedef VarMetric<T>          MetricType;
  typedef Mesher<T, MetricType> MesherType;
    
  typename MetricType::eInterpType interpol = (parent_t::mode.metric_interpol_type == parent_t::mode.LINEAR) ?
                                              MetricType::eInterpType::LINEAR : MetricType::eInterpType::GEOMETRIC;

  MetricType metric(*data.pos, parent_t::mode.hmin, parent_t::mode.hmax, interpol);

  //std::cout << data.metrics << std::endl;
  metric.init_metric(data.metrics, *data.pos, *data.connectB, data.hardNodes);

  parent_t::clear(); // clear container attributes
  parent_t::set(metric);

#if defined(DEBUG_MESHER)
  mesher.dbg_flag=dbg_flag;
#endif

  E_Int err = parent_t::run (data);

  return err;
}

}

#endif
