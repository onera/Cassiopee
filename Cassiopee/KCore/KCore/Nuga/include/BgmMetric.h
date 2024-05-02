/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __GENERATOR_BGM_METRIC_H__
#define __GENERATOR_BGM_METRIC_H__

#include "Metric.h"
#include "MeshData.h"
#include "MeshTool.h"
#include "Nuga/include/KdTree.h"

namespace DELAUNAY
{

  template <typename T>
  class BgmMetric : public VarMetric<T>
  {
  public:
    typedef     std::vector<T>          field_type;
    typedef     VarMetric<T>            parent_type;    

  public:

    BgmMetric(K_FLD::FloatArray& pos, const field_type& metric, const DELAUNAY::MeshData& bgm, const Interpolator<T>& interpol)
      :parent_type (pos, metric, interpol), _bgm(bgm), _node_detector(pos), _meshtool(bgm, _node_detector){}

    ~BgmMetric(void){}

    E_Float length (size_type Ni, size_type Nj);

    inline void computeMetric(size_type N, size_type Ni, size_type Nj, E_Float r);

    inline void setMetric(E_Int N, const T& m);

  private:

    void __setMetrics(const std::vector<size_type>& Qs, const std::vector<T>& Ms);

    void __setMetric(size_type N){/*fixme*/_metric.push_back(_metric[0]);}

  private:
    MeshData              _bgm;
    K_SEARCH::KdTree      _node_detector;
    MeshTool              _meshtool;

  };

  template <typename T>
  E_Float
    BgmMetric<T>::length (size_type Ni, size_type Nj){

      int_vector_type        Qs;
      std::vector<T>         Ms;

      _meshtool.getIntersectionNodes(Ni, Nj, Qs, Ms);

      assert (Qs.size() > 1); // at least Ni and Nj;

      __setMetrics(Qs, Ms);

      E_Float d = 0;
      size_type sz = (size_type)Qs.size() - 1;
      for (size_type i = 0; i < sz; ++i)
        d += parent_type::length(Qs[i], Qs[i+1]);

      return d;
  }

  template <typename T>
  void
    BgmMetric<T>::__setMetrics(const std::vector<size_type>& Qs, const std::vector<T>& Ms)
  {
    std::vector<T> ms;
    ms.resize(_pos.cols(), _metric[0]);//fixme
    _metric.resize(_pos.cols(), _metric[0]);//fixme

    size_type sz = (size_type)Qs.size();
    for (size_type i = 0; i < sz; ++i)
      _metric[Qs[i]] = ms[i];//fixme
  }

}

#endif
