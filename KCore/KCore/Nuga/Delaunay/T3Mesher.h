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

#ifndef __T3_MESHER_H__
#define __T3_MESHER_H__

#include "MeshData.h"
#include "MesherMode.h"
#include "Mesher.h"

namespace DELAUNAY
{
template <typename T>
class T3Mesher
{
public:
  T3Mesher()
#ifdef DEBUG_MESHER
   :dbg_flag(false)
#endif
  {}
  T3Mesher(MesherMode& mode):_mode(mode)
#ifdef DEBUG_MESHER
   ,dbg_flag(false)
#endif
  {}
  ~T3Mesher(void){}
  E_Int run (MeshData& data);
  
  MesherMode& _mode;
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
    
  typename MetricType::eInterpType interpol = (_mode.metric_interpol_type == _mode.LINEAR) ?
                                              MetricType::eInterpType::LINEAR : MetricType::eInterpType::GEOMETRIC;

  MetricType metric(*data.pos, _mode.hmin, _mode.hmax, interpol);

  //std::cout << data.metrics << std::endl;
  metric.init_metric(data.metrics, *data.pos, *data.connectB, data.hardNodes);

  MesherType mesher(metric, _mode);
  
  
#if defined(DEBUG_MESHER)
  mesher.dbg_flag=dbg_flag;
#endif

  E_Int err = mesher.run (data);
  
  _edge_errors = mesher._edge_errors;

  return err;
}

}

#endif
