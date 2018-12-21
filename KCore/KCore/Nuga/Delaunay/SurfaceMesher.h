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

#ifndef __SURFACE_MESHER_H__
#define __SURFACE_MESHER_H__

#include "MeshData.h"
#include "MesherMode.h"
#include "Mesher.h"
#include "GeomMetric.h"

namespace DELAUNAY
{
template <typename SurfaceType>
class SurfaceMesher
{
public:
  typedef K_CONT_DEF::size_type size_type;
public:
  SurfaceMesher(){}

  SurfaceMesher(const SurfaceMesherMode& mode):mode(mode)
  {
#ifdef DEBUG_MESHER
    dbg_flag = false;
#endif
  };
  ~SurfaceMesher(void){}
  E_Int run (SurfaceMeshData<SurfaceType>& data);
  
  SurfaceMesherMode  mode;

private:
  void __mapToSurface(const SurfaceType& surface, const K_FLD::FloatArray& pos2D, K_FLD::FloatArray& pos3D);
    
#ifdef DEBUG_MESHER
public:
    bool dbg_flag;
#endif
};

template <typename SurfaceType>
E_Int
SurfaceMesher<SurfaceType>::run(SurfaceMeshData<SurfaceType>& data)
{
  typedef GeomMetric<Aniso2D, SurfaceType>  MetricType;
  typedef Mesher<Aniso2D, MetricType>       MesherType;

  MetricType metric_aniso(*data.pos, data.surface, (typename MetricType::GMmode)mode.metric_mode,
                          mode.chordal_error, mode.hmin, mode.hmax, mode.growth_ratio);
  
  metric_aniso.set_pos2D(*data.pos);// hack to avoid to create an argument for init_metric (that is not required for other metric than GeomMetric

  //std::cout << "run 0" << std::endl;
  metric_aniso.init_metric(data.metrics, data.pos3D, *data.connectB, data.hardNodes);
  
  MesherType mesher(metric_aniso, mode);
  //std::cout << "run 1" << std::endl;
  
#ifdef DEBUG_MESHER
  mesher.dbg_flag = dbg_flag;
#endif

  E_Int err = mesher.run(data);

  //std::cout << "run 2" << std::endl;

  if (!err && (data.connectM.cols() != 0))//fixme : pos3D might have data upon entry so need to preserve them
  {
    //std::cout << "run 3" << std::endl;
    __mapToSurface(data.surface, *data.pos, data.pos3D);
    //std::cout << "run 4" << std::endl;

    // Replace the real contour coordinates.
    /*std::vector<E_Int> Bnodes;
    data.connectB.uniqueVals(Bnodes);
    E_Int Ni;
    for (E_Int i = 0; i < Bnodes.size(); ++i)
    {
      Ni = Bnodes[i];
      for (E_Int j = 0; j < 3; ++j)
        data.pos3D(j, Ni) = data.surface._pos(j,Ni);
    }*/
  }
  //std::cout << "run 5" << std::endl;

#ifdef WIN32
#ifdef E_DEBUG
  MIO::write("param.mesh", data.pos, data.connectM);
#endif
#endif

  return err;
}

template <typename SurfaceType>
void
SurfaceMesher<SurfaceType>::__mapToSurface
(const SurfaceType& surface, const K_FLD::FloatArray& pos2D, K_FLD::FloatArray& pos3D)
{
  E_Float          pt[3];
  size_type        COLS = pos2D.cols(), col0 = pos3D.cols();

  //pos3D.clear();

  pos3D.resize(3, COLS);

  for (size_type c = col0; c < COLS; ++c)
  {
    surface.point(pos2D(0,c), pos2D(1,c), pt);
    pos3D(0,c) = pt[0];
    pos3D(1,c) = pt[1];
    pos3D(2,c) = pt[2];
  }
}

}

#endif
