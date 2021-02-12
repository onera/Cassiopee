/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef _DELAUNAY_MESHER_MODE_H_
#define _DELAUNAY_MESHER_MODE_H_

#include "Nuga/include/defs.h"

namespace DELAUNAY
{
  struct MesherMode
  {
    enum eMeshMode {TRIANGULATION_MODE, REFINE_MODE};
    enum eInterpolType {LINEAR = 0, GEOMETRIC};

    MesherMode():mesh_mode(REFINE_MODE), hmin(NUGA::FLOAT_MAX), hmax(-NUGA::FLOAT_MAX),
                 do_not_shuffle(false), remove_holes(true), silent_errors(false), ignore_coincident_nodes(false), ignore_unforceable_edges(false),
                 metric_interpol_type(LINEAR), growth_ratio(1.2), nb_smooth_iter(0), symmetrize(false){}

    eMeshMode     mesh_mode;
    E_Float       hmin;
    E_Float       hmax;

    E_Bool        do_not_shuffle;
    E_Bool        remove_holes;
    E_Bool        silent_errors;
    E_Bool        ignore_coincident_nodes;
    E_Bool        ignore_unforceable_edges;
    
    eInterpolType metric_interpol_type;
    E_Float       growth_ratio;
    E_Int         nb_smooth_iter;
    E_Bool        symmetrize;
    

  };

  struct SurfaceMesherMode : public MesherMode
  {
    
    SurfaceMesherMode& operator=(const SurfaceMesherMode& rhs)
    {
      mesh_mode = rhs.mesh_mode;
      hmin = rhs.hmin;
      hmax = rhs.hmax;

      do_not_shuffle = rhs.do_not_shuffle;
      remove_holes = rhs.remove_holes;
      silent_errors = rhs.silent_errors;
      ignore_coincident_nodes = rhs.ignore_coincident_nodes;
      ignore_unforceable_edges = rhs.ignore_unforceable_edges;
      if (ignore_unforceable_edges)ignore_coincident_nodes=true;
    
      metric_interpol_type = rhs.metric_interpol_type;
      growth_ratio = rhs.growth_ratio;
      nb_smooth_iter = rhs.nb_smooth_iter;
      symmetrize = rhs.symmetrize;
    
      chordal_error = rhs.chordal_error;
      metric_mode = rhs.metric_mode;
      return *this;
    }
    
    enum GMmode
    {
      ISO_CST, ///< A constant size is specified to mesh the surface.
      ISO_RHO, ///< local minimum curvature radius is used to compute the metric.
      ANISO    ///< both local principal curvature radii are used to compute the metric.
    };

    SurfaceMesherMode():MesherMode(), chordal_error(0.01), metric_mode(ISO_RHO){}

    E_Float     chordal_error;
    GMmode      metric_mode;
  };
}

#endif

