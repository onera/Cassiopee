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

#include "CADviaOCC.h"
// IGES/STEP
#include "IGESControl_Reader.hxx" 
#include "STEPControl_Reader.hxx"

//Data structure
#include "TColStd_HSequenceOfTransient.hxx"
#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GCPnts_AbscissaPoint.hxx" 
#include "GCPnts_UniformDeflection.hxx"
#include "GCPnts_UniformAbscissa.hxx"
#include "TopExp_Explorer.hxx"

#include "ShapeUpgrade_ShapeConvertToBezier.hxx" //to deal with surface of Revolution
#include "ShapeUpgrade_ShapeDivideAngle.hxx" //to deal with surface of Revolution
#include "Geom_SurfaceOfRevolution.hxx"
#include "ShapeCustom.hxx"
#include "Geom_TrimmedCurve.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"

#include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/merge.h"
#include "Nuga/include/SurfaceMesher.h"

#include "OCCSurface.h"
#include "Nuga/include/KdTree.h"

#include "Nuga/include/ContourSplitter.h"
#include "Nuga/include/BARSplitter.h"
#include "Nuga/include/BbTree.h"
#include "Nuga/include/FittingBox.h"
#include "Nuga/include/MeshUtils1D.h"
#include <Precision.hxx>

// algo=2
// Brancher __eval_nb_points2
//#define DEBUG_CAD_READER
// si debug, copier libconverter.a dans Dist
// decommenter dans setup.scons et setupScons

#ifdef DEBUG_CAD_READER
#include "Nuga/include/medit.hxx"
#include <sstream>
#endif

// Parametrise les edges et appelle le mailleur par face
E_Int K_OCC::CADviaOCC::mesh_faces2
(const K_FLD::FloatArray& coords, const std::vector<K_FLD::IntArray>& connectBs, std::vector<K_FLD::FloatArray>& crds, std::vector<K_FLD::IntArray>& connectMs, bool aniso, bool do_join)
{
  E_Int err(0), nb_faces(_surfs.Extent());
  
  if (!nb_faces) return 0;

  std::vector<K_FLD::FloatArray> crds1(nb_faces);
  std::vector<K_FLD::IntArray> connectMs1(nb_faces);
  
  std::vector<E_Int> nodes, nids;
  K_FLD::FloatArray UVcontour, pos3D;
  
  DELAUNAY::SurfaceMesher<OCCSurface> mesher;
  
#ifdef DEBUG_CAD_READER
  E_Int faulty_id = 16;
#endif
  //size_t t;

//#pragma omp parallel for private(nodes, nids, UVcontour, pos3D, err, mesher, t)
  for (E_Int i=1; i <= nb_faces; ++i)
  {
    std::cout << "Processing face: " << i << " / "<< nb_faces << std::endl;
    fflush(stdout);
    
    const OCCSurface& F = *_faces[i];
    
    K_FLD::IntArray connectB = connectBs[i];
    
    if (connectB.cols() == 0)
    {
#ifdef DEBUG_CAD_READER
      std::cout << "ERROR Face: " << i << ": empty discretized contour!" << std::endl;
#endif
      continue;
    }
    
#ifdef DEBUG_CAD_READER
    if (i == faulty_id)
      MIO::write("connectB.mesh", coords, connectB);
#endif
    
    connectB.uniqueVals(nodes);
    
    if (nodes.size() == 2) continue;
    
    if (nodes.size() == 3)
    {
      crds1[i-1].pushBack(coords.col(nodes[0]), coords.col(nodes[0])+3);
      crds1[i-1].pushBack(coords.col(nodes[1]), coords.col(nodes[1])+3);
      crds1[i-1].pushBack(coords.col(nodes[2]), coords.col(nodes[2])+3);
      E_Int T[] = {crds1[i-1].cols()-3, crds1[i-1].cols()-2, crds1[i-1].cols()-1};
      connectMs1[i-1].pushBack(T, T+3);
      continue;
    }
    
    // compact to mesh
    nids.clear();
    pos3D = coords;
    NUGA::MeshTool::compact_to_mesh(pos3D, connectB, nids);
    
    // added shrink pour forcer les pts a l'interieur de la surface
    //_faces[i]->shrink(pos3D, 0.9);
  
#ifdef DEBUG_CAD_READER
    if (i == faulty_id)
      MIO::write("connectBcompacted.mesh", pos3D, connectB, "BAR");
#endif
    
    /* Projection de pos3D sur la surface */
    //_faces[i]->project(pos3D);
    
#ifdef DEBUG_CAD_READER
    //if (i == faulty_id)
    //  MIO::write("connectBprojected.mesh", pos3D , connectB, "BAR");
#endif
    
#ifdef DEBUG_CAD_READER
    if (i == faulty_id)
    {
      K_FLD::FloatArray surfc;
      K_FLD::IntArray con;
      _faces[i]->discretize(surfc, con, 30, 30);
      MIO::write("surface.plt", surfc, con, "QUAD");
    }
#endif
      
    // surface of revolution => duplicate, reverse and separate seams
    //bool is_of_revolution = ((E_Int)nodes.size() != connectB.cols());
    bool is_of_revolution = (nodes.size() != (size_t)connectB.cols());
#ifdef DEBUG_CAD_READER
    if (i == faulty_id && is_of_revolution) printf("Edge is not loop - Is of revolution=true\n");
#endif
    //is_of_revolution = false; // DBG
    
    std::map<E_Int, std::pair<E_Int, E_Int> > seam_nodes;
    
    if (is_of_revolution)
    {  
      _faces[i]->_isRevol = true; // trigger jump test and branch correction
      //__split_surface_of_revolution2(_faces[i], connectB, pos3D, seam_nodes);
      
      // Parcours la BAR pour n'avoir qu'une seule loop
      //_faces[i]->parcoursBAR(pos3D, connectB);
      
      K_FLD::FloatArray pos3Dorig = pos3D;
      K_FLD::IntArray connectBorig = connectB;
      K_FLD::IntArray switcha(E_Int(1),pos3D.cols(),E_Int(0));
      std::map<E_Int, E_Int> mirror;
      
      _faces[i]->dupBAR(pos3D, connectB, switcha, mirror);
      err = _faces[i]->parameters2(pos3D, connectB, UVcontour);
      
      if (err > 0)
      { 
        err = err-1;
#ifdef DEBUG_CAD_READER
        if (i == faulty_id) printf("Warning: switcha in=%d\n", err);
#endif
        if (err >= switcha.size()) err = mirror[err];
        switcha[err] = 1;
        pos3D = pos3Dorig; connectB = connectBorig; mirror.clear();
        _faces[i]->dupBAR(pos3D, connectB, switcha, mirror);
        err = _faces[i]->parameters2(pos3D, connectB, UVcontour);
      }
      
      if (err > 0)
      { 
        err = err-1;
#ifdef DEBUG_CAD_READER
        if (i == faulty_id) printf("Warning: switcha in=%d\n", err);
#endif
        if (err >= switcha.size()) err = mirror[err];
        switcha[err] = 1;
        pos3D = pos3Dorig; connectB = connectBorig; mirror.clear();
        _faces[i]->dupBAR(pos3D, connectB, switcha, mirror);
        err = _faces[i]->parameters2(pos3D, connectB, UVcontour);
      }
#ifdef DEBUG_CAD_READER
      if (err > 0) printf("still unresolved loop.\n");
#endif
      _faces[i]->_isRevol = false; // avoid jump check
      
#ifdef DEBUG_CAD_READER      
      if (i == faulty_id)
        MIO::write("connectBrevol.mesh", pos3D, connectB, "BAR");
#endif
    }
    //printf("before param\n"); fflush(stdout);
    
    // Up to 2 tries : first by asking OCC for params, Second by "hand" (sampling)
    err = 0;
    for (size_t t=0; t<1; ++t) // supp. la parametrisation discrete
    {
      if (t==0)
        err = _faces[i]->parameters2(pos3D, connectB, UVcontour);
      else
        err = _faces[i]->parametersSample(pos3D, UVcontour);
      
      
      if (!err) // Need to reorient holed surface.
        err = __reorient_holed_surface(connectB, UVcontour);
      
      
      if (err)
      {
        if (t==1)
          std::cout << "ERROR Face: " << i << ": cannot retrieve parametrization !" << std::endl;
        continue;
      }
      
#ifdef DEBUG_CAD_READER
      if (i==faulty_id && t==0)
      MIO::write("connectBUV1.mesh", UVcontour, connectB, "BAR");
      if (i==faulty_id && t==1)
      MIO::write("connectBUV2.mesh", UVcontour, connectB, "BAR");
#endif
    
      OCCSurface occ_surf(F);
      DELAUNAY::SurfaceMeshData<OCCSurface> data(UVcontour, pos3D, connectB, occ_surf);
    

#ifdef DEBUG_MESHER
      if (i == faulty_id)
        mesher.dbg_flag=true;
      else
        mesher.dbg_flag=false;
#endif
      DELAUNAY::SurfaceMesherMode mode;
      
      mode.chordal_error = _chordal_err;
      
      if (aniso) mode.metric_mode = mode.ANISO;

#ifndef DEBUG_CAD_READER
      //mode.silent_errors = true;
#endif
      if (_gr <= 0.) // unspecified == OLD MODE
      {
        mode.symmetrize = false;
        mode.hmax = _h;
        mode.growth_ratio = 0.;
      }
      else
      {
        mode.growth_ratio = std::max(_gr, 1.); //fixme : values in [1.; 1.25] might cause errors. Anyway do not allow bellow 1. since not useful
        E_Int MIN_NB = 20;
        if (_gr > 1. && (connectB.cols() > (MIN_NB * 4))) // check if growth ratio is applicable to this patch, i.e. it is symmetrizable
        {
          // Count the minimum number of edges on a boundary of the param space per direction (U,V)
          E_Float nu(0), nv(0);
          for (E_Int i=0; i < connectB.cols(); ++i)
          {
            E_Float ds[2];
            const E_Int& Ni = connectB(0,i);
            const E_Int& Nj = connectB(1,i);
            K_FUNC::diff<2>(UVcontour.col(Ni), UVcontour.col(Nj), ds);
            K_FUNC::normalize<2>(ds);
            
            if (::fabs(ds[0]) < ::fabs(ds[1])) ++nv;
            else ++nu;
          }
          
          E_Int n = std::min(nu, nv) / 2; // nu and nv are accumlation for "2" sides
          
          if (n > MIN_NB)
          {
            mode.symmetrize = true;
            mode.nb_smooth_iter = 2;
          }
        }
      }
      mesher.mode = mode;
      
      //printf("before mesher\n"); fflush(stdout);
      mesher.seed_random(1);
      err = mesher.run(data);
      
      if (err || (data.connectM.cols() == 0))
      {
        //if (t==0) continue; // pour lever l'erreur sur la param OCC
        if (err)
          std::cout << "ERROR Face : " << i << " : Geometric Mesher failed." << std::endl;
        else
          std::cout << "ERROR Face : " << i << " : Cannot retrieve parametrization (OCC Limitation) !" << std::endl;    
        continue;
      }
      
      // Join the seam for surface of revolution.
      if (!seam_nodes.empty())
      {
        std::vector<E_Int> nids;
        K_CONNECT::IdTool::init_inc(nids, data.pos3D.cols());
        
        for (auto &s : seam_nodes){
          assert (s.first != E_IDX_NONE);
          nids[s.second.first] = s.first;
          nids[s.second.second] = s.first;
        }
        
        K_FLD::IntArray::changeIndices(data.connectM, nids);
      }
    
      crds1[i-1] = data.pos3D;
      connectMs1[i-1] = data.connectM;
      
      if (!err) break; // done
      }

    } // End face loop

    if (do_join)
    {
      E_Int max_solid_id=0;
      for (E_Int i=1; i <= nb_faces; ++i)
        max_solid_id = std::max(_faces[i]->_parent, max_solid_id);
      
      crds.resize(max_solid_id+1);
      connectMs.resize(max_solid_id+1);
      
      for (E_Int i=1; i <= nb_faces; ++i)
      {
        const OCCSurface& F = *_faces[i];
        E_Int pid = F._parent;
        
        E_Int shft = crds[pid].cols();
        connectMs1[i-1].shift(shft);
        crds[pid].pushBack(crds1[i-1]);
        connectMs[pid].pushBack(connectMs1[i-1]);
      }
    
      std::vector<E_Int> nids;
      for (E_Int i=0; i <= max_solid_id; ++i)
      {
        if (crds[i].cols()==0) continue;

        K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(crds[i]);
        ::merge(crdA, _merge_tol, nids);
        K_FLD::IntArray::changeIndices(connectMs[i], nids);
        nids.clear();
        NUGA::MeshTool::compact_to_mesh(crds[i], connectMs[i], nids);
      }
    }
    else
    {
      connectMs = connectMs1;
      crds = crds1;
    }
  
  return 0;
}


// ma version
// Essai de decouper la surface periodique
void K_OCC::CADviaOCC::__split_surface_of_revolution2(const OCCSurface* face, K_FLD::IntArray& connectB, K_FLD::FloatArray& pos3D, std::map<E_Int, std::pair<E_Int, E_Int> >& seam_nodes)
{
  
#ifdef DEBUG_CAD_READER
  std::cout << "splitting surface ..." << std::endl;
#endif
  
  // 1. removing any duplicate
  std::set<K_MESH::NO_Edge> uedges;
  seam_nodes.clear();
  
  //assert ((face->_isUClosed && !face->_isVClosed) || (!face->_isUClosed && face->_isVClosed));
    
  for (E_Int i=0; i < connectB.cols(); ++i)
  {
    E_Int& N0 = connectB(0,i);
    E_Int& N1 = connectB(1,i);
    K_MESH::NO_Edge e(N0, N1);
    
    if (!uedges.insert(e).second) // already in => seam edge
    {
      E_Int N0 = e.node(0);
      if (seam_nodes.find(N0) == seam_nodes.end())
        __add_seam_node2(face, pos3D, N0, seam_nodes);
      E_Int N1 = e.node(1);
      if (seam_nodes.find(N1) == seam_nodes.end())
        __add_seam_node2(face, pos3D, N1, seam_nodes);
    }  
  }
  
  K_FLD::IntArray new_connect;
  
  for (auto& e : uedges)
  {
    E_Int N0 = e.node(0);
    E_Int N1 = e.node(1);
    
    auto it0 = seam_nodes.find(N0);
    auto it1 = seam_nodes.find(N1);
    
    if (it0 == seam_nodes.end() && it1 == seam_nodes.end()) //regular non-seam edge
      new_connect.pushBack(e.begin(), e.end());
    else if (it0 != seam_nodes.end() && it1 == seam_nodes.end()) // seam-connected edge
    {
      E_Float u,v,up=-1,vp,upp,vpp;
      face->parameters2(pos3D.col(N1), u, v, N1, up, vp, upp, vpp);
      
      E_Float* p = &u;
      if (face->_isVClosed && !face->_isUClosed) p = &v;
      
      E_Int e[] = {N1, E_IDX_NONE};
      if (*p < K_CONST::E_PI)
        e[1] = it0->second.second;
      else
        e[1] = it0->second.first;
      new_connect.pushBack(e, e+2);
    }
    else if (it0 == seam_nodes.end() && it1 != seam_nodes.end()) // seam-connected edge
    {
      E_Float u,v,up=-1,vp,upp,vpp;
      //E_Int err = 
      face->parameters2(pos3D.col(N0), u, v, N0, up, vp, upp, vpp);
      
      E_Float *p = &u;
      if (face->_isVClosed && !face->_isUClosed) p = &v;
      
      E_Int e[] = {N0, E_IDX_NONE};
      if (*p < K_CONST::E_PI)
        e[1] = it1->second.second;
      else
        e[1] = it1->second.first;
      new_connect.pushBack(e, e+2);
    }
    else //seam-edge
    {
      E_Int e[] = {it0->second.first, it1->second.first};
      new_connect.pushBack(e, e+2);
      E_Int e1[] = {it0->second.second, it1->second.second};
      new_connect.pushBack(e1, e1+2); 
    }
  }
  
  connectB = new_connect;

  // Ajout compactage
  /*
  std::vector<E_Int> nids;
  K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(pos3D);
  E_Int nb_merge = ::merge(crdA, _merge_tol, nids);
  K_FLD::IntArray::changeIndices(connectB, nids);
  printf("nbmerge=%d\n", nb_merge);
  nids.clear();
  NUGA::MeshTool::compact_to_mesh(pos3D, connectB, nids);
  */
  
  /*
  for (E_Int i = 0; i < npts; i++)
  {
    ind1 = newId[i];
    if (ind1 == i) // not merged
    {
      indirp[ind1] = np; indir2p[np] = ind1; np++;
    }
  }
  */
  
#ifdef DEBUG_CAD_READER
  MIO::write("revol.mesh", pos3D, connectB);
#endif  
}

void K_OCC::CADviaOCC::__add_seam_node2
(OCCSurface const *face, K_FLD::FloatArray& pos3D, E_Int N0,
 std::map<E_Int, std::pair<E_Int, E_Int> >& seam_nodes)
{ 
  //create the points
  E_Float u,v,up=-1,vp,upp,vpp;
  face->parameters2(pos3D.col(N0), u, v, N0, up, vp, upp, vpp);

  E_Float Pt[3];

  if (face->_isUClosed && !face->_isVClosed)
  {
    const E_Float& u0 = face->_U0+0.01;
    const E_Float& u1 = face->_U1-0.01;

    face->point(u0, v, Pt);
    pos3D.pushBack(Pt, Pt+3);
    E_Int Nright=pos3D.cols()-1;
    seam_nodes[N0].second = Nright;

    face->point(u1, v, Pt);
    pos3D.pushBack(Pt, Pt+3);
    E_Int Nleft=pos3D.cols()-1;
    seam_nodes[N0].first = Nleft;
  }
  else if (!face->_isUClosed && face->_isVClosed)
  {
    const E_Float& v0 = face->_V0+0.01;
    const E_Float& v1 = face->_V1-0.01;
    
    face->point(u, v0, Pt);
    pos3D.pushBack(Pt, Pt+3);
    seam_nodes[N0].second = pos3D.cols()-1;

    face->point(u, v1, Pt);
    pos3D.pushBack(Pt, Pt+3);
    seam_nodes[N0].first = pos3D.cols()-1;
  }
  else
  {
#ifdef DEBUG_CAD_READER
    //todo
    std::cout << "to implement !" << std::endl;
    //assert (false);
#endif
  }
}

