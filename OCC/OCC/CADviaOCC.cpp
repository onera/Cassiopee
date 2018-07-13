/*    
    Copyright 2013-2018 Onera.

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
//Author : Sâm Landier (sam.landier@onera.fr)


#include "CADviaOCC.h"

// IGES/STEP
#include "IGESControl_Reader.hxx" //TKIGES
//#include "STEPControl_Reader.hxx" // TKSTEP
//Data structure
#include "TColStd_HSequenceOfTransient.hxx"
#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GCPnts_AbscissaPoint.hxx" // TKGeomBase
#include "GCPnts_UniformDeflection.hxx"
#include "GCPnts_UniformAbscissa.hxx"
#include "TopExp_Explorer.hxx"

#include "ShapeUpgrade_ShapeConvertToBezier.hxx" //to deal with surface of Revolution
#include "ShapeUpgrade_ShapeDivideAngle.hxx" //to deal with surface of Revolution
#include "Geom_SurfaceOfRevolution.hxx"
#include "ShapeCustom.hxx"
#include "Geom_TrimmedCurve.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"

#include "Fld/ArrayAccessor.h"
#include "Connect/merge.h"
#include "Nuga/Delaunay/SurfaceMesher.h"

#include "OCCSurface.h"
#include "Search/KdTree.h"

#include "Connect/ContourSplitter.h"
#include "Connect/BARSplitter.h"
#include "Search/BbTree.h"
#include "Nuga/GapFixer/FittingBox.h"

/*
#include "TopoDS_Edge.hxx"
// Curve and Mesh
#include "GeomAdaptor_Curve.hxx"
#include "Geom_Line.hxx"
#include "gp_Lin.hxx"
*/

//
#ifdef DEBUG_CAD_READER
#include "IO/io.h"
#include <sstream>
#endif


K_OCC::CADviaOCC::CADviaOCC()
{
}


K_OCC::CADviaOCC::~CADviaOCC()
{
  for (size_t i=1; i < _faces.size(); ++i)
    delete _faces[i];
}

E_Int import_iges(const char* fname, TopoDS_Shape& sh)
{    
  // Read the file
  IGESControl_Reader reader;
  reader.ReadFile(fname);
    	
  // Transfer CAD faces (only) into a OCC list
  Handle(TColStd_HSequenceOfTransient) occ_list = reader.GiveList("iges-faces");
    
  if (occ_list.IsNull()) return 1;
    
  Standard_Integer nb_cad_faces = occ_list->Length();
  Standard_Integer nb_transfered_faces = reader.TransferList(occ_list);

#ifdef DEBUG_CAD_READER
  std::cout << "IGES Faces: " << nb_cad_faces << "   Transferred:" << nb_transfered_faces << endl;
#endif

  sh = reader.OneShape();
	
  return sh.IsNull();
}


int import_step(const char* fname, TopoDS_Shape& sh)
{
  // Read the file
  //STEPControl_Reader reader;
  //reader.ReadFile(fname);
	
//   // Transfer CAD faces (only) into a OCC list
//   Handle(TColStd_HSequenceOfTransient) occ_list = reader.GiveList("step-faces");
		
//   Standard_Integer nb_cad_faces = occ_list->Length();
//   Standard_Integer nb_transfered_faces = reader.TransferList(occ_list);

// #ifdef DEBUG_CAD_READER
//   cout << "STEP Faces: " << nb_cad_faces << "   Transferred:" << nb_transfered_faces << endl;
// #endif

//   sh = reader.OneShape();
	
//   return sh.IsNull();

  return 0;
}

//
E_Int K_OCC::CADviaOCC::import_cad(const char* fname, const char* format, E_Float h, E_Float chordal_err)
{
  _chordal_err = chordal_err;
  _h = h;
    
  E_Int err(1);
  if (::strcmp(format, "iges")==0)
    err = import_iges(fname, _occ_shape);
  /*else if (::strcmp(format, "step")==0)
    err = import_step(fname, _occ_shape);*/
  
  if (err) return err;
    
  return __build_graph(_occ_shape, _faces);
}

E_Int K_OCC::CADviaOCC::compute_h_sizing(K_FLD::FloatArray& coords, std::vector<E_Int>& Ns)
{
  E_Int err(0), nb_edges(_edges.Extent());
  
  if (!nb_edges) return 0;
  
  Ns.resize(nb_edges+1, 0);
  
  std::vector<E_Float> Ls;
  Ls.resize(nb_edges+1, 0.);
  
  _Lmin = K_CONST::E_MAX_FLOAT;
  _Lmean = _Lmax = 0.;
  
  for (E_Int i=1; i <= nb_edges; ++i)
  {
    __h_sizing(TopoDS::Edge(_edges(i)), Ls[i]);
    _Lmin = std::min(Ls[i], _Lmin);
    _Lmax = std::max(Ls[i], _Lmax);
    _Lmean += Ls[i];
  }
  
  _Lmean /= nb_edges;
  
#ifdef DEBUG_CAD_READER
  std::cout << "Lmin/Lmean/Lmax : " << _Lmin << "/" << _Lmean << "/" << _Lmax << std::endl;
#endif
  
  if (_h <= 0.) // undefined or badly defined
    _h = _Lmean/10.;
  
  for (E_Int i=1; i <= nb_edges; ++i)
  {
    Ns[i] = (E_Int)(Ls[i]/_h) + 1;
  }
  
  return err;
}

E_Int K_OCC::CADviaOCC::update_with_chordal_sizing(std::vector<E_Int>& Ns)
{
  E_Int err(0), nb_edges(_edges.Extent());
  
  if (!nb_edges) return 0;
  
  Ns.resize(nb_edges+1, 0);
  
  if (_chordal_err <= 0.) //undefined or badly defined
    _chordal_err = _h/50.;

  for (E_Int i=1; i <= nb_edges; ++i)
    __chord_sizing(TopoDS::Edge(_edges(i)), _chordal_err, Ns[i]);

  return err;
}

///
E_Int K_OCC::CADviaOCC::__h_sizing(const TopoDS_Edge& E, E_Float& L)
{  
  if (BRep_Tool::Degenerated (E))      // Exit if the edge is degenerated.
    return 1;
  
  BRepAdaptor_Curve C0(E);
  /*gp_Lin line = C0.Line();
  std::cout << "dir : " << line.Direction().X() << "/" << line.Direction().Y() << "/" << line.Direction().Z() << "/" << std::endl;*/
    
  GeomAdaptor_Curve geom_adap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  L = (E_Float) GCPnts_AbscissaPoint::Length(geom_adap, geom_adap.FirstParameter(), geom_adap.LastParameter());
  
  return 0;
}

///
E_Int K_OCC::CADviaOCC::__chord_sizing(const TopoDS_Edge& E, E_Float chordal_err, E_Int& nb_points)
{  
  if (BRep_Tool::Degenerated (E))      // Exit if the edge is degenerated.
    return 1;
  
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geom_adap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  
  GCPnts_UniformDeflection unif_defl(geom_adap, chordal_err);
  if (!unif_defl.IsDone())
    return 1;
  
  nb_points = std::max(nb_points, unif_defl.NbPoints()); // choose finer
   
#ifdef DEBUG_CAD_READER
  //std::cout << "nb poins : " << nb_points << std::endl; 
  //assert(nb_points > 1);
#endif
  
  return 0;
}

///
E_Int K_OCC::CADviaOCC::mesh_edges(K_FLD::FloatArray& coords, std::vector<K_FLD::IntArray>& connectEs)
{
  E_Int err(0), nb_edges(_edges.Extent());
  
  if (!nb_edges) return 0;
  
  std::vector<E_Int> Ns;
  compute_h_sizing(coords, Ns);   // compute the number of points per edge regarding h
  update_with_chordal_sizing(Ns); // take chordal error into account if it gives a finer discretization 
  
  connectEs.resize(nb_edges+1);
  
  for (E_Int i=1; i <= _surfs.Extent(); ++i)
  {
    //std::cout << " cad face nb : " << i << std::endl;
    
    if (!_faces[i]) // A priori Surface of revolution
      continue;
    
    const OCCSurface& F = *_faces[i];
    
    E_Int nb_edges = F._edges.size();
    
    for (size_t j=0; j < nb_edges; ++j)
    {
      E_Int id = ::abs(F._edges[j]);
      
      if (connectEs[id].cols())
        continue;//already meshed.
      
#ifdef DEBUG_CAD_READER
    std::cout << "Edge :  " << id <<  ". Nb points : " << Ns[id] << std::endl;
#endif
      __mesh_edge(TopoDS::Edge(_edges(id)), Ns[id], coords, connectEs[id]);
      
#ifdef DEBUG_CAD_READER
    if (id==41 || id==101)
      MIO::write("Ei.mesh", coords , connectEs[id], "BAR");
#endif
      
    }
  }
    
#ifdef DEBUG_CAD_READER
  K_FLD::IntArray tmp;
  for (size_t i=0; i <connectEs.size(); ++i)
    tmp.pushBack(connectEs[i]);
  MIO::write("wire.mesh", coords, tmp, "BAR");
#endif
  
  return err;
}

E_Int K_OCC::CADviaOCC::__remove_degenerated(K_FLD::IntArray& connectE)
{
  E_Int                       Si, COLS(connectE.cols()), ROWS(2);
  K_FLD::IntArray::iterator   pS;
  K_FLD::IntArray             connectOut;
  
  connectOut.reserve(ROWS, COLS);

  //
  for (Si = 0; Si < COLS; ++Si)
  {
    pS = connectE.col(Si);
    
    if (*pS != *(pS+1))
      connectOut.pushBack(pS, pS+ ROWS);
  }

  connectE = connectOut;
 
  return (COLS - connectE.cols());
}

//
E_Int K_OCC::CADviaOCC::build_loops
(K_FLD::FloatArray& coords, const std::vector<K_FLD::IntArray>& connectEs, std::vector<K_FLD::IntArray>& connectBs)
{
  E_Int nb_faces(_surfs.Extent());
  
  if (!nb_faces) return 0;
  
  connectBs.resize(nb_faces+1);
  
  K_FLD::IntArray tmp;
  K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(coords);
  std::vector<E_Int> end_nodes;
  E_Float tol2=-1.;
  
#ifdef DEBUG_CAD_READER
  E_Int faulty_id=69;
#endif
  
  for (E_Int i=1; i <= nb_faces; ++i)
  {
    //std::cout << " cad face nb : " << i << std::endl;
    
    if (!_faces[i]) // A priori Surface of revolution
      continue;
    
    const OCCSurface& F = *_faces[i];
    
    std::vector<E_Int> unods;
        
    // prepare contour  
    E_Int nb_edges = F._edges.size();
    //
    for (size_t j=0; j < nb_edges; ++j)
    {
      E_Int id = F._edges[j];
      
      unods.clear();
      
      if (id > 0)
      {
        if (connectEs[id].cols()==0)
          continue;
        connectBs[i].pushBack(connectEs[id]);
        connectEs[id].uniqueVals(unods);
        
      }
      else
      {
        tmp=connectEs[-id];
        if (tmp.cols()==0)
          continue;
        tmp.uniqueVals(unods);
        for (size_t k=0; k < tmp.cols(); ++k)std::swap(tmp(0,k), tmp(1,k));
        connectBs[i].pushBack(tmp);
      }
    }

#ifdef DEBUG_CAD_READER
    if (faulty_id==i)
      MIO::write("connectBi.mesh", coords , connectBs[i], "BAR");
#endif
    
    // clean it
    while (__clean(crdA, connectBs[i], tol2))
    ;
    
#ifdef DEBUG_CAD_READER
    if (faulty_id==i)
      MIO::write("connectBc.mesh", coords , connectBs[i], "BAR");
    //std::ostringstream o;
    //o << "clean_loop_" << i << ".mesh";
    //meshIO::write(o.str().c_str(), coords, connectBs[i]);
#endif
  }
    
  // Global pass to join the loops.  
  std::vector<E_Int> nids;
  _merge_tol = E_EPSILON;//::sqrt(tol2) + E_EPSILON;
  E_Int nb_merges = ::merge(crdA, _merge_tol, nids);
  
  if (!nb_merges)
    return 0;

  for (E_Int i=1; i <= nb_faces; ++i)
  {
    K_FLD::IntArray::changeIndices(connectBs[i], nids);
    __remove_degenerated(connectBs[i]);
    
#ifdef DEBUG_CAD_READER
    //meshIO::write("connectBf.mesh",coords , connectBs[i]);
#endif
  }
  
#ifdef DEBUG_CAD_READER
  /*{
  K_FLD::IntArray tmp;
  for (E_Int i=1; i <= nb_faces; ++i)
    tmp.pushBack(connectBs[i]);
  meshIO::write("cleaned_loops.mesh", coords, tmp);
  }*/
#endif
  
  return 0;
}

//
E_Int K_OCC::CADviaOCC::__clean(const K_FLD::ArrayAccessor<K_FLD::FloatArray>& crdA, K_FLD::IntArray& connectB, E_Float& tol2)
{ 
  //
  _end_nodes.clear();
  _enodes.clear();
  
  for (size_t i=0; i < connectB.cols(); ++i)
  {
    if (_enodes.find(connectB(0,i)) == _enodes.end())
      _enodes[connectB(0,i)]=1;
    else
      _enodes[connectB(0,i)]+=1;
    if (_enodes.find(connectB(1,i)) == _enodes.end())
      _enodes[connectB(1,i)]=1;
    else
      _enodes[connectB(1,i)]+=1;
  }
  
  for (std::map<E_Int, E_Int>::const_iterator it = _enodes.begin(); it != _enodes.end(); ++it)
    if (it->second == 1)
      _end_nodes.push_back(it->first);
  
  if (_end_nodes.empty()) // the contour is clean.
    return 0;
    
#ifdef DEBUG_CAD_READER
  /*K_FLD::FloatArray coords;
  K_FLD::IntArray dummy;//(2,end_nodes.size()/2);//an dummy edge to show the point cloud in medit
  E_Int sz = end_nodes.size()/2, count(0);
  for (size_t k=0; k < sz; k++)
  {
    E_Int e[]={count++,count++};
    dummy.pushBack(e,e+2);
  }
  coords.append_selection(crdA.array(), end_nodes);
  meshIO::write("merge_cloud.mesh", coords, dummy);*/
#endif
  
  // build the KdTree.
  K_SEARCH::KdTree<K_FLD::FloatArray> tree(crdA, _end_nodes);
  
  typedef std::vector<std::pair <E_Float, std::pair <E_Int, E_Int> > > palmares_t;
  palmares_t palma;
  E_Float d2;
  E_Int N;
  
  for (size_t i = 0; i < _end_nodes.size(); ++i)
  {
    const E_Int& N1 = _end_nodes[i];
    N = tree.getClosest(N1, d2);
    palma.push_back(std::make_pair(d2, std::make_pair(N1, N)));
  }
  
  std::sort(palma.begin(), palma.end());
  
  // move the node
  E_Int Fi, Mi, nb_merges(0);
  size_t psz = palma.size();
  std::vector<bool> mergeable(crdA.size(), true);
  
  std::vector<E_Int> nids;
  nids.resize(crdA.size());
  for (size_t i = 0; i < nids.size(); ++i) nids[i]=i;
  //
  for (size_t i = 0; i < psz; ++i)
  {
    Fi=palma[i].second.first;
    Mi=palma[i].second.second;

    // rule of merging : if target and moving has not been moved already
    if (mergeable[Fi] && mergeable[Mi])
    {
      nids[Mi]=Fi;
      ++nb_merges;
      mergeable[Fi]=mergeable[Mi]=false;
      tol2 = std::max(tol2, palma[i].first); // set the global tolerance
    }
  }
  
  //
  K_FLD::IntArray::changeIndices(connectB, nids);
  
  return 1; 
}

//
E_Int K_OCC::CADviaOCC::mesh_faces
(K_FLD::FloatArray& coords, const std::vector<K_FLD::IntArray>& connectBs, std::vector<K_FLD::FloatArray>& crds, std::vector<K_FLD::IntArray>& connectMs)
{
  E_Int err(0), nb_faces(_surfs.Extent());
  
  if (!nb_faces) return 0;
  
  std::vector<E_Int> nodes, nids;
  K_FLD::FloatArray UVcontour, pos3D;
  
  DELAUNAY::SurfaceMesherMode mode;
  mode.chordal_error=_chordal_err;
  //mode.hmin=_h;
  mode.hmax=_h;
  DELAUNAY::SurfaceMesher<OCCSurface> mesher(mode);
  
  E_Int max_solid_id=0;
  for (E_Int i=1; i <= nb_faces; ++i)
  {
    if (!_faces[i]) // A priori Surface of revolution
      continue;
    max_solid_id = std::max(_faces[i]->_parent, max_solid_id);
  }
  
  connectMs.resize(max_solid_id+1);
  crds.resize(max_solid_id+1);
  
#ifdef DEBUG_CAD_READER
  E_Int faulty_id = 5;
#endif
  
  for (E_Int i=1; i <= nb_faces; ++i)
  {
#ifdef DEBUG_CAD_READER
    std::cout << " face nb : " << i << std::endl;
#endif
    
    if (!_faces[i]) // A priori Surface of revolution
      continue;
    
    const OCCSurface& F = *_faces[i];
    
    K_FLD::IntArray connectB = connectBs[i];
    
    if (connectB.cols() == 0)
    {
#ifdef DEBUG_CAD_READER
      std::cout << "EEROR Face : " << i << " : empty discretized contour!" << std::endl;
#endif
      continue;
    }
    
#ifdef DEBUG_CAD_READER
    //if (i==faulty_id)
      //meshIO::write("connectB.mesh",coords , connectB);
#endif
    
    connectB.uniqueVals(nodes);
    
#ifdef DEBUG_CAD_READER
    //assert (nodes.size() == connectB.cols()); // cleaning has been done to close the loops
    if (nodes.size() != connectB.cols())
    {
      std::cout << "EEROR Face : " << i << " : unclosed contour!" << std::endl;
      //E_Float gap = ::sqrt(K_FUNC::sqrDistance(coords.col(1459), coords.col(1464), 3));
      //std:: cout << "gap : " << gap << std::endl;
      continue;
    }
#endif
    
    if (nodes.size() == 2)
      continue;
    
    if (nodes.size() == 3)
    {
      crds[F._parent].pushBack(coords.col(nodes[0]), coords.col(nodes[0])+3);
      crds[F._parent].pushBack(coords.col(nodes[1]), coords.col(nodes[1])+3);
      crds[F._parent].pushBack(coords.col(nodes[2]), coords.col(nodes[2])+3);
      E_Int T[] = {crds[F._parent].cols()-3, crds[F._parent].cols()-2, crds[F._parent].cols()-1};
      connectMs[F._parent].pushBack(T, T+3);
      continue;
    }
    
    // compact to mesh
    nids.clear();
    pos3D=coords;
    K_CONNECT::MeshTool::compact_to_mesh(pos3D, connectB, nids);
    
#ifdef DEBUG_CAD_READER
    //if (i==faulty_id)
      MIO::write("connectBcompacted.mesh",pos3D , connectB, "BAR");
#endif
      
    // Up to 2 tries : first by asking OCC for params, Second by "hand" (sampling)
    err = 0;
    for (size_t t=0; t<2; ++t)
    {
      if (t==0)
        err=_faces[i]->parameters(pos3D, UVcontour);
      else
        err=_faces[i]->parametersSample(pos3D, UVcontour);
    
      if (err)
      {
        if (t==1)
          std::cout << "ERROR Face : " << i << " : cannot retrieve parametrization !" << std::endl;
        continue;
      }
      
      // Need to reorient trimmed surface (i.e. with holes)
      {
        std::vector<K_FLD::IntArray> cntLoops;
        std::set<E_Int> dummy;
        ContourSplitter<K_MESH::Edge, E_Int>::splitConnectivity(connectB, dummy, cntLoops);
        E_Int nb_loops = cntLoops.size();
        if (nb_loops > 1)
        {
          std::vector<E_Int> indices;
          K_SEARCH::BBox3D boxOuter, box;
          K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd(UVcontour);
          E_Int outer=0;
          cntLoops[0].uniqueVals(indices);
          boxOuter.compute(acrd, indices);
          for (size_t i=1; i < nb_loops; ++i)
          {
            cntLoops[i].uniqueVals(indices);
            box.compute(acrd, indices);
            if (!box.is_included(boxOuter))
            {
              outer=i;
              boxOuter=box;
            }
          }
          
          //now sort and reorient
          E_Int o;
          Vector_t<E_Int> sorted_nodes;
          E_Int E[2];
          for (size_t i=0; i < nb_loops; ++i)
          {
            sorted_nodes.clear();
            BARSplitter::getSortedNodes(cntLoops[i], sorted_nodes);
            
            cntLoops[i].clear();
            size_t sz(sorted_nodes.size());
            for (size_t j=0; j < sz; ++j)
            {
              E[0]=sorted_nodes[j];
              E[1]=sorted_nodes[(j+1)%sz];
              
              cntLoops[i].pushBack(&E[0], &E[0]+2);
            }

            __computeOrient(UVcontour, cntLoops[i], o);
            
            if ( (i==outer && o==-1) || (i != outer && o==1) )
              for (size_t j=0; j < sz; ++j)
                std::swap(cntLoops[i](0,j), cntLoops[i](1,j));
          }
          
          //concatenate back to connectB
          connectB.clear();
          for (size_t i=0; i < nb_loops; ++i)
            connectB.pushBack(cntLoops[i]);
        }
      }
      
#ifdef DEBUG_CAD_READER
      if (i==faulty_id)
        MIO::write("connectBUV.mesh", UVcontour , connectB, "BAR");
      //std::cout << UVcontour << std::endl;
#endif
    
      //
      OCCSurface occ_surf(F._F);
      DELAUNAY::SurfaceMeshData<OCCSurface> data(UVcontour, pos3D, connectB, occ_surf);
    
#ifdef DEBUG_CAD_READER
    
      /*std::cout << "UV : " << UVcontour.rows() << "/"  << UVcontour.cols() << std::endl;
      std::cout << "pos3D : " << pos3D.rows() << "/"  << pos3D.cols() << std::endl;
      std::cout << "connectB : " << connectB.rows() << "/"  << connectB.cols() << std::endl;*/
    
#ifdef DEBUG_MESHER
      if (i==faulty_id)
        mesher.dbg_flag=true;
      else
        mesher.dbg_flag=false;
#endif
#endif
    
      err = mesher.run (data);
      if (err || (data.connectM.cols() == 0))
      {
        if (t==0)
          continue;
        if (err)
          std::cout << "ERROR Face : " << i << " : Geometric Mesher failed." << std::endl;
        else
          std::cout << "ERROR Face : " << i << " : Cannot retrieve parametrization (OCC Limitation) !" << std::endl;    
        continue;
      }
    
#ifdef DEBUG_CAD_READER
      /*{
      std::ostringstream o;
      o << "surfaceUV_" << i << ".mesh";
      meshIO::write(o.str().c_str(), data.pos, data.connectM);
      }*/
    {
      std::ostringstream o;
      o << "surface3D_" << i << ".mesh";
      MIO::write(o.str().c_str(), data.pos3D, data.connectM, "TRI");
      }
#endif
    
      E_Int shift = crds[F._parent].cols();
      crds[F._parent].pushBack(data.pos3D);
      data.connectM.shift(shift);
    
      connectMs[F._parent].pushBack(data.connectM);
      
      if (!err) // done
        break;
      }
    }

    //Final cleaning and compacting
    {
      std::vector<E_Int> nids;
      for (E_Int i=0; i <= max_solid_id; ++i)
      {
        if (crds[i].cols()==0)
          continue;
        K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(crds[i]);
        ::merge(crdA, E_EPSILON, nids);
        K_FLD::IntArray::changeIndices(connectMs[i], nids);
        nids.clear();
        K_CONNECT::MeshTool::compact_to_mesh(crds[i], connectMs[i], nids);
      }
    }
  
#ifdef DEBUG_CAD_READER
  {
  K_FLD::IntArray tmp;
  K_FLD::FloatArray crd;
  for (E_Int i=0; i <= max_solid_id; ++i)
  {
    if (connectMs[i].cols()*crds[i].cols())
      tmp.pushBack(connectMs[i]);crd.pushBack(crds[i]);
  }
  MIO::write("surfaceALL.mesh", crd, tmp, "TRI");
  }
  
#endif
  return 0;
}

void K_OCC::CADviaOCC::__computeOrient(const K_FLD::FloatArray crd2D, const K_FLD::IntArray& cnt, E_Int&o)
{
  o=0;
  E_Float Z=0., z=0.;
  for (size_t i=0; i < cnt.cols(); ++i)
  {
    K_FUNC::crossProduct<2>(crd2D.col(cnt(0,i)), crd2D.col(cnt(1,i)), &z);
    Z +=z;
  }
  
  o = (Z > 0.) ? 1 : -1;
}


///
void K_OCC::CADviaOCC::__traverse_face_edges(const TopoDS_Face& F, TopExp_Explorer& edge_expl, std::vector<E_Int>& edges)
{
  for (edge_expl.Init(F, TopAbs_EDGE); edge_expl.More(); edge_expl.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(edge_expl.Current());
     
    if (BRep_Tool::Degenerated (E))
      continue;

    // Get edge id in the flat list
    E_Int edge_idx = _edges.FindIndex(E);
    if (edge_idx == 0) //doesn' exist so add it (due to surface of revolution process)
    {
      _edges.Add(E);
      edge_idx = _edges.FindIndex(E);
    }

#ifdef DEBUG_CAD_READER
    //assert (E.IsSame (_edges(edge_idx)));
#endif

    // Take orientation into account
    if (E.Orientation() != _edges(edge_idx).Orientation())
      edge_idx = -edge_idx;
        
    edges.push_back(edge_idx);
  }
}

//
E_Int K_OCC::CADviaOCC::__build_graph(const TopoDS_Shape& occ_shape, std::vector<OCCSurface*>& vFG)
{
  int err = 0;
  
  _surfs.Clear();
  _edges.Clear();
  
  // Store all the faces and edges in a flat list.
  
  TopExp::MapShapes(occ_shape, TopAbs_FACE, _surfs);
  E_Int nb_faces = _surfs.Extent();
  
#ifdef DEBUG_CAD_READER
  std::cout << "nb of faces : " <<  nb_faces << endl;
#endif
  
  TopExp::MapShapes(occ_shape, TopAbs_EDGE, _edges);
  E_Int nb_edges = _edges.Extent();
  
#ifdef DEBUG_CAD_READER
  std::cout << "nb of edges : " << nb_edges << endl;
#endif
  
  // Now build the graph : for each Face in _surfs associate edges ids in _edges and stamp the solid id.
  
  vFG.resize(nb_faces+1, 0);
  
  TopExp_Explorer top_expl, edge_expl;
  TopTools_IndexedMapOfShape sol_surfs;
  E_Int nb_surfs, nb_solids(0);
  TopTools_IndexedMapOfShape toto;
  // Traverse the solids and stamp the faces as belonging to the current solid
  // And store ORIENTED edges ids. //fixme : orinetation doesn't seem to work..
  for (top_expl.Init(occ_shape, TopAbs_SOLID); top_expl.More(); top_expl.Next(), ++nb_solids)
  {
    sol_surfs.Clear();
    
    TopExp::MapShapes (top_expl.Current(), TopAbs_FACE, sol_surfs);
    
#ifdef DEBUG_CAD_READER
    E_Int nb_edges2=0;
#endif
    
    nb_surfs = sol_surfs.Extent();
    
    for (E_Int i = 1; i <= nb_surfs; ++i)
    {
      E_Int idx = (E_Int)_surfs.FindIndex(sol_surfs(i));
      const TopoDS_Face& F=TopoDS::Face(sol_surfs(i));
            
      if (BRep_Tool::Surface(F)->IsUClosed() || BRep_Tool::Surface(F)->IsVClosed()) // Surface of Revolution
      {
      // Anomaly #3848: Lecture IGES finocyl
      // 0. The bounds of some iges TRIMMED surfaces are not retrieved (OCC issue ? ) so the following split will cut a basis surface (mostly cylinder !)
      // 1. The follwoing has been commented because it's not necessary and is probably not robust
      //     as the finocyl.iges case fails with it : only the first one is output upon exit (don't know why)

      //  // fixme : hugly workaround to get a TopoDS_Face from a Geom_SurfaceOfRevolution
      //  TopoDS_Shape shor = ShapeCustom::ConvertToRevolution(_surfs(i));
      //  TopExp::MapShapes(shor, TopAbs_SHELL, toto);
      //  E_Int nb_faces = toto.Extent(); 
      //  const TopoDS_Face& Fb = TopoDS::Shell(toto(1)); 
       
        err = __split_surface_of_revolution(F, vFG, nb_solids);
      }
      else
      {
        /*vFG[idx] = new OCCSurface(F);
      // Traverse the edges
      __traverse_face_edges(F, edge_expl, vFG[idx]->_edges);*/

      vFG[idx] = new OCCSurface(F, nb_solids);
      // Traverse the edges
      __traverse_face_edges(F, edge_expl, vFG[idx]->_edges);
      
      Handle(Geom_Surface)  S = BRep_Tool::Surface(F); 
      S->Bounds(vFG[idx]->_U0,vFG[idx]->_U1, vFG[idx]->_V0, vFG[idx]->_V1);
        
#ifdef DEBUG_CAD_READER
      nb_edges2 +=vFG[idx]->_edges.size();
#endif
      }
    }
  }

#ifdef DEBUG_CAD_READER
    E_Int nb_edges2=0;
#endif
  // Traverse the orphan faces (not belonging to any solid)
  // And store ORIENTED edges ids //fixme : orinetation doesn't seem to work..
  for (top_expl.Init(occ_shape, TopAbs_FACE, TopAbs_SOLID); top_expl.More(); top_expl.Next())
  {
    E_Int idx = (E_Int)_surfs.FindIndex(top_expl.Current());   
    const TopoDS_Face& F=TopoDS::Face(_surfs(idx));
        
    if (BRep_Tool::Surface(F)->IsUClosed() || BRep_Tool::Surface(F)->IsVClosed()) // Surface of Revolution
    {
      // Anomaly #3848: Lecture IGES finocyl
      // 0. The bounds of some iges TRIMMED surfaces are not retrieved (OCC issue ? ) so the following split will cut a basis surface (mostly cylinder !)
      // 1. The follwoing has been commented because it's not necessary and is probably not robust
      //     as the finocyl.iges case fails with it : only the first one is output upon exit (don't know why)

      //  // fixme : hugly workaround to get a TopoDS_Face from a Geom_SurfaceOfRevolution
      //  TopoDS_Shape shor = ShapeCustom::ConvertToRevolution(_surfs(idx));
      //  TopExp::MapShapes(shor, TopAbs_SHELL, toto);
      //  E_Int nb_faces = toto.Extent(); 
      //  const TopoDS_Face& Fb = TopoDS::Shell(toto(1)); 
 
      err = __split_surface_of_revolution(F, vFG);
    }
    else
    {
      /*vFG[idx] = new OCCSurface(F);
      // Traverse the edges
      __traverse_face_edges(F, edge_expl, vFG[idx]->_edges);*/
        
      Handle(Geom_Surface)  S = BRep_Tool::Surface(F); 
      Standard_Real U[3], V[3]; // store {Umin, 0.5*Umax, Umax} if periodic, {Umin, Umax, Umax} otherwise. Idem for V. 
        
      vFG[idx] = new OCCSurface(F, nb_solids);
        
      // Traverse the edges
      __traverse_face_edges(F, edge_expl, vFG[idx]->_edges);
      
      S->Bounds(vFG[idx]->_U0,vFG[idx]->_U1, vFG[idx]->_V0, vFG[idx]->_V1);
        
#ifdef DEBUG_CAD_READER
      nb_edges2 +=vFG[idx]->_edges.size();
#endif
    }
  }

  return err;
}

///
E_Int K_OCC::CADviaOCC::__split_surface_of_revolution
(const TopoDS_Face& f, std::vector<OCCSurface*>& vFG, E_Int nb_solid)
{
//#ifdef DEBUG_CAD_READER
//      std::cout << " Surface " << idx << " is revolving, it is not supported yet." << std::endl;
//#endif
//      continue;
  
  Handle(Geom_Surface)  S = BRep_Tool::Surface(f); 
  
  Standard_Real U[3], V[3]; // store {Umin, 0.5*Umax, Umax} if periodic, {Umin, Umax, Umax} otherwise. Idem for V. 
  S->Bounds(U[0],U[1], V[0], V[1]);
  
  U[2]=std::max(U[0], U[1]); 
  V[2]=std::max(V[0], V[1]);
  U[0]=std::min(U[0], U[1]);
  V[0]=std::min(V[0], V[1]);
  
  E_Int umax(2), vmax(2);
  if (S->IsUClosed())
    U[1]=0.5*U[2];
  else
    umax=1;
  
  if (S->IsVClosed())
    V[1]=0.5*V[2];
  else
    vmax=1;
  
  // Up to 4 bits for a sphere, 2 bits for a classical surface of revolution
  for (size_t v=0; v < vmax; ++v)
  {
    for (size_t u=0; u < umax; ++u)
    {
      BRepBuilderAPI_MakeFace mf1(S, U[u], U[u+1], V[v], V[v+1], _merge_tol);
    
      TopoDS_Face* F1 = new TopoDS_Face(mf1.Face());
      _surfs.Add(*F1);
    
      E_Int idx = (E_Int)_surfs.FindIndex(*F1);  
      if (idx >= vFG.size()) vFG.resize(idx+1, 0);
        vFG[idx] = new OCCSurface(*F1);
        
      vFG[idx]->_U0=U[u];
      vFG[idx]->_U1=U[u+1];
      vFG[idx]->_V0=V[v];
      vFG[idx]->_V1=V[v+1];
          
      // Traverse the edges
      TopExp_Explorer edge_expl;
      __traverse_face_edges(*F1, edge_expl, vFG[idx]->_edges);
    }
  }
  
  return 0;
}

E_Int K_OCC::CADviaOCC::__mesh_edge(const TopoDS_Edge& E, E_Int& nb_points, K_FLD::FloatArray& coords, K_FLD::IntArray& connectE)
{
  connectE.clear();
  
  if (BRep_Tool::Degenerated (E))      // Exit if the edge is degenerated.
    return 1;
  
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geom_adap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
       
  // fixme ? : uniform because cannot identify 2 edges with different orientations currently.
  // consequently, a non uniform meshing will give 2 different result depending on which extremity is used to start
  // it would then lead to a mismatch between shared edge discretizations.
  Standard_Real u0=geom_adap.FirstParameter();
  Standard_Real u1=geom_adap.LastParameter();
  GCPnts_UniformAbscissa unif_abs (geom_adap, nb_points, u0, u1);
  if (!unif_abs.IsDone())
    return 1;
   
  nb_points = unif_abs.NbPoints();// just in case the number of constructed points is different from what was asked.
    
  gp_Pnt Pt;
  E_Float P[3];
  E_Int Ei[2], sz(coords.cols());
   
  // Insert new points
  for (Standard_Integer i = 1; i <= nb_points; ++i) //in case NbPoints() != nb_points)
  {
    C0.D0 (unif_abs/*unif_defl*/.Parameter(i), Pt);
    P[0]=Pt.X();P[1]=Pt.Y();P[2]=Pt.Z();
      
    coords.pushBack(P, P+3);
  }
   
  // Insert new edges
  for (Standard_Integer i=0; i < nb_points-1; ++i)
  {
    Ei[0]=sz+i; Ei[1]=sz+i+1;
    connectE.pushBack(Ei, Ei+2);
  }
  return 0;
}




