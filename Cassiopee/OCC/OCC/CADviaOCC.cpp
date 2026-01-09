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
//Author : Sam Landier (sam.landier@onera.fr)

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
#include "String/kstring.h"
/*
#include "TopoDS_Edge.hxx"
// Curve and Mesh
#include "GeomAdaptor_Curve.hxx"
#include "Geom_Line.hxx"
#include "gp_Lin.hxx"
*/

#ifdef DEBUG_CAD_READER
#include "Nuga/include/medit.hxx"
#include <sstream>
#include <iostream>
#endif

K_OCC::CADviaOCC::CADviaOCC()
{
}


K_OCC::CADviaOCC::~CADviaOCC()
{
  for (size_t i=1; i < _faces.size(); ++i) delete _faces[i];
}

E_Int import_iges(const char* fname, TopoDS_Shape& sh)
{    
  // Read the file
  IGESControl_Reader reader;
  reader.ReadFile(fname);
    	
  // Transfer CAD faces (only) into a OCC list
  Handle(TColStd_HSequenceOfTransient) occ_list = reader.GiveList("iges-faces");
    
  if (occ_list.IsNull()) return 1;

#ifdef DEBUG_CAD_READER
  Standard_Integer nb_transfered_faces =
#endif 
  reader.TransferList(occ_list);
  
#ifdef DEBUG_CAD_READER
  Standard_Integer nb_cad_faces = occ_list->Length();
  std::cout << "IGES Faces: " << nb_cad_faces << "   Transferred:" << nb_transfered_faces << std::endl;
#endif

  sh = reader.OneShape();
  return sh.IsNull();
}


int import_step(const char* fname, TopoDS_Shape& sh)
{
  // Read the file
  STEPControl_Reader reader;
  reader.ReadFile(fname);
  reader.TransferRoots();
  sh = reader.OneShape();
  return sh.IsNull();
}

// Lit le fichier CAD et retourne les entites openscascade
E_Int K_OCC::CADviaOCC::import_cad(const char* fname, const char* format, E_Float h, E_Float chordal_err, E_Float gr /*groqth ratio*/)
{
  _chordal_err = chordal_err;
  _h = h;
  _gr = gr;
  E_Int err(1);
  if (::strcmp(format, "fmt_iges")==0) err = import_iges(fname, _occ_shape);
  else if (::strcmp(format, "fmt_step")==0) err = import_step(fname, _occ_shape);
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
  
  _Lmin = NUGA::FLOAT_MAX;
  _Lmean = _Lmax = 0.;
  int nb_valid_edge{ 0 };
  
  for (E_Int i=1; i <= nb_edges; ++i)
  {
    int er = __h_sizing(TopoDS::Edge(_edges(i)), Ls[i]);

    if (er == 1) continue; // i-th edge is degen

    _Lmin = std::min(Ls[i], _Lmin);
    _Lmax = std::max(Ls[i], _Lmax);
    _Lmean += Ls[i];
    ++nb_valid_edge;
  }
  
  _Lmean /= nb_valid_edge;
  
#ifdef DEBUG_CAD_READER
  std::cout << "Lmin/Lmean/Lmax : " << _Lmin << "/" << _Lmean << "/" << _Lmax << std::endl;
#endif
  
  if (_h <= 0.) // undefined or badly defined
  {
    _h = _Lmean/10.;
  }
  
  std::cout << "OCC: size h: " << _h << std::endl;
  
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
    _chordal_err = 0.02;

  for (E_Int i=1; i <= nb_edges; ++i)
  {
    __chord_sizing(TopoDS::Edge(_edges(i)), _chordal_err, Ns[i]);
  }

  std::cout << "OCC: chordal_error h: " << _chordal_err << std::endl;
  
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

// Retourne le nbre de points pour recuperer le bon chordal_err (fait un max)
E_Int K_OCC::CADviaOCC::__chord_sizing(const TopoDS_Edge& E, E_Float chordal_err, E_Int& nb_points)
{  
  if (BRep_Tool::Degenerated (E))      // Exit if the edge is degenerated.
    return 1;
  
  BRepAdaptor_Curve C0(E);
  GeomAdaptor_Curve geom_adap(C0.Curve()); // Geometric Interface <=> access to discretizations tool
  E_Float u0 = geom_adap.FirstParameter();
  E_Float u1 = geom_adap.LastParameter();
  
  //E_Float L = (E_Float) GCPnts_AbscissaPoint::Length(geom_adap, u0, u1);
  
  E_Int nb_pts_defl = 0;
  __eval_nb_points(C0, u0, u1, chordal_err, nb_pts_defl);
  nb_pts_defl += 1; //give the number of split => +1 to have the nb of points
  nb_pts_defl = std::max(nb_pts_defl, E_Int(3)); // at least 3 points
  nb_points = std::max(nb_pts_defl, nb_points); // max des nbre de pts
  return 0;
}

// Calcul le nbre de pts pour avoir la bonne erreur de corde sur la courbe C
E_Int K_OCC::CADviaOCC::__eval_nb_points(const BRepAdaptor_Curve& C, E_Float u0, E_Float u1, E_Float chordal_err, E_Int& nb_points)
{

  GeomAdaptor_Curve geom_adap(C.Curve()); // Geometric Interface <=> access to discretizations tool
  E_Float L = (E_Float)GCPnts_AbscissaPoint::Length(geom_adap, u0, u1);

  E_Float dm;
  __eval_chordal_error(C, u0, u1, dm);
  if (dm < chordal_err * L)
  {
    nb_points = 1; return 0;
  }
  E_Int ns = 1; // nbre de splits
  E_Float du = u1-u0;
  
  while (ns < 50)
  {
    E_Float cm1 = 0.;
    for (E_Int i = 0; i <= ns; i++)
    {
      __eval_chordal_error(C, u0+i*du/(ns+1), u0+(i+1)*du/(ns+1), dm);
      E_Float Li = (E_Float)GCPnts_AbscissaPoint::Length(geom_adap, u0 + i * du / (ns + 1), u0 + (i + 1)*du / (ns + 1));
      cm1 = std::max(cm1, dm / Li);
    }
    if (cm1 < chordal_err) { nb_points = ns; return 0; }
    ns += 1;
  }
  nb_points = ns;
  return 0;
}

// Fonction recursive : parfois boucle a l'infini
E_Int K_OCC::CADviaOCC::__eval_nb_points2(const BRepAdaptor_Curve& C, E_Float u0, E_Float u1, E_Float dmax, E_Int& nb_points)
{
  E_Float dm;
  __eval_chordal_error(C, u0, u1, dm);
  if (dm <= dmax)
  {
    ++nb_points; return 0;
  }
  
  __eval_nb_points(C, u0, 0.5*(u1+u0), dmax, nb_points);
  __eval_nb_points(C, 0.5*(u1+u0), u1, dmax, nb_points);
  
  return 0;
}

E_Int K_OCC::CADviaOCC::__eval_chordal_error(const BRepAdaptor_Curve& C, E_Float u0, E_Float u1, E_Float& dmax)
{
  gp_Pnt pu0;
  C.D0 (u0, pu0);
  gp_Pnt pu1;
  C.D0 (u1, pu1);
    
  gp_Pnt pu, P0, P1;
  E_Float Pu[3], Pm[3];
  
  C.D0 (u0, P0);
  C.D0 (u1, P1);
  
  dmax = -1;
  
  // 4 samples
  for (size_t n=0; n < 3; ++n)
  {
    E_Float u = u0 + 0.25 * (n+1) * (u1-u0);
    
    C.D0 (u, pu);
    Pu[0] = pu.X(); Pu[1] = pu.Y(); Pu[2] = pu.Z();
    
    Pm[0] = P0.X() + 0.25 * (n+1) * (P1.X() - P0.X());
    Pm[1] = P0.Y() + 0.25 * (n+1) * (P1.Y() - P0.Y());
    Pm[2] = P0.Z() + 0.25 * (n+1) * (P1.Z() - P0.Z());
    
    E_Float dm = ::sqrt(NUGA::sqrDistance(Pm, Pu, 3));
    
    dmax = (dm > dmax) ? dm : dmax;
  }
  return 0;
}

// maille les edges (uniforme)
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
   
    const OCCSurface& F = *_faces[i];
    
    E_Int nb_edges = F._edges.size();
    
    for (E_Int j=0; j < nb_edges; ++j)
    {
      E_Int id = ::abs(F._edges[j]);
      
      if (connectEs[id].cols()) continue; //already meshed.
    
#ifdef DEBUG_CAD_READER
    std::cout << "Edge :  " << id <<  ". Nb points : " << Ns[id] << std::endl;
#endif
      __mesh_edge(TopoDS::Edge(_edges(id)), Ns[id], coords, connectEs[id]);
      
#ifdef DEBUG_CAD_READER
    if (id==41 || id==101)
      medith::write("Ei.mesh", coords, connectEs[id], "BAR");
#endif
    }
  }
    
#ifdef DEBUG_CAD_READER
  K_FLD::IntArray tmp;
  for (size_t i=0; i <connectEs.size(); ++i)
    tmp.pushBack(connectEs[i]);
  medith::write("wire.mesh", coords, tmp, "BAR");
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
(K_FLD::FloatArray& coords, const std::vector<K_FLD::IntArray>& connectEs, 
  std::vector<K_FLD::IntArray>& connectBs, E_Float merge_tol)
{
  E_Int nb_faces(_surfs.Extent());
  
  if (!nb_faces) return 0;
  
  connectBs.resize(nb_faces+1);
  
  K_FLD::IntArray tmp;
  K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(coords);
  std::vector<E_Int> end_nodes;
  E_Float tol2 = -1.;
  
#ifdef DEBUG_CAD_READER
  E_Int faulty_id=69;
#endif
  
  for (E_Int i=1; i <= nb_faces; ++i)
  {
    //std::cout << " cad face nb : " << i << std::endl;
    
    const OCCSurface& F = *_faces[i];
    
    std::vector<E_Int> unods;
        
    // prepare contour  
    E_Int nb_edges = F._edges.size();
    //
    for (E_Int j=0; j < nb_edges; ++j)
    {
      E_Int id = F._edges[j];
      
      unods.clear();
      
      if (id > 0)
      {
        if (connectEs[id].cols()==0) continue;
        connectBs[i].pushBack(connectEs[id]);
        connectEs[id].uniqueVals(unods); 
      }
      else
      {
        tmp=connectEs[-id];
        if (tmp.cols()==0) continue;
        tmp.uniqueVals(unods);
        for (E_Int k=0; k < tmp.cols(); ++k) std::swap(tmp(0,k), tmp(1,k));
        connectBs[i].pushBack(tmp);
      }
    }

#ifdef DEBUG_CAD_READER
    if (faulty_id==i)
      medith::write("connectBi.mesh", coords, connectBs[i], "BAR");
#endif
    
    // clean it
    while (__clean(crdA, connectBs[i], tol2))
    ;
    
#ifdef DEBUG_CAD_READER
    if (faulty_id==i)
      medith::write("connectBc.mesh", coords, connectBs[i], "BAR");
    //std::ostringstream o;
    //o << "clean_loop_" << i << ".mesh";
    //meshIO::write(o.str().c_str(), coords, connectBs[i]);
#endif
  }
  
  // merge tol
  _merge_tol = ::sqrt(tol2); 
  _merge_tol = std::max(_merge_tol, 1.e-4*_Lmean);
  if (merge_tol > 0.) _merge_tol = merge_tol;
  printf("merge tol = %g\n", _merge_tol);
  
  // Global pass to join the loops.  
  std::vector<E_Int> nids;
  E_Int nb_merges = ::merge(crdA, _merge_tol, nids);
  
  if (!nb_merges) return 0;

  for (E_Int i=1; i <= nb_faces; ++i)
  {
    K_FLD::IntArray::changeIndices(connectBs[i], nids);
    __remove_degenerated(connectBs[i]);
    
#ifdef DEBUG_CAD_READER
    //medith::write("connectBf.mesh",coords , connectBs[i]);
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

// Ressort aussi la tolerance
E_Int K_OCC::CADviaOCC::__clean(const K_FLD::ArrayAccessor<K_FLD::FloatArray>& crdA, K_FLD::IntArray& connectB, E_Float& tol2)
{ 
  //
  _end_nodes.clear();
  _enodes.clear();
  
  // calcul la valence
  for (E_Int i=0; i < connectB.cols(); ++i)
  {
    if (_enodes.find(connectB(0,i)) == _enodes.end())
      _enodes[connectB(0,i)] = 1;
    else
      _enodes[connectB(0,i)] += 1;
    if (_enodes.find(connectB(1,i)) == _enodes.end())
      _enodes[connectB(1,i)] = 1;
    else
      _enodes[connectB(1,i)] += 1;
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
  
  K_FLD::IntArray::changeIndices(connectB, nids);
  
  return 1; 
}

// Parametrise les edges et appelle le mailleur par face
E_Int K_OCC::CADviaOCC::mesh_faces
(const K_FLD::FloatArray& coords, const std::vector<K_FLD::IntArray>& connectBs, std::vector<K_FLD::FloatArray>& crds, std::vector<K_FLD::IntArray>& connectMs, bool aniso, bool do_join)
{
  E_Int nb_faces{_surfs.Extent()};
  
  if (!nb_faces) return 0;

  std::vector<K_FLD::FloatArray> crds1(nb_faces);
  std::vector<K_FLD::IntArray> connectMs1(nb_faces);
  
  std::vector<E_Int> nodes, nids;
  K_FLD::FloatArray UVcontour, pos3D;
  
  DELAUNAY::SurfaceMesher<OCCSurface> mesher;
  
#ifdef DEBUG_CAD_READER
  E_Int faulty_id = 3;
#endif
  E_Int t;

#ifndef DEBUG_CAD_READER // omp still disabled because not good perf (mem concurrency?) BUT now result is multithread-independant
//#pragma omp parallel default(none) shared(coords, connectBs, aniso, nb_faces, crds1, connectMs1) private (t, nodes, nids, pos3D, UVcontour, mesher)
#endif
  for (E_Int i=1; i <= nb_faces; ++i)
  {
    std::cout << "Processing face: " << i << " / "<< nb_faces << std::endl;
#ifdef DEBUG_CAD_READER
    std::cout << "Processing face: " << i << " / "<< nb_faces << std::endl;
#endif

    const OCCSurface& F = *_faces[i];
    
    K_FLD::IntArray connectB = connectBs[i];
    
    if (connectB.cols() == 0)
    {
#ifdef DEBUG_CAD_READER
      std::cout << "ERROR Face : " << i << " : empty discretized contour!" << std::endl;
#endif
      continue;
    }
    
#ifdef DEBUG_CAD_READER
    if (i==faulty_id)
      medith::write("connectB.mesh", coords, connectB);
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
      medith::write("connectBcompacted.mesh", pos3D, connectB, "BAR");
#endif
    
#ifdef DEBUG_CAD_READER
    //if (i == faulty_id)
    {
      K_FLD::FloatArray surfc;
      K_FLD::IntArray con;
      _faces[i]->discretize(surfc, con, 30, 30);
      std::ostringstream o;
      o << "discretized_surf_" << i;
      //medith::write(o.str().c_str(), surfc, con, "QUAD");
    }
#endif
      
    // surface of revolution => duplicate, reverse and separate seams
    //bool is_of_revolution = ((E_Int)nodes.size() != connectB.cols());
    bool is_of_revolution = (nodes.size() != (size_t)connectB.cols());
    
    std::map<E_Int, std::pair<E_Int, E_Int> > seam_nodes;
    
    if (is_of_revolution)
    {
      _faces[i]->_normalize_domain = false; // fixme : currently normalizing not working with revol surfaces.
      __split_surface_of_revolution(_faces[i], connectB, pos3D, seam_nodes);
      
    }

    E_Int nb_loops = 1;
    std::vector<K_FLD::IntArray> cntLoops;
    {
      std::set<E_Int> dummy;
      ContourSplitter<K_MESH::Edge, E_Int>::splitConnectivity(connectB, dummy, cntLoops);
      nb_loops = cntLoops.size();
    }
    
    // Up to 2 tries : first by asking OCC for params, Second by "hand" (sampling)
    E_Int err = 0;
    for (t=0; t<2; ++t) // supp. la parametrisation discrete
    {
      if (t==0)
        err = _faces[i]->parameters(pos3D, connectB, UVcontour);
      else if (nb_loops == 1)
        err = _faces[i]->parametersSample(pos3D, UVcontour);
      else
      {
        // todo : try to mesh in the contour mean plane
        err = 1;
      }
      
      if (!err)
      {
        // check if there are spikes in the contour :  == angular node equal to 0 == overlapping edges
        err = __check_for_spikes(cntLoops, UVcontour);
      }

      if (!err) // Need to reorient holed surface.
      {
        err = __reorient_holed_surface(cntLoops, UVcontour);
        //concatenate back to connectB
        connectB.clear();
        for (E_Int c = 0; c< nb_loops; ++c) connectB.pushBack(cntLoops[c]);
      }
      
      if (err)
      {
#ifdef DEBUG_CAD_READER
        if (t==1)
          std::cout << "ERROR Face : " << i << " : cannot retrieve parametrization !" << std::endl;
#endif
        continue;
      }

#ifdef DEBUG_CAD_READER
    //if (i == faulty_id)
    {
      E_Int nj = 50;//connectB.cols() /  2;
      E_Int ni = 100;//2 * nj;
      K_FLD::FloatArray surfc;
      K_FLD::IntArray con;
      _faces[i]->discretize(surfc, con, ni, nj);
      
      std::ostringstream o;
      o << "patch_" << i ;
      crds1[i-1] = surfc;
      connectMs1[i-1] = con;
      medith::write(o.str().c_str(), crds1[i-1], connectMs1[i-1], "QUAD");
    }
#endif
      
#ifdef DEBUG_CAD_READER
      if (/*i==faulty_id&&*/ t==0)
        medith::write("connectBUV1.mesh", UVcontour, connectB, "BAR");
      else if (/*i==faulty_id&&*/ t==1)
        medith::write("connectBUV2.mesh", UVcontour, connectB, "BAR");
      /*else if (t == 2)
        medith::write("connectBUV3.mesh", UVcontour, connectB, "BAR");*/
      //std::cout << UVcontour << std::endl;
#endif
    
      OCCSurface occ_surf(F);
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

      auto& mode = mesher.mode;
      
      mode.chordal_error = _chordal_err;
      mode.symmetrize = false; //reset 

      if (aniso) mode.metric_mode = mode.ANISO;

#ifndef DEBUG_CAD_READER
      mode.silent_errors = true;
#endif
      if (_gr <= 0.) // unspecified == OLD MODE
      {
        mode.hmax = _h;
        mode.growth_ratio = 0.;
      }
      else
      {
        mode.growth_ratio = std::max(_gr, 1.); //fixme : values in [1.; 1.25] might cause errors. Anyway do not allow bellow 1. since not useful
        E_Int MIN_NB = 50;
        if (_gr >= 1. && (connectB.cols() > (MIN_NB * 4))) // check if growth ratio is applicable to this patch, i.e. it is symmetrizable
        {
          // Count the minimum number of edges on a boundary of the param space per direction (U,V)
          E_Float nu(0), nv(0);
          for (E_Int i=0; i < connectB.cols(); ++i)
          {
            E_Float ds[2];
            const E_Int& Ni = connectB(0,i);
            const E_Int& Nj = connectB(1,i);
            NUGA::diff<2>(UVcontour.col(Ni), UVcontour.col(Nj), ds);
            NUGA::normalize<2>(ds);
            
            if (::fabs(ds[0]) < ::fabs(ds[1])) ++nv;
            else ++nu;
          }
          
          E_Int nmin = std::min(nu, nv) / 2; // nu and nv are accumulation for "2" sides
          E_Int nmax = std::max(nu, nv) / 2; // nu and nv are accumulation for "2" sides
          
          if (nmin > 0.5*nmax) // roughly iso domain
          {
            mode.symmetrize = true;
            mode.nb_smooth_iter = 2;
          }
          if (_gr == 1.) mode.hmax = _h;
        }
      }
      printf("selected sym=%d grading=%g hmax=%g hmin=%g smooth=" SF_D_ "\n", 
        mode.symmetrize, mode.growth_ratio, mode.hmax, mode.hmin, mode.nb_smooth_iter);

      mesher.clear();
      mesher.seed_random(1);
      err = mesher.run(data);

      if (err || (data.connectM.cols() == 0))
      {
        if (t==0) continue;

#ifdef DEBUG_CAD_READER
        if (err)
          std::cout << "ERROR Face : " << i << " : Geometric Mesher failed." << std::endl;
        else
          std::cout << "ERROR Face : " << i << " : Cannot retrieve parametrization (OCC Limitation) !" << std::endl; 
#endif   
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
    
#ifdef DEBUG_CAD_READER
      {
      std::ostringstream o;
      o << "surfaceUV_" << i << ".mesh";
      medith::write(o.str().c_str(), *data.pos, data.connectM);
      }
      {
      std::ostringstream o;
      o << "surface3D_" << i << ".mesh";
      medith::write(o.str().c_str(), data.pos3D, data.connectM, "TRI");
      }
#endif
    
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

        /*K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(crds[i]);
        ::merge(crdA, _merge_tol, nids);
        K_FLD::IntArray::changeIndices(connectMs[i], nids);
        nids.clear();
        NUGA::MeshTool::compact_to_mesh(crds[i], connectMs[i], nids);*/
      }
    }
    else
    {
      connectMs = connectMs1;
      crds = crds1;
    }
  
#ifdef DEBUG_CAD_READER
  /*{ manque un shift dans avant d'ajouter connectMs[i]
    K_FLD::IntArray tmp;
    K_FLD::FloatArray crd;
    for (E_Int i=0; i < connectMs.size(); ++i)
    {
      if (connectMs[i].cols()*crds[i].cols())
        tmp.pushBack(connectMs[i]);crd.pushBack(crds[i]);
    }
    medith::write("surfaceALL.mesh", crd, tmp, "TRI");
  }*/
#endif
  return 0;
}

void K_OCC::CADviaOCC::__computeOrient(const K_FLD::FloatArray crd2D, const K_FLD::IntArray& cnt, E_Int&o)
{
  o=0;
  E_Float Z=0., z=0.;
  for (E_Int i=0; i < cnt.cols(); ++i)
  {
    NUGA::crossProduct<2>(crd2D.col(cnt(0,i)), crd2D.col(cnt(1,i)), &z);
    Z += z;
  }
  
  o = (Z > 0.) ? 1 : -1;
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
  std::cout << "INFO: nb of faces: " <<  nb_faces << std::endl;
#endif
  
  TopExp::MapShapes(occ_shape, TopAbs_EDGE, _edges);
  
#ifdef DEBUG_CAD_READER
  E_Int nb_edges = _edges.Extent();
  std::cout << "INFO: nb of edges: " << nb_edges << std::endl;
#endif
  
  // Now build the graph: for each Face in _surfs associate edges ids in _edges and stamp the solid id.
  vFG.resize(nb_faces+1, 0);
  
  TopExp_Explorer top_expl;
  TopTools_IndexedMapOfShape sol_surfs;
  E_Int nb_surfs, nb_solids(0);
  TopTools_IndexedMapOfShape toto;
  // Traverse the solids and stamp the faces as belonging to the current solid
  // And store ORIENTED edges ids. //fixme : orientation doesn't seem to work..
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

      vFG[idx] = new OCCSurface(F, _edges, nb_solids);
  
#ifdef DEBUG_CAD_READER
      nb_edges2 += vFG[idx]->_edges.size();
#endif
    }
  }

#ifdef DEBUG_CAD_READER
    E_Int nb_edges2=0;
#endif

  // Traverse the orphan faces (not belonging to any solid)
  // And store ORIENTED edges ids //fixme : orientation doesn't seem to work..
  for (top_expl.Init(occ_shape, TopAbs_FACE, TopAbs_SOLID); top_expl.More(); top_expl.Next())
  {
    E_Int idx = (E_Int)_surfs.FindIndex(top_expl.Current());   
    const TopoDS_Face& F = TopoDS::Face(_surfs(idx));

    vFG[idx] = new OCCSurface(F, _edges, nb_solids);

#ifdef DEBUG_CAD_READER
    nb_edges2 +=vFG[idx]->_edges.size();
#endif
  }

  return err;
}

///
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
  Standard_Real u0 = geom_adap.FirstParameter();
  Standard_Real u1 = geom_adap.LastParameter();
  GCPnts_UniformAbscissa unif_abs(geom_adap, int(nb_points), u0, u1);
  if (!unif_abs.IsDone()) return 1;
   
  nb_points = unif_abs.NbPoints(); // just in case the number of constructed points is different from what was asked.
    
  gp_Pnt Pt;
  E_Float P[3];
  E_Int Ei[2], sz(coords.cols());
   
  // Insert new points
  for (Standard_Integer i = 1; i <= nb_points; ++i) //in case NbPoints() != nb_points)
  {
    C0.D0(unif_abs/*unif_defl*/.Parameter(i), Pt);
    P[0]=Pt.X(); P[1]=Pt.Y(); P[2]=Pt.Z();
    coords.pushBack(P, P+3);
  }
   
  // Insert new edges
  for (Standard_Integer i=0; i < nb_points-1; ++i)
  {
    Ei[0] = sz+i; Ei[1] = sz+i+1;
    connectE.pushBack(Ei, Ei+2);
  }
  return 0;
}

E_Int K_OCC::CADviaOCC::__reorient_holed_surface(K_FLD::IntArray& connectB, const K_FLD::FloatArray& UVcontour)
{
  std::vector<K_FLD::IntArray> cntLoops;
  std::set<E_Int> dummy;
  ContourSplitter<K_MESH::Edge, E_Int>::splitConnectivity(connectB, dummy, cntLoops);
  E_Int nb_loops = cntLoops.size();

  if (nb_loops == 1) return 0;

  E_Int err =  K_OCC::CADviaOCC::__reorient_holed_surface(cntLoops, UVcontour);
  if (err) printf("Warning: reorient_holed_surface: fail to reorient surface.\n");

  //concatenate back to connectB
  connectB.clear();
  for (E_Int i = 0; i < nb_loops; ++i) connectB.pushBack(cntLoops[i]);

  return 0;
}

E_Int K_OCC::CADviaOCC::__reorient_holed_surface(std::vector<K_FLD::IntArray>& cntLoops, const K_FLD::FloatArray& UVcontour)
{
  E_Int nb_loops = cntLoops.size();

  if (nb_loops == 1) return 0;

  std::vector<E_Int> indices;
  K_SEARCH::BBox2D boxOuter, box;
  K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd(UVcontour);
  E_Int outer=0;
  
  cntLoops[0].uniqueVals(indices);
  boxOuter.compute(acrd, indices);
  
  for (E_Int i=1; i < nb_loops; ++i)
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
  for (E_Int i=0; i < nb_loops; ++i)
  {
    sorted_nodes.clear();
    int err = BARSplitter::getSortedNodes(cntLoops[i], sorted_nodes);
    // the following test is added to catch a getSortedNodes error. Not added inside it for efficiency (generally works fine).
    if (std::find(sorted_nodes.begin(), sorted_nodes.end(), E_IDX_NONE) != sorted_nodes.end()) err = 1;
    if (err) return err;

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

  return 0;
}

E_Int K_OCC::CADviaOCC::__check_for_spikes(const std::vector<K_FLD::IntArray>& cntLoops, const K_FLD::FloatArray& UVcontour)
{
  E_Int err = 0;
  
  Vector_t<E_Int> sorted_nodes;
  E_Int nb_loops = cntLoops.size();

  for (E_Int l = 0; (l < nb_loops) && !err; ++l)
  {
    sorted_nodes.clear();
    err = BARSplitter::getSortedNodes(cntLoops[l], sorted_nodes);
    int nnodes = sorted_nodes.size();
    // the following test is added to catch a getSortedNodes error. Not added inside it for efficiency (generally works fine).
    if (std::find(sorted_nodes.begin(), sorted_nodes.end(), E_IDX_NONE) != sorted_nodes.end()) err = 1;
    for (size_t n = 0; n < sorted_nodes.size() && !err; ++n)
    {
      // detect spikes

      int Nim1 = sorted_nodes[n];
      int Ni = sorted_nodes[(n + 1) % nnodes];
      int Nip1 = sorted_nodes[(n + 2) % nnodes];

      double Pim1Pi[] = { UVcontour(0, Ni) - UVcontour(0, Nim1) , UVcontour(1, Ni) - UVcontour(1, Nim1) };
      NUGA::normalize<2>(Pim1Pi);
      double Norm1[] = { -Pim1Pi[1], Pim1Pi[0] , 0.};  // {-b, a}

      double PiPip1[] = { UVcontour(0, Nip1) - UVcontour(0, Ni) , UVcontour(1, Nip1) - UVcontour(1, Ni) };
      NUGA::normalize<2>(PiPip1);
      double Norm2[] = { -PiPip1[1], PiPip1[0], 0. };  // {-b, a}

      double spiky = ::fabs(NUGA::normals_angle(Norm1, Norm2) - NUGA::PI);

      err = (spiky < ZERO_M);
    }
  }

  return err;
}

void K_OCC::CADviaOCC::__split_surface_of_revolution(const OCCSurface* face, K_FLD::IntArray& connectB, K_FLD::FloatArray& pos3D, std::map<E_Int, std::pair<E_Int, E_Int> >& seam_nodes)
{
  
#ifdef DEBUG_CAD_READER
  std::cout << "splitting surface ..." << std::endl;
#endif
  
  //1. removing any duplicate
  std::set<K_MESH::NO_Edge> uedges;  
  seam_nodes.clear();

  double midway = face->_normalize_domain ? NUGA::PI : 0.5;
  
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
        __add_seam_node(face, pos3D, N0, seam_nodes);
      E_Int N1 = e.node(1);
      if (seam_nodes.find(N1) == seam_nodes.end())
        __add_seam_node(face, pos3D, N1, seam_nodes);
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
      E_Float u,v;
      face->parameters(pos3D.col(N1), u, v);
      
      E_Float * p = &u;
      if (face->_isVClosed && !face->_isUClosed) p = &v;
      
      E_Int e[] = {N1, E_IDX_NONE};
      if (*p < midway)
        e[1] = it0->second.second;
      else
        e[1] = it0->second.first;
      new_connect.pushBack(e, e+2);
    }
    else if (it0 == seam_nodes.end() && it1 != seam_nodes.end()) // seam-connected edge
    {
      E_Float u,v;
      //E_Int err = 
      face->parameters(pos3D.col(N0), u, v);
      
      E_Float *p = &u;
      if (face->_isVClosed && !face->_isUClosed) p = &v;
      
      E_Int e[] = {N0, E_IDX_NONE};
      if (*p < midway)
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

#ifdef DEBUG_CAD_READER
  medith::write("revol.mesh", pos3D, connectB);
#endif  
}

void K_OCC::CADviaOCC::__add_seam_node
(OCCSurface const *face, K_FLD::FloatArray& pos3D, E_Int N0,
 std::map<E_Int, std::pair<E_Int, E_Int> >& seam_nodes)
{ 
  //create the points
  E_Float u,v;
  //E_Int err = 
  face->parameters(pos3D.col(N0), u, v);

  E_Float Pt[3];
  E_Float eps = 1.e-12;

  if (face->_isUClosed && !face->_isVClosed)
  {
    E_Float u0 = face->_U0 + eps;
    E_Float u1 = face->_U1 -eps;

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
    E_Float v0 = face->_V0 + eps;
    E_Float v1 = face->_V1 - eps;
    
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
