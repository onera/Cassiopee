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
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef __NGON_DEBUG_H__
#define	__NGON_DEBUG_H__

#include "Fld/ngon_t.hxx"
#include "IO/io.h"
#include <sstream>
#include "Nuga/Delaunay/Triangulator.h"
#include "Connect/MeshTool.h"
#include "TRI_debug.h"

#define TEMPLATE_COORD_CONNECT template <typename Coordinate_t, typename Connectivity_t>
#define NGON_DBG_CLASS NGON_debug<Coordinate_t,Connectivity_t>

#define COL_DEFAULT 0
#define COL_RED     1
#define COL_GREEN   2
#define COL_NEXT    5
#define COL_NEXT2   7

TEMPLATE_COORD_CONNECT
class NGON_debug
{
//
public:
  typedef ngon_t<Connectivity_t>               ngon_type;    
  typedef K_FLD::ArrayAccessor<Coordinate_t>   ACoordinate_t;
  typedef K_FLD::ArrayAccessor<Connectivity_t> AConnectivity_t;
//    
public:
  static void write(const char* fname, const ACoordinate_t& coord, const ngon_type& ng, const std::vector<bool>*keep=0, bool clean=false);
  static void write_external_phs(const char* fname, const ACoordinate_t& coord, const ngon_type& ng);
  static void write(const char* fname, const ACoordinate_t& coord, const ngon_unit& pgs);
  
  static void enabling_write(const char* fname, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect, const char* elt_type);
  
  static void draw_PHT3(const K_FLD::FloatArray& coord, const ngon_type& ng, E_Int PHi);
  static void draw_PHT3(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::map<E_Int, Vector_t<E_Int> >& PHT3s, E_Int PHi, bool localid=false);
  static void draw_PHT3(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3, E_Int PHi);
  static void draw_PHT3s(const char* fname, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::map<E_Int, Vector_t<E_Int> >& PHT3s, Vector_t<E_Int>& colors, const Vector_t<bool>* flag=0);
  static void draw_PGT3s(const K_FLD::FloatArray& coord,const ngon_unit& PGs);
  static void draw_PGT3(const char* fname, const K_FLD::FloatArray& coord,const ngon_unit& PGs, E_Int PGi, K_FLD::FloatArray* PGnormals=0);
  static void draw_PH(const char* fname, const K_FLD::FloatArray& coord, const ngon_type& ng, E_Int i);
  static void draw_wired_PH(const char* fname, const K_FLD::FloatArray& crd, const ngon_type& ng, E_Int ith, E_Int index_start);
  static void draw_PHs(const char* fname, const K_FLD::FloatArray& coord, const ngon_type& ng, const Vector_t<E_Int>& ids);
  static void draw_PG (const K_FLD::FloatArray& crd, const ngon_unit& PGs, E_Int PGi, bool localid=false);
  static void draw_PGs(const char* fname, const K_FLD::FloatArray& crd, const ngon_unit& PGs, const Vector_t<E_Int>& PGis, bool localid = false);
  static void draw_PGs(const char* fname, const K_FLD::FloatArray& crd, const ngon_unit& PGs, bool localid = false);
  static void draw_PG_to_T3(E_Int PGi, const Vector_t<E_Int>& nT3_to_PG, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3);
  
  static void draw_wired_PG(const char* fname, const K_FLD::FloatArray& coord, const ngon_unit& PGs, E_Int ith, E_Float *normal);
  
  template <typename Triangulator_t>
  static void highlight_PH(const ngon_type& ng, const K_FLD::FloatArray& coord, E_Int PHi);
  
  static void get_PHT3_neighbors(E_Int PHT3i, const std::map<E_Int, Vector_t<E_Int> >& PHT3s,
                                 const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, K_FLD::IntArray& neighbors, bool both_orient=false);
  
  static void get_PHT3_points (E_Int PHT3i, const K_FLD::IntArray& connectT3, 
                                const std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s, 
                                const K_FLD::FloatArray& icrd, K_FLD::FloatArray& ocrd);
  
  static void __get_historical_PHs (const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::map<E_Int, Vector_t<E_Int> >& PHT3s, E_Int PHi,
                                    E_Int shift, E_Int nb_pgs1, const K_FLD::IntArray& F2E, K_FLD::IntArray& anc_PH, std::vector<E_Int>&nT3_to_oPG,
                                    std::set<E_Int>& PHs1, std::set<E_Int>& PHs2);
  
  static void extract_pgs_of_type(E_Int type, const char* fname, const ngon_type& ng, const K_FLD::FloatArray& crd);
  
};

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::write
(const char* fname, const ACoordinate_t& coord, const ngon_type& ng, const std::vector<bool>*keep, bool clean)
{
  if (ng.PGs.size() == 0)
    return;

  Connectivity_t cNGON;
  ngon_type ngtmp = ng;
  Coordinate_t crd(coord.array());
  if (clean)
  {
    ngon_type::clean_connectivity(ngtmp, crd);
    ngon_type::compact_to_used_nodes(ngtmp.PGs, crd);
  }

  if (keep == 0)
  {
    cNGON.pushBack(ngtmp.PGs._NGON);
    cNGON.pushBack(ngtmp.PHs._NGON);
  }
  else
  {
    ngon_type NGout;
    Vector_t<E_Int> npgids;
    ngtmp.select_phs(ngtmp, *keep, npgids, NGout);//fixme : ngtmp is refered as object AND argument
    cNGON.pushBack(NGout.PGs._NGON);
    cNGON.pushBack(NGout.PHs._NGON);
  }

  if (cNGON.getSize() == 0)
    return;

  //std::cout << cNGON << std::endl;

  std::ostringstream fullname;
  fullname << fname;

#ifdef WIN32
  fullname << ".tp";
#else
  fullname << ".plt";
#endif
  
  MIO::write(fullname.str().c_str(), crd, cNGON, "NGON");
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::write_external_phs
(const char* fname, const ACoordinate_t& coord, const ngon_type& ng)
{
  ngon_unit PHex;
  std::vector<E_Int> oIds;
  ng.PHs.extract_of_type(INITIAL_SKIN, PHex, oIds);
  
  if (oIds.empty()) return;
  
  ngon_type ngex;
  ngex.PHs=PHex;
  ngex.PGs=ng.PGs;
  
  Connectivity_t cNGON; 
  cNGON.pushBack(ngex.PGs._NGON);
  cNGON.pushBack(ngex.PHs._NGON);
  
  if (cNGON.getSize() == 0)
    return;
  
  //std::cout << cNGON << std::endl;
  
  MIO::write(fname, coord.array(), cNGON, "NGON");
}
  
///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::write
(const char* fname, const ACoordinate_t& coord, const ngon_unit& pgs)
{
  if (pgs._NGON.size() == 0)
  {
    std::cout << "input pgs is empty !!!" << std::endl;
    return;
  }
  
  pgs.updateFacets();
  
  Connectivity_t cNGON; 
  cNGON.pushBack(pgs._NGON);
    
  if (cNGON.getSize() == 0)
    return;
  size_t sz = pgs.size();
  // create fictive elements, 1 by pgs
  Vector_t<E_Int> PHs(2, 0);
  PHs[0]=pgs.size();
  for (size_t i = 0; i < sz; ++i)
  {
    PHs.push_back(1);
    PHs.push_back(i+1);
  }
  PHs[1]=PHs.size();
  
  cNGON.pushBack(PHs);
  
  //std::cout << cNGON << std::endl;
  
  std::ostringstream fullname;
  fullname << fname;

#ifdef WIN32
  fullname << ".tp";
#else
  fullname << ".plt";
#endif
  
  MIO::write(fullname.str().c_str(), coord.array(), cNGON, "NGON");
}

/*
///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::write
(const char* fname, const ACoordinate_t& coord, const ngon_type& ng, const Vector_t<bool>& face_flag)
{
  Connectivity_t cNGON;
  
  ngon_unit pgs;
  //fictive elements (1pg -> 1 element) to make the IO work
  Vector_t<E_Int> phs(2,0);
  
  size_t sz = ng.PGs.size();
  for (size_t i = 0; i < sz; ++i)
  {
    if (face_flag[i])
      pgs.add(ng.PGs, i);
  }
  
  sz = pgs.size();
  for (size_t i = 0; i < sz; ++i)
  {
    phs.push_back(1);
    phs.push_back(i+1);
  }
  
  phs[0]=sz;
  phs[1]=phs.size()-2;
      
  cNGON.pushBack(pgs._NGON);
  cNGON.pushBack(phs);
  
  size_t z = pgs.size();
  
  MIO::write(fname, coord.array(), cNGON, "NGON");
}*/

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::enabling_write
(const char* fname, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect, const char* elt_type)
{
  //if (_enabled)
    //MIO::write(fname, coord, connect, elt_type);
}

#define dABS(x) ((x<0)? -x : x)
///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PHT3
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::map<E_Int, Vector_t<E_Int> >& PHT3s, E_Int PHi, bool localid)
{
  K_FLD::IntArray cT3;
  std::map<E_Int, Vector_t<E_Int> >::const_iterator it=PHT3s.find(PHi);
  if (it == PHT3s.end())
    return;
  
  for (size_t i=0; i < it->second.size(); ++i)
  {
    E_Int Ti = dABS(it->second[i]);
    cT3.pushBack(connectT3.col(Ti), connectT3.col(Ti)+3);
  }
  
  std::ostringstream o;
  o << "PHT3_" << PHi << ".mesh";
    
  if (localid){
    
    K_FLD::FloatArray tmpCrd (coord);
    K_FLD::IntArray tmpCnt (cT3);
    
    std::vector<E_Int> nids;
    K_CONNECT::MeshTool::compact_to_mesh(tmpCrd, tmpCnt, nids);
    
    MIO::write(o.str().c_str(), tmpCrd, tmpCnt, "TRI");
  }
  else
    MIO::write(o.str().c_str(), coord, cT3, "TRI");
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PHT3
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3, E_Int PHi)
{
  std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >::const_iterator itPH = PH_to_PGT3.find(PHi);
  if (itPH == PH_to_PGT3.end())
    return;
  
  std::vector<E_Int> colors;
  K_FLD::IntArray cT3;
  std::map<E_Int, Vector_t<E_Int> >::const_iterator itPG;
  for (itPG= itPH->second.begin(); itPG != itPH->second.end(); ++itPG)
  {
    for (size_t j=0; j < itPG->second.size(); ++j)
      cT3.pushBack(connectT3.col(itPG->second[j]), connectT3.col(itPG->second[j])+3);
    colors.resize(cT3.cols(), itPG->first);
  }
  
  std::ostringstream o;
  o << "PHT3col_" << PHi << ".mesh";
  MIO::write(o.str().c_str(), coord, cT3, "TRI", 0, &colors);
  
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PHT3 (const K_FLD::FloatArray& coord, const ngon_type& ng, E_Int PHi)
{
  ng.PGs.updateFacets();
  ng.PHs.updateFacets();
  
//  if (PHi < ng.PHs.size())
//    return;
  
  K_FLD::IntArray connectT3;
  DELAUNAY::Triangulator t;
  std::vector<E_Int> colors;
  E_Int err = K_MESH::Polyhedron<UNKNOWN>::triangulate(t, ng.PGs, ng.PHs.get_facets_ptr(PHi), ng.PHs.stride(PHi), coord, connectT3, colors, true, false); // PH -> PHT3
  
  std::ostringstream o;
  o << "PHT3ized_" << PHi << ".mesh";
  MIO::write(o.str().c_str(), coord, connectT3, "TRI");
}


///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PHT3s
(const char* fname, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::map<E_Int, Vector_t<E_Int> >& PHT3s, Vector_t<E_Int>& colors, const Vector_t<bool>* flagPHT3)
{
  K_FLD::IntArray cT3;
  std::map<E_Int, Vector_t<E_Int> >::const_iterator it=PHT3s.begin();
  Vector_t<E_Int> cols;
  size_t count(0);
  for (; it != PHT3s.end(); ++it, ++count)
  {
    if (flagPHT3 && (*flagPHT3)[count]==false)
      continue;
    
    for (size_t i=0; i < it->second.size(); ++i)
    {
    E_Int Ti = dABS(it->second[i]);
    cT3.pushBack(connectT3.col(Ti), connectT3.col(Ti)+3);
    }
    cols.resize(cT3.cols(), colors[count]);
  }
  
  MIO::write(fname, coord, cT3, "TRI", 0, &cols);
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PGT3s(const K_FLD::FloatArray& coord,const ngon_unit& PGs)
{
  K_FLD::IntArray connectT3;
  Vector_t<E_Int> T3_to_nPG;
  E_Int err = ngon_type::template triangulate_pgs<DELAUNAY::Triangulator>(PGs, coord, connectT3, T3_to_nPG, true, false);
  if (err) std::cout << "failed on draw_PGT3s" << std::endl;
  if (err) return;
  Vector_t<E_Int> colors = T3_to_nPG;
  for (size_t i = 0; i< colors.size(); ++i)colors[i]+=1;//to use 0 for distinguish externals
  
  for (size_t i = 0; i < T3_to_nPG.size(); ++i)
  {
    E_Int pg = T3_to_nPG[i];
    if (!PGs._type.empty() && PGs._type[pg] == INITIAL_SKIN)
      colors[i]=0;
  }
  Vector_t<bool> keep(colors.size(), true);
  MIO::write("ngT3.mesh", coord, connectT3, "TRI", &keep, &colors);
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PGT3
(const char* fname, const K_FLD::FloatArray& coord,const ngon_unit& PGs, E_Int PGi, K_FLD::FloatArray* PGnormals)
{
  K_FLD::IntArray connectT3;
    
  PGs.updateFacets();
  DELAUNAY::Triangulator t;
  t.run(coord, PGs.get_facets_ptr(PGi), PGs.stride(PGi), 1, connectT3, true/*do ot shuffle*/, false);
     
  MIO::write(fname, coord, connectT3, "TRI");
  
  if (PGnormals)
  {
    K_FLD::FloatArray norms;
    for (size_t i=0; i< connectT3.cols(); ++i)
      norms.pushBack(PGnormals->col(PGi), PGnormals->col(PGi)+3);
    
    TRI_debug::write_wired("WngT3.mesh", coord, connectT3, norms, 0,0,true/*localid*/);
  }
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PG
(const K_FLD::FloatArray& crd, const ngon_unit& PGs, E_Int PGi, bool localid)
{
  PGs.updateFacets();

  ngon_unit pg;
  pg.add(PGs.stride(PGi), PGs.get_facets_ptr(PGi));

  ngon_type ng(pg);
  K_FLD::IntArray cnt;
  ng.export_to_array(cnt);

  const K_FLD::FloatArray *pCrd(&crd);
  K_FLD::FloatArray tmpCrd;
  //if (localid){
  //  tmpCrd = crd;
  //  ngon_type::compact_to_used_nodes(ng.PGs, tmpCrd);
  //  pCrd = &tmpCrd;
  //  //if (colors) K_CONNECT::IdTool::compact(*colors, nids);
  //}

  std::ostringstream o;
#ifndef WIN32
  o << "one_pg_" << PGi << ".plt";
#else
  o << "one_pg_" << PGi << ".tp";
#endif
  MIO::write(o.str().c_str(), *pCrd, cnt, "NGON");
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PGs(const char* fname, const K_FLD::FloatArray& crd, const ngon_unit& PGs, const Vector_t<E_Int>& PGis, bool localid)
{
  ngon_unit pg;
  for (size_t i = 0; i < PGis.size(); ++i)
    pg.add(PGs.stride(PGis[i]), PGs.get_facets_ptr(PGis[i]));
    
  ngon_type ng(pg);
  K_FLD::IntArray cnt;
  ng.export_to_array(cnt);
  std::ostringstream o;
#ifndef WIN32
  o << fname  << ".plt";
#else
  o << fname << ".tp";
#endif
  
  MIO::write(o.str().c_str(), crd, cnt, "NGON");
}
///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PGs(const char* fname, const K_FLD::FloatArray& crd, const ngon_unit& PGs, bool localid)
{    
  ngon_type ng(PGs);
  K_FLD::IntArray cnt;
  ng.export_to_array(cnt);
  std::ostringstream o;
  
#ifndef WIN32
  o << fname  << ".plt";
#else
  o << fname << ".tp";
#endif
  
  MIO::write(o.str().c_str(), crd, cnt, "NGON");
}
/*
///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PH
(const char* fname, const K_FLD::FloatArray& coord, const ngon_type& ng, E_Int i)
{
  ng.PGs.updateFacets();
  ng.PHs.updateFacets();
  
  E_Int nbf = ng.PHs.stride(i);
  const E_Int* pt = ng.PHs.get_facets_ptr(i);
  
  Vector_t<E_Int> indices, oIds;
  indices.insert(indices.end(), pt, pt+nbf);
  for (size_t i=0; i < indices.size(); ++i)indices[i]-=1;
  
  ngon_unit pgs, ph;
  ng.PGs.extract(indices, pgs, oIds);
  
  Vector_t<E_Int> molec_ph;
  molec_ph.push_back(nbf);
  for (size_t i=0; i < nbf; ++i)
    molec_ph.push_back(i+1);
  
  ph.add(molec_ph);
  
  Connectivity_t cNGON; 
  cNGON.pushBack(pgs._NGON);
  cNGON.pushBack(ph._NGON);
  
  if (cNGON.getSize() == 0)
    return;
  
  //std::cout << cNGON << std::endl;
  
  MIO::write(fname, coord, cNGON, "NGON");
  
}*/

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PH
(const char* fname, const K_FLD::FloatArray& coord, const ngon_type& ng, E_Int i)
{
  ng.PGs.updateFacets();
  ng.PHs.updateFacets();
  
  ngon_unit ph;
  ph.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
  
  ngon_type one_ph(ng.PGs, ph);
  Vector_t<E_Int> pgnids, phnids;
  one_ph.remove_unreferenced_pgs(pgnids, phnids);
  
  Coordinate_t crd(coord);
  ngon_type::compact_to_used_nodes(one_ph.PGs, crd);
  
  Connectivity_t cnt;
  one_ph.export_to_array(cnt);
  
  MIO::write(fname, crd, cnt, "NGON");
  
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_wired_PH
(const char* fname, const K_FLD::FloatArray& crd, const ngon_type& ng, E_Int ith, E_Int index_start)
{
  K_FLD::IntArray cnt;
  
  E_Int nb_nodes, E[2];
  E_Int nb_pgs(ng.PHs.stride(ith));
  const E_Int* ptrF = ng.PHs.get_facets_ptr(ith);
  
  for (size_t i=0; i < nb_pgs; ++i)
  {
    E_Int Fi = *(ptrF+i) -1;
    nb_nodes = ng.PGs.stride(Fi);
    
    for (size_t j=0; j < nb_nodes; ++j)
    {
      E[0]=ng.PGs.get_facet(i,j)-index_start;
      E[1]=ng.PGs.get_facet(i,(j+1)%nb_nodes)-index_start;
      cnt.pushBack(E, E+2);
    }
  }
  
  
  MIO::write(fname, crd, cnt, "BAR");
  
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_wired_PG(const char* fname, const K_FLD::FloatArray& coord, const ngon_unit& PGs, E_Int ith, E_Float *normal)
{
  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(coord);
  K_FLD::IntArray connectE;
  E_Int n0, n1;
  E_Float P0[3], P1[3], Lmin(K_CONST::E_MAX_FLOAT), L2;
  
  E_Int nb_nodes, E[2];
  const E_Int* pNi = PGs.get_facets_ptr(ith);
  nb_nodes = PGs.stride(ith);
  
  for (size_t j=0; j < nb_nodes; ++j)
  {
    E[0]=*(pNi+j)-1;
    E[1]=*(pNi+(j+1)%nb_nodes)-1;
    connectE.pushBack(E, E+2);
    L2 = K_FUNC::sqrDistance(coord.col(E[0]), coord.col(E[1]), 3);
    Lmin = (L2 < Lmin) ? L2 : Lmin;
  }
  
  Lmin = 0.5*::sqrt(Lmin);
  
  K_MESH::Polygon::iso_barycenter<acrd_t, 3 >(acrd, pNi, nb_nodes, 1, P0);
  
  K_FLD::FloatArray crd(coord);
  
  K_FUNC::sum<3>(1., P0, Lmin, normal, P1);
  crd.pushBack(P0, P0+3);
  E[0]=crd.cols()-1;
  crd.pushBack(P1, P1+3);
  E[1]=crd.cols()-1;
  connectE.pushBack(E, E+2);
  
  MIO::write(fname, crd, connectE, "BAR");
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PHs(const char* prefixname, const K_FLD::FloatArray& coord, const ngon_type& ng, const Vector_t<E_Int>& ids)
{
  ng.PGs.updateFacets();
  ng.PHs.updateFacets();
  
  ngon_unit ph;
  
  //
  for (size_t i=0; i< ids.size(); ++i)
    ph.add(ng.PHs.stride(ids[i]), ng.PHs.get_facets_ptr(ids[i]));
  
  //
  ngon_type phs(ng.PGs, ph);
  
  Vector_t<E_Int> pgnids, phnids;
  phs.remove_unreferenced_pgs(pgnids, phnids);
  
  Coordinate_t crd(coord);
  ngon_type::compact_to_used_nodes(phs.PGs, crd);
  
  Connectivity_t cnt;
  phs.export_to_array(cnt);
  
  MIO::write(prefixname, crd, cnt, "NGON"); 
}

///
TEMPLATE_COORD_CONNECT
template <typename Triangulator_t>
void NGON_DBG_CLASS::highlight_PH
(const ngon_type& ng, const K_FLD::FloatArray& coord, E_Int PHi)
{
  E_Int s = ng.PHs.stride(PHi);
  const E_Int* ptr = ng.PHs.get_facets_ptr(PHi);
  std::set<E_Int> pgis;
  for (size_t i = 0; i < s; ++i) pgis.insert(ptr[i]);
  
  K_FLD::IntArray connectT3;
  Vector_t<E_Int> T3_to_nPG;
  E_Int err = ngon_type::template triangulate_pgs<Triangulator_t>(ng.PGs, coord, connectT3, T3_to_nPG);
  
  Vector_t<E_Int> colors(connectT3.cols(), COL_DEFAULT);
    
  for (size_t i = 0; i < T3_to_nPG.size(); ++i)
  {
    E_Int pg = T3_to_nPG[i];
    if (pgis.find(pg+1) != pgis.end())
    {
      colors[i]=COL_GREEN;
      if (ng.PGs._external[pg])
        colors[i]=COL_RED;
    }
  }
  MIO::write("higlightPH.mesh", coord, connectT3, "TRI", 0, &colors);
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::get_PHT3_neighbors
(E_Int PHT3i, const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const K_FLD::FloatArray& coord,
 const K_FLD::IntArray& connectT3, K_FLD::IntArray& neighbors, bool both_orient)
{
  K_FLD::IntArray connect;
  std::map<E_Int, Vector_t<E_Int> >::const_iterator it = PHT3s.find(PHT3i);
  
  if (it == PHT3s.end())
    return;
  
  const std::vector<E_Int>& T3s = it->second;
  
  for (size_t i=0; i < T3s.size(); ++i)
  {
    E_Int Ti = T3s[i];
    std::ostringstream o;
    o << "neigh_" << Ti << ".mesh";
    TRI_debug::get_T3_neighbors(o.str().c_str(), Ti, coord, connectT3, neighbors, both_orient); 
  }
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::draw_PG_to_T3
(E_Int PGi, const Vector_t<E_Int>& nT3_to_PG, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3)
{
  K_FLD::IntArray connect;
  for (size_t i=0; i< nT3_to_PG.size(); ++i)
  {
    if (nT3_to_PG[i]==PGi)
      connect.pushBack(connectT3.col(i), connectT3.col(i)+3);
  }
  
  std::ostringstream o;
  o << "PG_" << PGi << "_toT3.mesh";
  MIO::write(o.str().c_str(), coord, connect, "TRI");
}

///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::__get_historical_PHs
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::map<E_Int, Vector_t<E_Int> >& PHT3s, E_Int PHi,
 E_Int shift, E_Int nb_pgs1, const K_FLD::IntArray& F2E, K_FLD::IntArray& anc_PH, std::vector<E_Int>&nT3_to_oPG,
 std::set<E_Int>& PHs1, std::set<E_Int>& PHs2)
{
  //
  typedef std::map<E_Int, Vector_t<E_Int> > map_t;
  typedef Vector_t<E_Int> vec_t;
  typedef std::set<E_Int> set_t;
  
  map_t::const_iterator it = PHT3s.find(PHi);  
  if (it == PHT3s.end()) return;
  
  const vec_t& T3s = it->second;
  
  set_t* PHs = 0;
  
  E_Int nb_t3s = T3s.size();
  
  for (size_t i = 0; (i< nb_t3s); ++i)
  {
    const E_Int& t = T3s[i];
    E_Int T = (t<shift) ? t : t - shift;
    //E_Int I = (t<shift) ? 1 : 0;
    
    if (T >= nT3_to_oPG.size()) continue;
    
    E_Int wPG = nT3_to_oPG[T];
    if (wPG >= F2E.cols()) continue;
    
    //E_Int hPH = F2E(I, wPG);
    
    PHs = (wPG < nb_pgs1) ? &PHs1 : &PHs2;
    
    if (F2E(0, wPG) != E_IDX_NONE) PHs->insert(F2E(0, wPG));
    if (F2E(1, wPG) != E_IDX_NONE) PHs->insert(F2E(1, wPG)); 
  }
}


///
TEMPLATE_COORD_CONNECT
void NGON_DBG_CLASS::extract_pgs_of_type(E_Int type, const char* fname, const ngon_type& ng, const K_FLD::FloatArray& crd)
{
  ngon_unit pg_ext;
  Vector_t<E_Int> oids;
  //
  ng.PGs.extract_of_type(type, pg_ext, oids);
  write(fname, ACoordinate_t(crd), pg_ext);
  draw_PGT3s(crd, pg_ext);
}


#endif
