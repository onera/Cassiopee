/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_XCELLN_HXX
#define NUGA_XCELLN_HXX

//#define DEBUG_XCELLN

#include "Nuga/include/collider.hxx"
#ifdef DEBUG_XCELLN
#include "Nuga/include/medit.hxx"
#endif
#include "Nuga/include/mesh_t.hxx"
#include "Nuga/include/polygon.hxx"

namespace NUGA
{

  enum eClassify {AMBIGUOUS=-1, IN = 0, X = 1, OUT = 2, UPPER_COL=3};

  template <typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  void __compact_to_box(zmesh_t const & z_mesh, 
  	                    std::vector< bound_mesh_t*> & mask_bits,
                        std::vector< std::vector<E_Int> > &mask_wall_ids)
  {
  	//using vmesh_t = NUGA::vmesh_t<bound_mesh_t::ARG1, bound_mesh_t::ARG2>; //i.e same type as boundaries of zmesh but as a pure view
  	// zone reduction
  	K_SEARCH::BBox3D z_box;
  	z_mesh.bbox(z_box);

  	int nmasks = mask_bits.size();
  	for (size_t m=0; m < nmasks; ++m)
  	{
      if (mask_bits[m] == nullptr) continue;
      
  	  int nbcells = mask_bits[m]->ncells();

  	  std::vector<bool> keep(nbcells, false);

  	  for (int i=0; i< nbcells; ++i)
  	  {
  		  K_SEARCH::BBox3D bxi;
  		  mask_bits[m]->bbox(i, bxi);

  		  if (K_SEARCH::BbTree<3>::boxesAreOverlapping(&z_box, &bxi, E_EPSILON))
  			  keep[i] = true;
  	  }
      
      mask_bits[m]->compress(keep);

//      std::vector<E_Int> nids;
//  	  //K_FLD::IntArray::compact(*mask_cnts[m], keep, nids);
//
//      // propagate new ids, NONE can be in so compress
//      std::vector<E_Int> tmp;
//  	  for (size_t u=0; u < mask_wall_ids[m].size(); ++u)
//      {
//        E_Int nwid = nids[mask_wall_ids[m][u]];
//        if (nwid != E_IDX_NONE)tmp.push_back(nwid);
//      }
//      mask_wall_ids[m] = tmp;
//      
//      //compact coordinates
//      std::vector<E_Int> pnids;
//      //K_CONNECT::MeshTool::compact_to_mesh(*mask_crds[m], *mask_cnts[m], pnids);
  	}
  }

  template <typename bound_mesh_t>
  void __build_mask_bits(const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
                         std::vector< std::vector<E_Int> > &mask_wall_ids,
                         const std::vector<E_Int>& z_priorities, E_Int rank_wnp,
                         std::vector< bound_mesh_t*> & mask_bits)
  {

  	mask_bits.clear();
  	//int nb_comps = mask_crds.size();
    
    //std::cout << "__build_mask_bits : 1 mask_wall_ids size : " <<  mask_wall_ids.size() << std::endl;

    // grabbing OP (WP are discarded) and WNP
  	for (size_t i=0; i <z_priorities.size(); ++i)
  	{
  	  int compi = z_priorities[i];
      
      //std::cout << "__build_mask_bits : comp " << compi << std::endl;
  	  
  	  bool is_prior = (i < rank_wnp);

      //std::cout << "is prior ? " << is_prior << " rank/rank_wnp " << i << "/" << rank_wnp << std::endl;

  	  if (!is_prior && mask_wall_ids[compi].empty()) continue; // no WNP
      
      //std::cout << "__build_mask_bits : 1 "  << std::endl;

  	  bound_mesh_t* bit = new bound_mesh_t(*mask_crds[compi], *mask_cnts[compi]);
      
  	  int nbcells = bit->ncells();
      
      //std::cout << "__build_mask_bits : 3 : nb of cells in mask : " << nbcells  << std::endl;

      // completely removed by __compact_to_box or only walls in it (WP are not required)
  	  bool discard_bit = ( (nbcells == 0) || ( is_prior && (nbcells == mask_wall_ids[compi].size()) ) );

  	  if (discard_bit) 
  	  {
        std::cout << "mask bit " << i << " is discarded" << std::endl;
  	  	delete bit; continue;
  	  }
      
      //
      
      //std::cout << "__build_mask_bits : 4 "  << std::endl;

#ifdef DEBUG_XCELLN
    {
      std::ostringstream o;
      o << "mask_0_" << i;
      medith::write<>(o.str().c_str(), bit->crd, bit->cnt);
    }
#endif
      
      //std::cout << "__build_mask_bits : 5 "  << std::endl;

      // reduce to OP (discard WP) when prior comp, keep only WNP otherwise
  	  std::vector<bool> keep(nbcells, is_prior);
  	  for (size_t u=0; u<mask_wall_ids[compi].size(); ++u )
        keep[mask_wall_ids[compi][u]]=!is_prior;
      
      //std::cout << "__build_mask_bits : 6 "  << std::endl;
  	  
  	  bit->compress(keep);
  	  if (bit->ncells() == 0) // completely gone
  	  {
  	  	delete bit; continue;
  	  }
      
      //std::cout << "__build_mask_bits : 7 "  << std::endl;
      
      if (!is_prior) // reverse WNPs
      {
        //std::cout << "reversing " << i << std::endl;
        bit->reverse_orient();
      }

  	  mask_bits.push_back(bit);
      
      //std::cout << "mask rank : " << mask_bits.size() -1 << std::endl;

#ifdef DEBUG_XCELLN
    {
      //std::cout << "ouput mask_1_ " << i << std::endl;
      std::ostringstream o;
      o << "mask_1_" << i;
      medith::write<>(o.str().c_str(), bit->crd, bit->cnt);
    }
#endif
  	}
    //std::cout << "__build_mask_bits : exit" << std::endl;
  }

  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  void __process_overlapping_boundaries(zmesh_t & z_mesh, std::vector< bound_mesh_t*> & mask_bits, E_Int rank_wnp, E_Float RTOL)
  {
  	// WARNING : assume consistent orientation

    // get zmesh boundary (not using the constructor to do so because we need the ancesort info)
    bound_mesh_t zbound;

    std::vector<E_Int> ancestor;
  	z_mesh.get_boundary(zbound, ancestor);

#ifdef DEBUG_XCELLN
    {
      std::ostringstream o;
      o << "m";
      medith::write(o.str().c_str(), z_mesh.crd, z_mesh.cnt);
    }
    {
      std::ostringstream o;
      o << "bound";
      medith::write(o.str().c_str(), zbound.crd, zbound.cnt);
    }
#endif

  	std::vector<COLLIDE::eOVLPTYPE> is_x1, is_x2;

    bool has_abut_ovlp = false;
    for (size_t m=0; m < mask_bits.size(); ++m)
    {
      is_x1.clear();
      is_x2.clear();

#ifdef DEBUG_XCELLN
    {
      std::ostringstream o;
      o << "mask_" << m;
      medith::write<>(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt);
    }
#endif
      
      bound_mesh_t& maski = *mask_bits[m];
      COLLIDE::compute_overlap(zbound.crd, zbound.cnt, 
      	                       maski.crd, maski.cnt, *(maski.get_localizer()), 
                               is_x1, is_x2, RTOL);

      // flag the attached PG of z_mesh as IN
      std::vector<E_Int> ids;
      for (size_t u=0; u< is_x1.size(); ++u) 
      {
      	bool appendit = (is_x1[u] == COLLIDE::ABUTTING && m >=  rank_wnp); // test on rank : ie. is a WNP
      	appendit     |= (is_x1[u] == COLLIDE::OVERSET  && m <   rank_wnp); // test on rank : ie. is not a WNP, it is as OVLP
      	if (appendit) ids.push_back(ancestor[u]);
      }
      if (!ids.empty()) 
      {
        z_mesh.set_type(IN, ids);
        has_abut_ovlp = true;
      }
      
      // discard overlap elts in masks
      int nbcells = mask_bits[m]->ncells();
      assert (is_x2.size() == nbcells);
      //std::cout << "nbcells/is_x2 size : " << nbcells << "/" << is_x2.size() << std::endl;
      std::vector<bool> keep(nbcells, true);
  	  for (size_t u=0; u<nbcells; ++u )
  	   	if (is_x2[u] == COLLIDE::ABUTTING || is_x2[u] == COLLIDE::OVERSET) keep[u] = false;
  	  
  	  mask_bits[m]->compress(keep);
  	  if (mask_bits[m]->ncells() == 0) // completely gone
  	  {
  	  	delete mask_bits[m]; mask_bits[m]=nullptr;
  	  }
    }

    // compress (so keep the sorting) the masks to the remaining ones
    std::vector< bound_mesh_t*> remaining_masks;
    for (size_t m=0; m < mask_bits.size(); ++m)
      if (mask_bits[m] != nullptr)
      	remaining_masks.push_back(mask_bits[m]);
    
   mask_bits = remaining_masks;

#ifdef DEBUG_XCELLN
    if (has_abut_ovlp)
    {
      std::cout << "has ABUTTING/OVERSET" << std::endl;
      std::ostringstream o;
      o << "m_in";
      medith::write(o.str().c_str(), z_mesh.crd, z_mesh.cnt, z_mesh.e_type.empty() ? nullptr : &z_mesh.e_type);
    }
#endif
  }

  ///
  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  void PREP_build_structures_and_reduce_to_zone(zmesh_t & z_mesh, 
  	                                            const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
                                                std::vector< std::vector<E_Int> > &mask_wall_ids,
                                                const std::vector<E_Int>& z_priorities, E_Int rank_wnp,
                                                std::vector< bound_mesh_t*> & mask_bits)
  {
    //std::cout << "PREP_build_structures_and_reduce_to_zone : 1" << std::endl;
  	
    // Fast return : if prioritary and there are no walls => return
  	// if (z_priorities.empty())
  	// {
  	//   bool has_walls = false;
  	//   for (size_t i=0; (i < mask_wall_ids.size()) && !has_walls; ++i) has_walls |= !mask_wall_ids[i].empty();
  	//   if (!has_walls) return;
   //  }
    
    //std::cout << "PREP_build_structures_and_reduce_to_zone : 2" << std::endl;
  
  	// build mask data structures (mesh object) : WP are discarded. Putting first decreasing OP, then remaining WNP
  	__build_mask_bits(mask_crds, mask_cnts, mask_wall_ids, z_priorities, rank_wnp, mask_bits);
    
#ifdef DEBUG_XCELLN
    for (size_t m=0; m < mask_bits.size(); ++m){
      std::ostringstream o;
      o << "mask1a_" << m;
      medith::write(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt);
    }
    std::cout << "PREP_build_structures_and_reduce_to_zone : 3" << std::endl;
#endif

    // reduce masks to pieces in zone box : TO REMOVE IF DONE BEFORE IN THE PYTHON
  	__compact_to_box(z_mesh, mask_bits, mask_wall_ids);
    
#ifdef DEBUG_XCELLN
    for (size_t m=0; m < mask_bits.size(); ++m){
      std::ostringstream o;
      o << "mask1b_" << m;
      medith::write(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt);
    }
    //std::cout << "PREP_build_structures_and_reduce_to_zone : 4" << std::endl;
#endif    

  	// 1. detecting overlaps, 
  	// 2. marking on these overlap regions zmesh' border elts as IN, 
  	// 3. discarding overlap boundaries of mask (not required anymore, neither to blank, nor to clip)
  	//__process_overlapping_boundaries(z_mesh, mask_bits, rank_wnp, RTOL);
    
#ifdef DEBUG_XCELLN
    std::cout << "PREP_build_structures_and_reduce_to_zone : 5" << std::endl;
#endif
  	// now we have reduced the meshes to the useful part, create & add localizers
  	for (size_t m=0; m < mask_bits.size(); ++m)
  	  mask_bits[m]->build_localizer();
  }

  ///
  inline double __project(double const * plane_pt, double const * plane_dir, double const * pt, double const * dir, double* proj_pt)
  {
    double PPt[3];
    K_FUNC::diff<3>(pt, plane_pt, PPt);
  
    double k =  - K_FUNC::dot<3>(PPt, plane_dir);
    double c = K_FUNC::dot<3>(plane_dir, dir);
  
    //assert(SIGN(c, EPSILON) != 0); //fixme
    k /= c;
  
    K_FUNC::sum<3>(1., pt, k, dir, proj_pt); //project

    return k;
  }


  ///
  template <typename aelt_t, typename bound_mesh_t>
  bool __are_colliding(const aelt_t& e1, const bound_mesh_t& mask_bit, const std::vector<E_Int>& cands, E_Int idx_start, E_Int & cid, double RTOL);
  
  ///
  template <> inline
  bool __are_colliding<NUGA::aPolygon, edge_mesh_t>
  (const NUGA::aPolygon& ae1, const edge_mesh_t& mask_bit, const std::vector<E_Int>& cands, E_Int idx_start, E_Int & cid, double RTOL)
  {
  	bool isx=false;
    cid = E_IDX_NONE;

  	// reduce mask to candidates
  	edge_mesh_t lmask = mask_bit;
    std::vector<bool> keep(lmask.ncells(), false);

  	for (size_t u = 0; u < cands.size(); ++u)
  	{
  		//std::cout << cands[u] - idx_start << std::endl;
  		keep[cands[u] - idx_start] = true;
  	}

  	lmask.compress(keep); // WARNING : must be 1-based here

  	double normal1[3];
  	ae1.normal<3>(normal1);
  	const double * plane_pt = ae1.m_crd.col(0); // choosing first node
    
#ifdef DEBUG_XCELLN
    //medith::write("cutter_front_before_projection", lmask.crd, lmask.cnt);
#endif

  	// project candidates on e1's plane => 2D problem
    STACK_ARRAY(double, lmask.crd.cols(), signed_dists);
  	for (int i=0; i < lmask.crd.cols(); ++i)
  	{
  	  // orthogonal projection on a plane (projection dir is the normal to the plane)
      // overwrite it
  	  signed_dists[i] = __project(plane_pt, normal1, lmask.crd.col(i), normal1, lmask.crd.col(i));
  	}

    // check if not too far
    bool is_close=false;
    double l21 = ae1.L2ref();
    double l22 = lmask.L2ref();
    double ATOL(RTOL * ::sqrt(std::min(l21, l22))), min_d(K_CONST::E_MAX_FLOAT);

    for (int i=0; (i < lmask.crd.cols()) && !is_close; ++i)
    {
      if (i < lmask.crd.cols()-1) is_close = signed_dists[i]*signed_dists[i+1] < 0.; // means crossing
      min_d = std::min(min_d, ::fabs(signed_dists[i]));
    }

    is_close |= (min_d < ATOL);
    if (!is_close) 
    {
      //std::cout << "lr1/lr2/ATOL : " << l21 << "/" << l22 << "/" << ATOL <<  " and min_d : " << min_d << std::endl;
      return false;
    }

    // close enough so go to 2D for real test

#ifdef DEBUG_XCELLN
    //medith::write("cutter_front_projected", lmask.crd, lmask.cnt);
#endif

  	// compute collision between e1 and each candidate until founding one collision

  	// a. triangulate e1, check if any lmask point falls into the triangulation
  	DELAUNAY::Triangulator dt;
  	ae1.triangulate(dt);

  	// b. go 2D both crd (a copy of it)  and lmask.crd

  	// Computes the fitting box coordinate system optimizing the view over the contour
    K_FLD::FloatArray P(3,3), iP(3,3);

    FittingBox::computeAFrame(normal1, P);
    iP = P;
    K_FLD::FloatArray::inverse3(iP);
    
    K_FLD::FloatArray crd2D(ae1.m_crd);
    FittingBox::transform(crd2D, iP);// Now we are in the fitting coordinate system.
    crd2D.resize(2, crd2D.cols()); // e1 is 2D now.
    FittingBox::transform(lmask.crd, iP);// Now we are in the fitting coordinate system.
    lmask.crd.resize(2, lmask.crd.cols()); // lmask is 2D now.

#ifdef DEBUG_XCELLN
    {
      K_FLD::FloatArray c1(crd2D), c2(lmask.crd);
      c1.resize(3, c1.cols(), 0.);
      c2.resize(3, c2.cols(), 0.);
      medith::write("subj2D", c1, ae1.begin(), ae1.nb_nodes(), 0);
      medith::write("cutter_front2D", c2, lmask.cnt);
    }
#endif
  	
  	// c. detect collisions

    int T1[3];
  	for (E_Int i=0; (i < ae1.nb_tris()) && !isx; ++i)
    {
      ae1.triangle(i, T1);
      const E_Float* P1 = crd2D.col(T1[0]);
      const E_Float* Q1 = crd2D.col(T1[1]);
      const E_Float* R1 = crd2D.col(T1[2]);

      for (int j=0; (j<lmask.ncells()) && !isx; ++j)
      {
      	K_MESH::Edge e = lmask.element(j);
      	const E_Float* P2 = lmask.crd.col(e.node(0));
        const E_Float* Q2 = lmask.crd.col(e.node(1));

        isx = K_MESH::Triangle::fast_intersectE2_2D(P1, Q1, R1, P2, Q2, ATOL);

        if (!isx) //deeper check
        {
          E_Float u00, u01;
          E_Int tx;
          E_Bool overlap;
          isx = K_MESH::Triangle::intersect<2>(P1,Q1,R1, P2, Q2, ATOL, true, u00, u01, tx, overlap);
        }
        if (isx) cid = j;
      }
    }

  	return isx;
  }

  ///
  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  bool __flag_collided_zone_cells(zmesh_t const & z_mesh, bound_mesh_t const & mask_bit, std::vector<E_Float>& z_xcelln, double RTOL)
  {
  	bool has_X = false;
    
#ifdef DEBUG_XCELLN
      medith::write<>("currentz", z_mesh.crd, z_mesh.cnt);
      medith::write<>("currentm", mask_bit.crd, mask_bit.cnt);
#endif
  
    E_Int nbcells = z_mesh.ncells();

    using loc_t = typename zmesh_t::loc_t;
    const loc_t& mask_loc = *(mask_bit.get_localizer());

    mask_bit.get_nodal_tolerance();
    z_mesh.get_nodal_tolerance();
  
    std::vector<E_Int> cands;
    
    for (E_Int i = 0; i < nbcells; ++i)
    {
      //std::cout << i << " over " << nbcells << std::endl;
      if (z_xcelln[i] != OUT) continue;

      // autonomous element
      auto ae1 = z_mesh.aelement(i);

      cands.clear();
      mask_loc.get_candidates(ae1, ae1.m_crd, cands, 1); //return as 0-based (fixme for volumic, was 1-based)
      if (cands.empty()) continue;
   
#ifdef DEBUG_XCELLN
     //NGDBG::draw_PGs("shape", crdp, *skinp, candsp);
      medith::write<ngon_type>("subj", z_mesh.crd, z_mesh.cnt, i);
      medith::write("cutter_front", mask_bit.crd, mask_bit.cnt, cands, 1);
#endif

      E_Int cid;
      bool is_x = __are_colliding(ae1, mask_bit, cands, 1, cid, RTOL);
      
      if (is_x)
      {
        z_xcelln[i] = -X; // minus to mark as new X
        z_mesh.set_flag(i, cands[cid]);
      }
      
      has_X |= is_x;

    }
    
#ifdef DEBUG_XCELLN
    std::vector<E_Int> xs;
    
    for (size_t i=0; i < z_xcelln.size(); ++i)
    {
      if (z_xcelln[i] == -X) xs.push_back(i);
    }
    if (xs.empty())
      std::cout << "NO COLLISIONS WITH CURRENT MASK" << std::endl;
    else
    {
      medith::write<ngon_type>("colliding_set", z_mesh.crd, z_mesh.cnt, xs);
      medith::write("flag_collided_zone_cells", z_mesh.crd, z_mesh.cnt, &z_xcelln);
    }
#endif

    return has_X;
  }

//  ///
//  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
//  void __flag_hidden_subzones(zmesh_t const & z_mesh, bound_mesh_t const & mask_bit, std::vector<E_Float>& z_xcelln)
//  {
//    // 1. compute neighborhood (lazy mode)
//    auto neighorz = z_mesh.get_neighbors();
//    
//  	// 2. incremental coloring, starting from UPPER_COL
//    
//    // 2.1 : current field with only IN cells and new Xs (-X)
//    std::vector<E_Float> cur_xcelln(z_xcelln.size(), OUT);
//    for (size_t i=0; i < z_xcelln.size(); ++i)
//    {
//      if (z_xcelln[i] == IN)
//        cur_xcelln[i] = IN;
//      else if (z_xcelln[i] == -X) // new X
//      {
//        cur_xcelln[i] = z_xcelln[i] = X; //reset on z_xcelln
//      }
//    }
//    // 2.2 : incremental coloring => new INs, some X
//    assert (UPPER_COL > OUT);
//    K_CONNECT::EltAlgo<typename zmesh_t::elt_t>::coloring(*neighorz, cur_xcelln, (E_Float)OUT, (E_Float)UPPER_COL);
//    // 2.3 : update z_xcelln with IN & X: RULE : IN > X > OUT
//    for (size_t i=0; i < z_xcelln.size(); ++i)z_xcelln[i] = std::min(z_xcelln[i], cur_xcelln[i]);
//    // 2.4 : check nb of subzones. STOP if < 2
//    E_Int nb_subzones = *std::max_element(ALL(cur_xcelln)) - UPPER_COL + 1;
//    if (nb_subzones <= 1) // no new subzones, just return
//      return;
//    
//    // 2.5 : classify subzones
//    // 2.5.1 extract "X river"
//    std::vector<bool> keep(z_mesh.ncells(), false);
//    for (size_t i=0; i < cur_xcelln.size(); ++i) if (cur_xcelln[i] == X)keep[i] = true;
//    zmesh_t xriver(z_mesh, true/*with neighbors*/);
//    xriver.compress(keep);
//
//#ifdef DEBUG_XCELLN
//    medith::write("xriver", xriver.crd, xriver.cnt);
//#endif
//
//    // 2.5.2 separate river boundaries by color AND CONNEXITY
//    // replace neigbors ids by their color
//    std::vector<E_Int> nids(cur_xcelln.size(), E_IDX_NONE);
//    for (size_t i=0; i < nids.size(); ++i)nids[i] = (E_Int)cur_xcelln[i];
//    auto neigh_river = xriver.get_neighbors();
//    neigh_river->change_indices(nids, 0 /*neigh_river is 0-based*/);
//    // get the boundaries using this coloring
//    std::map<E_Int, bound_mesh_t> color_to_riversides;
//    xriver.get_boundaries_by_color(color_to_riversides);
//   
//    color_to_riversides.erase(X);          //discard inner (inside the river)
//    color_to_riversides.erase(E_IDX_NONE); //discard domain's borders
//    
//#ifdef DEBUG_XCELLN
//    for (auto it : color_to_riversides)
//    {
//      std::ostringstream o;
//      o << "riverside_" << it.first;
//      medith::write(o.str().c_str(), it.second.crd, it.second.cnt);
//    }
//#endif
//    // 2.5.3 compute rough directions for the river itself and its boundaries
//    E_Float mdir[3], norm[3], origin[] = {0.,0.,0.};
//    
//    K_MESH::Polygon::compute_dir(mask_bit, origin, mdir);
//    // 2.5.4 classify to detect IN zones
//    STACK_ARRAY(eELTYPE, nb_subzones, z_color);
//    for (size_t c=0; c < nb_subzones; ++c)z_color[c]=OUT;
//    for (auto it : color_to_riversides)
//    {
//      E_Int color = it.first;
//      if (color == E_IDX_NONE) continue;
//      if (color == X) continue;
//      E_Float rs_dir[3];
//      K_MESH::Polygon::compute_dir(it.second, origin, rs_dir);
//      
//      E_Float ps = K_FUNC::dot<3>(mdir, rs_dir);
//      if (ps >= 0.) z_color[color] = IN;
//    }
//    // now update z_xcelln
//    for (size_t i=0; i < z_xcelln.size(); ++i)
//    {
//      if (z_color[cur_xcelln[i]] == IN) z_xcelln[i] = IN;
//    }
//
//#ifdef DEBUG_XCELLN
//    medith::write("coloring", z_mesh.crd, z_mesh.cnt, &z_xcelln);
//#endif
//  }
  
  template <typename elt_t>
  eClassify __classify(const elt_t& e1, const elt_t& e2);
  
  template <> inline
  eClassify __classify<K_MESH::aEdge>(const K_MESH::aEdge& e1, const K_MESH::aEdge& e2)
  {
    E_Float threshold = 0.75;//decide only with roughly colinear
    
    E_Float V1[3], V2[3];
    K_FUNC::diff<3>(e1.v2, e1.v1, V1);
    K_FUNC::diff<3>(e2.v2, e2.v1, V2);

    K_FUNC::normalize<3>(V1);
    K_FUNC::normalize<3>(V2);
    
    E_Float ps = K_FUNC::dot<3>(V1, V2);
    if (::fabs(ps) < threshold) return AMBIGUOUS;

    return (ps < 0.) ? IN : OUT;
  }
  
  ///
  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  void __flag_hidden_subzones(zmesh_t const & z_mesh, bound_mesh_t const & mask_bit, std::vector<E_Float>& z_xcelln)
  {
    // 1. compute neighborhood (lazy mode)
    auto neighborz = z_mesh.get_neighbors();
    
  	// 2. incremental coloring, starting from UPPER_COL
    
    // 2.1 : current field with only IN cells and new Xs (-X)
    std::vector<E_Float> cur_xcelln(z_xcelln.size(), OUT);

    for (size_t i=0; i < z_xcelln.size(); ++i)
    {
      if (z_xcelln[i] == IN)
        cur_xcelln[i] = IN;
      else if (z_xcelln[i] == -X) // new X
      {
        cur_xcelln[i] = z_xcelln[i] = X; //reset on z_xcelln
      }
    }
    
#ifdef DEBUG_XCELLN
    medith::write("initial_field_coloring", z_mesh.crd, z_mesh.cnt, &cur_xcelln);
#endif

    // 2.2 : incremental coloring => new INs, some X
    assert (UPPER_COL > OUT);
    K_CONNECT::EltAlgo<typename zmesh_t::elt_t>::coloring(*neighborz, cur_xcelln, (E_Float)OUT, (E_Float)UPPER_COL);
    
#ifdef DEBUG_XCELLN
    medith::write("colored_field_coloring", z_mesh.crd, z_mesh.cnt, &cur_xcelln);
#endif

    // 2.3 : update z_xcelln with IN & X: RULE : IN > X > OUT
    for (size_t i=0; i < z_xcelln.size(); ++i)z_xcelln[i] = std::min(z_xcelln[i], cur_xcelln[i]);
    // 2.4 : check nb of subzones. STOP if < 2
    E_Int nb_subzones = *std::max_element(ALL(cur_xcelln)) - UPPER_COL + 1;

#ifdef DEBUG_XCELLN
    // std::set<int> all_cols(ALL(cur_xcelln));
    // for (auto c : all_cols)
    //   std::cout << "colors : " << c << std::endl;
    // std::cout << "nb of supposed sub zones: " <<nb_subzones << std::endl;
#endif

    if (nb_subzones <= 1) // no new subzones, just return
      return;
    
    // 2.5 : classify subzones
    STACK_ARRAY(eClassify, nb_subzones, z_color);
    for (size_t i=0; i < nb_subzones; ++i)z_color[i] = AMBIGUOUS;
    
#ifdef DEBUG_XCELLN
    std::cout << "NB SUZONES : " << nb_subzones << std::endl;
#endif
    
    int ncell = z_mesh.ncells();
    int missing_col=nb_subzones;
    
    
    for (size_t i=0; (i < ncell) && missing_col; ++i)
    {
      if (cur_xcelln[i] != X) continue;
      
      int nneighs = neighborz->stride(i);
      const int* pneighs = neighborz->begin(i);

      for (size_t j=0; (j<nneighs); ++j)
      {
        if (pneighs[j] == E_IDX_NONE) continue;
        
        E_Int subid = cur_xcelln[pneighs[j]];
        
        if (subid == E_IDX_NONE)     continue;
        if (subid == X)              continue;
        
        subid -= UPPER_COL; // to be 0-based
        assert(subid > -1 && subid < nb_subzones);
        
        //std::cout << "subid " <<  subid << ": z_color[subid] : " << z_color[subid] << std::endl;
        
        if (z_color[subid] != AMBIGUOUS) continue;
        
        // subid is a unvisited subzone
        
        // get the boundary (edge/face)
        typename bound_mesh_t::elt_t bj;
        z_mesh.get_boundary(i, j, bj);
        E_Float L2rf = bj.L2ref(z_mesh.get_nodal_tolerance());
        //make it autonomous
        typename bound_mesh_t::aelt_t abj(bj, z_mesh.crd, L2rf);
        
        E_Int mid = z_mesh.flag[i]; //colliding mask element (1-based)
        if (mid == E_IDX_NONE) continue;
        --mid; // make it 0-based

        typename bound_mesh_t::aelt_t mask_e = mask_bit.aelement(mid);
        
#ifdef DEBUG_XCELLN
        K_FLD::FloatArray cc;
        cc.pushBack(abj.v1, abj.v1+3);
        cc.pushBack(abj.v2, abj.v2+3);
        K_FLD::IntArray ct(2,1,0);ct(1,0)=1;
        std::ostringstream o;
        o << "subid_" << subid;
        medith::write(o.str().c_str(), cc, ct, "BAR");
        std::vector<E_Int> toto;
        toto.push_back(mid);
        o.str("");
        o << "cutter_front_" << subid;
        medith::write(o.str().c_str(), mask_bit.crd, mask_bit.cnt, toto, 0);
        o.str("");
        o << "maskbit_" << subid;
        medith::write(o.str().c_str(), mask_bit.crd, mask_bit.cnt);
        
        medith::write<ngon_type>("ei", z_mesh.crd, z_mesh.cnt, i); 
#endif
        
        z_color[subid] = __classify(abj, mask_e);
        if (z_color[subid] != AMBIGUOUS)
        {
          --missing_col; //otherwise cary_on
        
#ifdef DEBUG_XCELLN
          std::cout << "subid : " << subid << " is in ? " << (z_color[subid] == IN) << std::endl;
#endif
        }
      }
    }
    
    // now update z_xcelln
    for (size_t i=0; i < z_xcelln.size(); ++i)
    {
      E_Int subid = cur_xcelln[i];
      
      if (subid == E_IDX_NONE)     continue;
      if (subid == X)              continue;
      
      subid -= UPPER_COL; // to be 0-based
      
      if (z_color[subid] == IN) z_xcelln[i] = IN;
    }

#ifdef DEBUG_XCELLN
    medith::write("flag_hidden_subzones", z_mesh.crd, z_mesh.cnt, &z_xcelln);
#endif
  }
  
  
  
  

  ///
  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  void MOVLP_zone_XcellN(zmesh_t const & z_mesh,
  	                     std::vector< bound_mesh_t*> const & mask_bits,
  	                     std::vector<E_Float>& z_xcelln, E_Float col_X, E_Float RTOL)
  {
  	// initialization of the ouput field
  	z_xcelln.clear();
  	z_xcelln.resize(z_mesh.ncells(), OUT);

  	// loop on mask bits
  	size_t nbits = mask_bits.size();
  	for (size_t i=0; i < nbits; ++i)
  	{
  	  // append z_xcelln with X
  	  bool has_X = __flag_collided_zone_cells(z_mesh, *(mask_bits[i]), z_xcelln, RTOL);
  	  if (!has_X) continue;

  	  // if there are at least 2 zones => mark as IN the hidden ones.
  	  __flag_hidden_subzones(z_mesh, *(mask_bits[i]), z_xcelln);
  	}

    // replace with output colors
    for (size_t i=0; i < z_xcelln.size(); ++i)
    {
      z_xcelln[i] = (z_xcelln[i] == X) ? col_X : (z_xcelln[i] == IN) ? 0. : 1./*(z_xcelln[i] == OUT) ? 1. : z_xcelln[i]*/;
    }

  }
}

#endif
