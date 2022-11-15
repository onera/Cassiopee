/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_SUPERMESH_HXX
#define NUGA_SUPERMESH_HXX

//#define XMATCH_DBG

#include "Nuga/include/clipper.hxx"


#ifdef XMATCH_DBG
#include "Nuga/include/medit.hxx"
#endif


namespace NUGA
{

struct field
{
  double* f;
  double* gradf[3];
  field() :f(nullptr) { gradf[0] = gradf[1] = gradf[2] = nullptr; }
};

template <typename zmesh_t> inline
void xmatch(const zmesh_t& m0, const zmesh_t& m1, double ARTOL, std::vector<E_Int>& anc0, std::vector<E_Int>& anc1, zmesh_t& xm)
{
  using aelt_t = typename zmesh_t::aelt_t;

  xm.clear();
  anc0.clear();
  anc1.clear();

  auto loc1 = m1.get_localizer();

  m0.get_nodal_metric2();
  m1.get_nodal_metric2();

  ngon_unit glob_edge_ids0;
  m0.build_global_edge_ids(glob_edge_ids0);
  ngon_unit glob_edge_ids1;
  m1.build_global_edge_ids(glob_edge_ids1);

 
  std::vector<E_Int> cands;
  std::vector<aelt_t> bits;
  size_t n, k;
  E_Int i,i2;
  DELAUNAY::Triangulator dt;

  // got 2D (in sub ref frame)
  K_FLD::FloatArray P(3, 3), iP(3, 3);
  bool ref2Dcomputed = false;
  K_FLD::FloatArray crd2D1, crd2D2;

  std::map<K_MESH::Edge, E_Int> key_to_id; // to ease gluing (see 'add' calls below)
  xm.crd = m0.crd;               // to ease gluing (see 'add' calls below)
  xm.crd.pushBack(m1.crd);       // to ease gluing (see 'add' calls below)

  E_Int nbcells0 = m0.ncells();
  E_Int nbpts0 = m0.crd.cols();
  E_Int nbptsI = xm.crd.cols();

  std::set<std::pair<E_Int, E_Int>> edge_to_node;
  std::vector<std::pair<double, E_Int>> lambda_to_node;
  std::vector<E_Int> xbit_ids;
  std::vector<E_Int> orient;
  std::set<K_MESH::Edge> w_oe_set;
  std::map<E_Int, E_Int> w_n_map;
  std::vector<long> poids;

  for (i = 0; i < nbcells0; ++i)
  {

    //if (i != 16106) continue;
    //if (i !=  142) continue;
    //if (i > 100) continue;
    auto ae0 = m0.aelement(i);
    int nnodes = ae0.nb_nodes();

    double surf0 = ae0.metrics();
    const double *norm0 = ae0.get_normal();

    cands.clear();
    loc1->get_candidates(ae0, ae0.m_crd, cands, 1, 1.e-2/*RTOL*/); //return as 1-based
    if (cands.empty()) continue;

#ifdef SUPERMESH_DBG
    medith::write<DELAUNAY::Triangulator>("cands", m1.crd, m1.cnt, &cands, 1);
    medith::write("ae0", ae0.m_crd, &ae0.m_nodes[0], ae0.m_nodes.size(), ae0.shift());
#endif

    ref2Dcomputed = false;
    edge_to_node.clear();
    xbit_ids.clear();

    for (n = 0; n < cands.size(); ++n)
    {
      //
      i2 = cands[n] - 1;
      auto ae1 = m1.aelement(i2);

      K_CONNECT::IdTool::shift(ae1.m_poids, nbpts0); //since xrd contains both crd0 & crd1

      double surfc = ae1.metrics();
      const double *normc = ae1.get_normal();

      double ps = NUGA::dot<3>(norm0, normc);
      if (::fabs(ps) < 0.9) continue; // must be nearly colinear

      //double toto = (180. / NUGA::PI) * NUGA::normals_angle(norm0, normc);
      //if (toto > 5. && toto != 180.) continue;

      bool ai1_reversed = false;
      if (ps < 0.) // revert one to have both with same orient (for NUGA::INTERSECT::INTERSECTION logic)
      {
        ae1.reverse_orient();
        ps = -ps;
        ai1_reversed = true;
      }

      bool check_for_equality = (::fabs(surf0 - surfc) < EPSILON) && (ps > 1. - EPSILON);

#ifdef SUPERMESH_DBG
      std::ostringstream o; o << "ae1_" << n;
      medith::write(o.str().c_str(), ae1.m_crd, &ae1.m_nodes[0], ae1.m_nodes.size(), ae1.shift());
#endif

      if (check_for_equality)
      {
        if (ae0 == ae1)
        {
          xm.add(ae0, false/*do capitalize crds*/);
          anc0.push_back(i);
          anc1.push_back(i2);
          break; // pure match found on a conformal mesh => no more candidate to check
        }
      }

      bits.clear();
      bool true_clip(false);
      aPolygon subj(ae0); // restart from a clean face as isolated_clip might erase subj and we are not in an incremental-boolean loop
      int err = NUGA::CLIP::isolated_clip<aelt_t, aelt_t>(subj, ae1, NUGA::INTERSECT::INTERSECTION, ARTOL, bits, true_clip);

      if (err) std::cout << "clipping erreor : " << i << "-th cell with " << n << "-th candidates" << std::endl;
      assert(!err);

      // check that bits are nearly colinear with subj : fixme : not found yet a better way to detect parasite bits

      for (size_t b = 0; b < bits.size(); ++b)
      {
        const double *normb = bits[b].get_normal();
        double ps = NUGA::dot<3>(norm0, normb);
        if (::fabs(ps) < 0.9) // must be nearly colinear
        {
          bits.clear(); 
          true_clip = false;
          break;
        }
      }

      if (true_clip)
      {
        E_Int le0, le1;
        auto ge0s = glob_edge_ids0.get_facets_ptr(i);
        auto ge1s = glob_edge_ids1.get_facets_ptr(i2);

        for (k = 0; k < bits.size(); ++k)
        {
          
          for (size_t p = 0; p < bits[k].m_poids.size(); ++p)
          {
            if (bits[k].m_poids[p] < 0) //intersection point : convert local to global szudzic
            {
              E_Int szudzik_val = -bits[k].m_poids[p] - 1;
              NUGA::szudzik_unpairing(szudzik_val, le0, le1);

              if (ai1_reversed)
              {
                E_Int n0 = ae1.m_crd.cols();
                if (le1 <= n0 - 2)
                  le1 = (n0 - 2) - le1;
              }
 
              K_MESH::Edge key(ge0s[le0], ge1s[le1]);
              auto it = key_to_id.find(key);

              if (it == key_to_id.end())
              {
                key_to_id[key] = xm.crd.cols();
                bits[k].m_poids[p] = xm.crd.cols();
                xm.crd.pushBack(bits[k].m_crd.col(p), bits[k].m_crd.col(p) + 3);
              }
              else
                bits[k].m_poids[p] = it->second;

              edge_to_node.insert(std::make_pair(le0, bits[k].m_poids[p]));
            }
          }
          
          xbit_ids.push_back(anc0.size());
          xm.add(bits[k], false/*do capitalize crds*/);
          anc0.push_back(i);
          anc1.push_back(i2);
          

#ifdef SUPERMESH_DBG
          std::ostringstream o;
          o << "bit_" << n << "_" << k;
          medith::write(o.str().c_str(), bits[k].m_crd, &bits[k].m_nodes[0], bits[k].m_nodes.size(), bits[k].shift());
#endif
        }

        continue;
      }

      // check for fully-in case : is ae0 fully inside ae1 ?

      if (!ref2Dcomputed)
      {
        ref2Dcomputed = true;
        NUGA::computeAFrame(norm0, P);
        iP = P;
        K_FLD::FloatArray::inverse3(iP);
        crd2D1 = ae0.m_crd;
        NUGA::transform(crd2D1, iP);
      }

      double Lref2 = ae0.Lref2();
      double ABSTOL = ::sqrt(Lref2) * 1.e-2;

      crd2D2 = ae1.m_crd;
      NUGA::transform(crd2D2, iP);

      aPolygon a0(std::move(crd2D1));
      aPolygon a1(std::move(crd2D2));

      if (ai1_reversed)
        a1.reverse_orient();

      NUGA::eClassify c = NUGA::CLASSIFY::classify2D(a0, a1, ABSTOL);
      assert(c != AMBIGUOUS);

      //
      if (c == IN)
      {
        xm.add(ae0, false/*do capitalize crds*/);
        xbit_ids.push_back(anc0.size());
        anc0.push_back(i);
        anc1.push_back(i2);

#ifdef SUPERMESH_DBG
        std::ostringstream o;
        o << "bit__ae0_in_" << n;
        medith::write(o.str().c_str(), ae0.m_crd, &ae0.m_nodes[0], ae0.m_nodes.size(), ae0.shift());
#endif

        break; // ae0 cannot be in several PG
      }
      else if (c == OUT)
      {
        // empty
      }
      else if (c == IN_1)
      {
        // project ae1 on ae0 first
        double zmean = 0;
        for (size_t u = 0; u < a0.m_crd.cols(); ++u) zmean += a0.m_crd(2, u);
        zmean /= a0.m_crd.cols();
        for (size_t u = 0; u < a1.m_crd.cols(); ++u)a1.m_crd(2, u) = zmean;
        NUGA::transform(a1.m_crd, P); // back to original ref frame  
        a1.m_poids = ae1.m_poids;

        xm.add(a1, false/*do capitalize crds*/);
        xbit_ids.push_back(anc0.size());
        anc0.push_back(i);
        anc1.push_back(i2);

#ifdef SUPERMESH_DBG
        std::ostringstream o;
        o << "bit__" << n << "_in_ae0";
        medith::write(o.str().c_str(), ae1.m_crd, &ae1.m_nodes[0], ae1.m_nodes.size(), ae1.shift());
#endif
      }

      crd2D1 = std::move(a0.m_crd);//give it back
      crd2D2 = std::move(a1.m_crd);//give it back
    }

    
    // COMPUTE THE RESIDUAL BIT : modified a0 minus {bits}
    
    // impacted ae0 : update ae0 with refining points
    if (!edge_to_node.empty())
    {
      crd2D1.clear();
      lambda_to_node.clear();

      for (size_t n = 0; n < nnodes; ++n)
        lambda_to_node.push_back(std::make_pair(double(n), ae0.m_poids[n]));

      double P0PX[3], P0P1[3];
      for (auto& e2n : edge_to_node)
      {
        int n0 = e2n.first;

        const double* P0 = xm.crd.col(ae0.m_poids[n0]);
        const double* PX = xm.crd.col(e2n.second);
        const double* P1 = xm.crd.col(ae0.m_poids[(n0+1)%nnodes]);

        NUGA::diff<3>(P1, P0, P0P1);
        NUGA::diff<3>(PX, P0, P0PX);

        double Lx = ::sqrt(NUGA::sqrNorm<3>(P0PX));
        double L = ::sqrt(NUGA::sqrNorm<3>(P0P1));

        assert(Lx <= L);

        double lambda = (Lx / L) + double(n0);
        
        lambda_to_node.push_back(std::make_pair(lambda, e2n.second));
      }

      std::sort(ALL(lambda_to_node));

      // from lambda_to_node to ae0

      poids.clear();
      
      double lambda_prev = -1.;
      for (auto& l2n : lambda_to_node)
      {
        double& lambda = l2n.first;
        if (lambda - lambda_prev > ZERO_M)
        {
          crd2D1.pushBack(xm.crd.col(l2n.second), xm.crd.col(l2n.second) + 3);
          poids.push_back(l2n.second);
        }
        lambda_prev = lambda;
      }

      ae0 = aPolygon(std::move(crd2D1));
      ae0.m_poids = poids;
      

#ifdef SUPERMESH_DBG
      std::ostringstream o;
      o << "modified_ae0";
      medith::write(o.str().c_str(), ae0.m_crd, &ae0.m_nodes[0], ae0.m_nodes.size(), ae0.shift());
      medith::write<DELAUNAY::Triangulator>("bits", xm.crd, xm.cnt, &xbit_ids, 0);
#endif
    }

    // boolean diff ae0 \ bound{bits} : add missing bits if any
    std::vector<std::deque<E_Int>> PGbs;
    orient.resize(xbit_ids.size(), 1);

    int err = K_MESH::Polygon::get_boundary(xm.crd, xm.cnt, xbit_ids/*0 based*/, PGbs, orient, w_oe_set, w_n_map);

    if (err) continue;

    assert(err == 0);
    assert(!PGbs.empty());

    std::vector<aelt_t> bits, tmpbits;

    bits.push_back(std::move(ae0));
    
    // incremental-boolean loop
    for (size_t p = 0; (p < PGbs.size()) && !bits.empty(); ++p)
    {

      poids.clear();
      crd2D2.clear();

      auto & PGb = PGbs[p];

      for (auto& N : PGb)
      {
        crd2D2.pushBack(xm.crd.col(N - 1), xm.crd.col(N - 1) + 3);
        poids.push_back(N - 1);
      }

      NUGA::aPolygon ae1f = std::move(crd2D2);
      ae1f.m_poids = poids;
      int nbits = bits.size();
      for (int b = 0; b < nbits; ++b)
      {
        if (bits[b].empty()) continue;

        if (bits[b] == ae1f)
        {
          // exact match for a diff => empty answer
          bits.clear();
          break;
        }

#ifdef SUPERMESH_DBG
        medith::write("ae1f", ae1f.m_crd, &ae1f.m_nodes[0], ae1f.m_nodes.size(), 0);
        {
          std::ostringstream o;
          o << "bitcur_" << i << "_" << b;
          medith::write(o.str().c_str(), bits[b].m_crd, &bits[b].m_nodes[0], bits[b].m_nodes.size(), bits[b].shift());
        }
#endif

        //ae1f.reverse_orient();

        tmpbits.clear();
        bool true_clip(false);
        int err = NUGA::CLIP::isolated_clip<aelt_t, aelt_t>(bits[b], ae1f, NUGA::INTERSECT::DIFFERENCE, ARTOL, tmpbits, true_clip);

        if (!true_clip) continue;

        // COMPRESS STRATEGY (with MOVE SEMANTICS) for bits => b can be reduced of 1 to treat the replaced at next iter 
        // if tmpbits is empty (IN) => compress (erase bits if single, put the last instead of current otherwise)
        // else replace the first bit, append the other bits  . 
        NUGA::CLIP::__replace_append_and_next_iter(bits, b, tmpbits);
        if ((E_Int)bits.size() != nbits) nbits = bits.size(); //update iff the last have replaced the current 
      }
    }

    if (!bits.empty())
    {
      for (size_t b = 0; b < bits.size(); ++b)
      {
        if (bits[b].empty()) continue;
        //fixme : m_poids might be inconsitent or inexistent for the residual bit => so force to be in 'coordinate appending' mode in add
        xm.add(bits[b], true/*append vertices*/);
        anc0.push_back(i);

#ifdef SUPERMESH_DBG
        std::ostringstream o;
        o << "remain_bit_" << i << "_" << b;
        medith::write(o.str().c_str(), bits[b].m_crd, &bits[b].m_nodes[0], bits[b].m_nodes.size(), bits[b].shift());
#endif
      }
      //std::cout << i << " has a remaining bit" << std::endl;
    }
  }

  xm.cnt.updateFacets();

#ifdef SUPERMESH_DBG
  K_FLD::IntArray cnto;
  ngon_type ngo(xm.cnt, false);
  ngo.export_to_array(cnto);
  tp::write("D:\\slandier\\DATA\\tmp\\interpol\\SURF\\tmp\\xm.tp", xm.crd, cnto, "NGON");
#endif

}

template <typename zmesh_t>
void interpol_coeffs_for_first(
  const zmesh_t& m0,
  const zmesh_t& m1, 
  double RTOL,
  std::vector<int>& dindices,
  std::vector<double>& dcoeffs,
  std::vector<int>& xr, bool do_omp=false)
{
  dindices.clear();
  dcoeffs.clear();
  xr.clear();

  using aelt_t = typename zmesh_t::aelt_t;

  auto loc1 = m1.get_localizer();

  m0.get_nodal_metric2();
  m1.get_nodal_metric2();

  E_Int nbcells0 = m0.ncells();
  std::vector<E_Int> cands;

  std::vector<aelt_t> bits;

  std::vector<std::vector<int>>    dindices_per_cell(nbcells0);
  std::vector<std::vector<double>> dcoeffs_per_cell(nbcells0);
  
  size_t n, k;
  E_Int i, i2;

#pragma omp parallel for private(cands, bits, n, k, i, i2) if(do_omp)
  for (i = 0; i < nbcells0; ++i)
  {
    auto ae0 = m0.aelement(i);
    
    cands.clear();
    loc1->get_candidates(ae0, ae0.m_crd, cands, 1, RTOL); //return as 1-based
    if (cands.empty()) continue;

    for (n = 0; n < cands.size(); ++n)
    {
      i2 = cands[n] - 1;
      auto ae1 = m1.aelement(i2);
      bits.clear();
      bool true_clip = false;
      NUGA::CLIP::isolated_clip<aelt_t, aelt_t>(ae0, ae1, NUGA::INTERSECT::INTERSECTION, RTOL, bits, true_clip);

      if (bits.empty()) continue;

      double si2 = 0.;
      for (k = 0; k < bits.size(); ++k)
        si2 += bits[k].extent();

      dindices_per_cell[i].push_back(i2);
      dcoeffs_per_cell[i].push_back(si2);
    }
  }

  // normalize coeffs with total covered surface
  for (size_t i = 0; i < nbcells0; ++i)
  {
    double stot = 0.;
    for (size_t k = 0; k < dcoeffs_per_cell[i].size(); ++k)
      stot += dcoeffs_per_cell[i][k];

    //double s0 = m0.aelement(i).extent();
    //if (stot / s0 < 0.99) std::cout << "cell " << i << " has " << stot / s0 << std::endl;

    for (size_t k = 0; k < dcoeffs_per_cell[i].size(); ++k)
      dcoeffs_per_cell[i][k] /= stot;
  }

  // concatenate info for exit
  xr.push_back(0);
  for (size_t i = 0; i < nbcells0; ++i)
  {
    dindices.insert(dindices.end(), ALL(dindices_per_cell[i]));
    dcoeffs.insert(dcoeffs.end(), ALL(dcoeffs_per_cell[i]));

    xr.push_back(dindices.size());

  }


}

inline double transfer_mass(aPolyhedron<0>& pbit, aPolyhedron<0>& pdon, double Vbit, double rho, 
                            double* gradx = nullptr, double* grady = nullptr, double* gradz = nullptr)
{
  double m = Vbit * rho;

  if (gradx == nullptr) // order1
    return m;

  double gradf[] = { *gradx, *grady, *gradz};
  const double* Gbit = pbit.get_centroid();
  const double* Gdon = pdon.get_centroid();
  double GdonGbit[3];
  NUGA::diff<3>(Gbit, Gdon, GdonGbit);

  m += NUGA::dot<3>(GdonGbit, gradf) * Vbit;
  return m;
}

template <typename zmesh_t>
int interpolate(
  const zmesh_t& mrec,
  const zmesh_t& mdon,
  double RTOL,
  const std::vector<field> & don_fields,
  std::vector<std::vector<double>>& rec_fields,
  bool do_omp = false)
{
  rec_fields.clear();

  int nfields = don_fields.size();

  rec_fields.resize(nfields);
  for (size_t i = 0; i < don_fields.size(); ++i)
    rec_fields[i].resize(mrec.ncells(), 0.);

  std::vector<double> covered_vol(mrec.ncells(), 0.);

  //if (mrec.oriented == 0) return;
  //if (mdon.oriented == 0) return;

  using aelt_t = typename zmesh_t::aelt_t;

  auto loc1 = mdon.get_localizer();

  mrec.get_nodal_metric2();
  mdon.get_nodal_metric2();

  E_Int nbcells0 = mrec.ncells();
  std::vector<E_Int> cands;
  std::vector<aelt_t> bits;

  size_t n, k, f, b;
  E_Int i, i2;

#pragma omp parallel for private(cands, bits, n, k, i, i2, f, b) if(do_omp)
  for (i = 0; i < nbcells0; ++i)
  {
    auto ae0 = mrec.aelement(i);

    cands.clear();
    loc1->get_candidates(ae0, ae0.m_crd, cands, 1, RTOL); //return as 1-based
    if (cands.empty()) continue;

#ifdef SUPERMESH_DBG
    medith::write<>("ae0", ae0.m_crd, ae0.m_pgs);
#endif

    for (n = 0; n < cands.size(); ++n)
    {
      i2 = cands[n] - 1;
      auto ae1 = mdon.aelement(i2);


#ifdef SUPERMESH_DBG
      medith::write<>("ae1", ae1.m_crd, ae1.m_pgs);
#endif

      bits.clear();
      bool just_io = !NUGA::CLIP::compute(ae0, ae1, NUGA::INTERSECT::INTERSECTION, bits); //robust_clip returns true if true clip
                                                                                         // IO : current bit does not intersect front (or only along its boundaries)
      if (just_io)
      {
        assert(bits.empty());
        NUGA::eClassify wher = NUGA::CLASSIFY::classify(ae0, ae1, true);

        if (wher == OUT) continue;
        
        if (wher == IN) // ae0 is fully inside ae1
        {
          double Vbit = ae0.extent();
          // transfer fields
          for (f = 0; f < nfields; ++f)
          {
            double* gradf[] = { nullptr, nullptr, nullptr };
            if (don_fields[f].gradf[0] != nullptr)//this field has gradients
            {
              gradf[0] = &don_fields[f].gradf[0][i2];
              gradf[1] = &don_fields[f].gradf[1][i2];
              gradf[2] = &don_fields[f].gradf[2][i2];
            }

            rec_fields[f][i] += transfer_mass(ae0, ae1, Vbit, don_fields[f].f[i2], gradf[0], gradf[1], gradf[2]);
              covered_vol[i] += Vbit;
          }
          break; // fully IN, so stop checking other candidates
        }
        else // IN_1 : ae1 is fully inside ae0
        {
          //gradients doesnt matter here since Vbit = Vdonnor
          double Vbit = ae1.extent();
          // transfer fields
          for (f = 0; f < nfields; ++f)
          {
            rec_fields[f][i] += transfer_mass(ae0, ae1, Vbit, don_fields[f].f[i2]);
            covered_vol[i] += Vbit;
          }
        }
      }
      else // true clip
      {
        for (b = 0; b < bits.size(); ++b)
        {
          double Vbit = bits[b].extent();
          // transfer fields
          for (f = 0; f < nfields; ++f)
          {
            double* gradf[] = { nullptr, nullptr, nullptr };
            if (don_fields[f].gradf[0] != nullptr)//this field has gradients
            {
              gradf[0] = &don_fields[f].gradf[0][i2];
              gradf[1] = &don_fields[f].gradf[1][i2];
              gradf[2] = &don_fields[f].gradf[2][i2];
            }

            rec_fields[f][i] += transfer_mass(bits[b], ae1, Vbit, don_fields[f].f[i2], gradf[0], gradf[1], gradf[2]);
            covered_vol[i] += Vbit;
          }
        }
      }
    }
  }
  
#pragma omp parallel for private(i) if (do_omp)
  for (i = 0; i < nbcells0; ++i)
  {
    for (f = 0; f < nfields; ++f)
    {
      if (covered_vol[i] != 0.) // discard non interpolated cells
        rec_fields[f][i] /= covered_vol[i];
    }
  }
  return 0;
}


}

#endif // NUGA_SUPERMESHHXX
