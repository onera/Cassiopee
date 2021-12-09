/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __NUGA_SYNC_FACES_HXX__
#define __NUGA_SYNC_FACES_HXX__

#include "Nuga/include/mesh_t.hxx"

// match perio stuff ////
using ngon_type = ngon_t<K_FLD::IntArray>;
using crd_t = K_FLD::FloatArray;
using acrd_t = K_FLD::ArrayAccessor<crd_t>;

///
void detect_async_modified_faces(NUGA::ph_mesh_t& vmesh, const double* center, const double * axis, double angle, const double* translation,
  double ARTOL, std::map<int, std::vector<int>>& glob_face_to_bits)
{
  DELAUNAY::Triangulator dt;

  NUGA::pg_smesh_t m;
  std::vector<int> oids, ancestors;

  vmesh.get_boundary(m, oids, ancestors);

  std::map<int, std::vector<int>> loc_face_to_bits;

  glob_face_to_bits.clear();

  // 0. TOL computation
  m.get_nodal_metric2();

  if (ARTOL == 0.) ARTOL = -0.01;

  double TOL = ARTOL;
  if (ARTOL < 0.) // relative
  {
    E_Float Lmin, Lmax;
    E_Int imin, imax;
    ngon_type::edge_length_extrema(m.cnt, m.crd, Lmin, imin, Lmax, imax);
    TOL = -ARTOL * Lmin;
    assert(TOL > 0.); // no degen in m
  }

  E_Int npgs = m.cnt.size();

  // 1. compute centroids
  K_FLD::FloatArray Cs(3, npgs);
  for (E_Int i = 0; i < npgs; ++i)
  {
    const E_Int* nodes = m.cnt.get_facets_ptr(i);
    E_Int nb_nodes = m.cnt.stride(i);
    //K_MESH::Polygon::simplified_centroid(crd, nodes, nb_nodes, 1, Cs.col(i));
    //K_MESH::Polygon::centroid<3>(crd, nodes, nb_nodes, 1, Cs.col(i));
    K_MESH::Polygon::iso_barycenter(m.crd, nodes, nb_nodes, 1, Cs.col(i));
  }

  // 2. apply transfo
  K_FLD::FloatArray CsR(Cs);
  if (axis)
    NUGA::axial_rotate(CsR, center, axis, angle);
  else if (translation){
    for (size_t k = 0; k < CsR.cols(); ++k)
    {
      CsR(0, k) += translation[0];
      CsR(1, k) += translation[1];
      CsR(2, k) += translation[2];
    }
  }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }

  ////////////////////////////////////////////////////////////////////////////

  acrd_t aCs(Cs);
  K_SEARCH::KdTree<> tree(aCs, EPSILON);

  std::vector<E_Int> nids(npgs, IDX_NONE);
  std::vector<double> d2s(npgs, IDX_NONE);

  ////////////////////////////////////////////////////////////////////////////

  // 3. append both
  K_FLD::FloatArray Cs_glog{ Cs };
  Cs_glog.pushBack(CsR);

#ifdef DEBUG_MATCH_PERIO
  {
    K_FLD::IntArray cntt(2, 1, 0);
    cntt(1, 0) = 1;
    medith::write("cloud", Cs_glog, cntt, "BAR");
  }
#endif

  // 4. Filtrate exact-matches

  K_FLD::ArrayAccessor<K_FLD::FloatArray> ca(Cs_glog);
  E_Int nmerges = ::merge(ca, TOL, nids);
  //std::cout << "nmerges : " << nmerges << std::endl;

  // check wrong auto_match
  for (size_t i = 0; i <npgs; ++i)
    assert(nids[i] == i);

  std::vector<E_Int> left, right, remain;
  for (size_t i = npgs; i <nids.size(); ++i)
  {
    if (nids[i] == i)
    {
      //remain.push_back(i-npgs);
      continue;
    }

    left.push_back(i - npgs);
    right.push_back(nids[i]);
    E_Int nid = nids[i];
    nids[i - npgs] = nid;
    nids[nid] = i - npgs;
  }

  nids.resize(npgs);

#ifdef DEBUG_MATCH_PERIO
  medith::write("left", m.crd, m.cnt, &left, 0);
  medith::write("right", m.crd, m.cnt, &right, 0);
  //medith::write("remain", crd, PGS, &remain, 0);
#endif

  // 5. Filtrate non-perio by boxes
  auto loc = m.get_localizer();
  std::vector<int> cands;

  K_FLD::FloatArray crdR(m.crd);
  if (axis)
    NUGA::axial_rotate(crdR, center, axis, angle);
  else if (translation)
    for (size_t k = 0; k < crdR.cols(); ++k)
    {
      crdR(0, k) += translation[0];
      crdR(1, k) += translation[1];
      crdR(2, k) += translation[2];
    }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }

  acrd_t acrdR(crdR);
  acrd_t acrd(m.crd);

  for (size_t i = 0; i < npgs; ++i)
  {
    if (nids[i] != i) continue; // exact match already found

    auto face = m.element(i);

    K_SEARCH::BBox3D b;
    face.bbox(crdR, b); // compute box in rotated frame

    loc->get_candidates(b, cands);
    if (cands.empty()) continue;

    /*{
      K_FLD::FloatArray cc;
      K_FLD::IntArray ccnt;
      ngon_type ng;
      b.convert2NG(cc, ng);
      ng.export_to_array(ccnt);
      tp::write("D:\\slandier\\DATA\\tmp\\match_perio\\box.tp", cc, ccnt, "NGON");

      medith::write("face", crdR, m.cnt, i);
    }*/

    double n[3], nc[3];
    face.normal<K_FLD::FloatArray, 3>(crdR, n);

    double G[3];
    face.iso_barycenter<acrd_t, 3>(acrdR, G);

    for (size_t c = 0; c< cands.size(); ++c)
    {
      E_Int cand = cands[c];

      auto faceC = m.element(cand);
      faceC.normal<K_FLD::FloatArray, 3>(m.crd, nc);

      // true candidates must be colinear and opposedly oriented, so ps should near -1
      double ps = NUGA::dot<3>(n, nc);
      if (ps > -0.99) continue;

      // if G is in faceC, face is a piece of faceC
      bool is_in1;
      int err = faceC.fast_is_in_pred<DELAUNAY::Triangulator>(dt, m.crd, G, is_in1);

      //reciprocal test
      double GC[3];
      faceC.iso_barycenter<acrd_t, 3>(acrd, GC);
      bool  is_in2;
      err = face.fast_is_in_pred<DELAUNAY::Triangulator>(dt, crdR, GC, is_in2);

      if (!is_in1 && !is_in2) continue;

      if (is_in1 && !is_in2)
      {
        loc_face_to_bits[-(cand+1)].push_back(i); //neg storing to apply -Transfo
        continue;
      }

      if (is_in2 && !is_in1)
      {
        loc_face_to_bits[i].push_back(-(cand + 1)); //pos storing for i to apply +Transfo
        continue;
      }

      // if mutual inclusion : surface test
      double s1 = face.surface<3>(crdR);
      double s2 = faceC.surface<3>(m.crd);

      if (s1 < s2)
      {
        loc_face_to_bits[-(cand + 1)].push_back(i);
        continue;
      }
      else
      {
        loc_face_to_bits[i].push_back(-(cand + 1));
        continue;
      }
    }
  }

  // FILTER : ONLY ONE FACE PER PH PER ITER
  std::map<int, std::vector<int>> tmp;
  std::set<int> involved_phs;

  for (auto ii : loc_face_to_bits)
  {
    int id = ii.first;
    if (id < 0) id = -(id + 1);

    if (!involved_phs.insert(ancestors[id]).second) // already in
      continue;

    tmp.insert(std::make_pair(ii.first, ii.second));
  }

  loc_face_to_bits = tmp;

  // indirection to refer to volume mesh ids
  for (auto i : loc_face_to_bits)
  {
    int lface = i.first;
    int gface = 0;
    if (lface < 0)
    {
      lface = -(lface + 1);
      gface = oids[lface];
      gface = -(gface + 1);
    }
    else
      gface = oids[lface];

    auto gbits = i.second; //copy
    for (auto& k : gbits)
    {
      if (k < 0)
      {
        k = oids[-k - 1];
        k = -(k + 1);
      }
      else
        k = oids[k];
    }
    glob_face_to_bits[gface] = gbits;
  }

  //make it reciprocal (ONLY IF bits belong to the same cell)
  std::map<int, std::vector<int>> tmp1;
  for (auto i : glob_face_to_bits)
  {
    int master = i.first;
    auto& bits = i.second;
    tmp1[master] = bits;

    bool ancPH = ancestors[bits[0]];
    bool samePH = true;
    for (size_t k = 0; (k < bits.size() - 1) && samePH; ++k)
    {
      int bk = bits[k];
      int bkp1 = bits[k + 1];
      bk = (bk < 0) ? -(bk + 1) : bk;
      bkp1 = (bkp1 < 0) ? -(bkp1 + 1) : bkp1;

      samePH = (ancestors[bk] == ancestors[bkp1]);
    }
    if (!samePH) continue;

    for (auto& b : bits)tmp1[b].push_back(master);
  }
  glob_face_to_bits = tmp1;
}

///
void duplicate_and_move_period_faces
(NUGA::ph_mesh_t& m, const double* center, const double * axis, double angle, const double* translation,
 std::map<int, std::vector<int>>& face_to_bits)
{
  std::set<int> leftF; //apply +Transfo
  std::set<int> rightF;//apply -Transfo

  // 1. update ids (i.e. remove neg vals), separate ids in these two groups to apply appropriate transformation
  std::map<int, std::vector<int>> new_face_to_bits;
  for (auto i : face_to_bits)
  {
    int master = i.first;
    auto& bits = i.second;

    if (master < 0) {
      master = -master - 1;
      rightF.insert(master);
    }
    else
      leftF.insert(master);

    for (auto& b : bits)
    {
      if (b < 0) {
        b = -b - 1;
        rightF.insert(b);
      }
      else
        leftF.insert(b);
    }

    new_face_to_bits[master] = bits;                  // no more neg vals
  }

  face_to_bits = new_face_to_bits;//no more neg vals

  //2. extract these faces and apply Transfo

  std::vector<int> oids1;
  ngon_unit pgsL;
  std::vector<E_Int> lF(ALL(leftF));
  m.cnt.PGs.extract(lF, pgsL, oids1);

  //modify history such moved bits receive target location history (no sense for split master that will vanish anyway)
  //reset first
  pgsL._ancEs.clear();
  pgsL._ancEs.resize(2, pgsL.size(), IDX_NONE);
  K_CONNECT::IdTool::init_inc(pgsL._ancEs, 0, pgsL.size());
  int k = -1;
  for (auto i : lF)
  {
    ++k;
    auto tgtid = face_to_bits.find(i);
    assert(tgtid != face_to_bits.end());
    if (tgtid->second.size() != 1) continue; //split master
    pgsL._ancEs(0, k) = tgtid->second[0];
  }


  K_FLD::FloatArray crdL(m.crd);
  ngon_type::compact_to_used_nodes(pgsL, crdL);

  // reverse orient
  /*{
    NUGA::pg_smesh_t toto;
    toto.crd = crdL;
    toto.cnt = pgsL;
    toto.oriented = 0;
    toto.reverse_orient();
    pgsL = toto.cnt;
  }*/

  //
  if (axis)
    NUGA::axial_rotate(crdL, center, axis, angle);
  else if (translation)
    for (size_t k = 0; k < crdL.cols(); ++k)
    {
      crdL(0, k) += translation[0];
      crdL(1, k) += translation[1];
      crdL(2, k) += translation[2];
    }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }


#ifdef DEBUG_MATCH_PERIO
  {
    ngon_type ngo(pgsL, false);
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tp::write("D:\\slandier\\DATA\\tmp\\match_perio\\pgsL.tp", crdL, cnto, "NGON");
  }
#endif

  std::vector<int> oids2;
  ngon_unit pgsR;
  std::vector<E_Int> rF(ALL(rightF));
  m.cnt.PGs.extract(rF, pgsR, oids2);

  //modify history such moved bits receive target location history (no sense for split master that will vanish anyway)
  //reset first
  pgsR._ancEs.clear();
  pgsR._ancEs.resize(2, pgsR.size(), IDX_NONE);
  K_CONNECT::IdTool::init_inc(pgsR._ancEs, 0, pgsR.size());
  k = -1;
  for (auto i : rF)
  {
    ++k;
    auto tgtid = face_to_bits.find(i);
    assert(tgtid != face_to_bits.end());
    if (tgtid->second.size() != 1) continue; //split master
    pgsR._ancEs(0, k) = tgtid->second[0];
  }

  K_FLD::FloatArray crdR(m.crd);
  ngon_type::compact_to_used_nodes(pgsR, crdR);

  // reverse orient
  /*{
    NUGA::pg_smesh_t toto;
    toto.crd = crdR;
    toto.cnt = pgsR;
    toto.oriented = 0;
    toto.reverse_orient();
    pgsR = toto.cnt;
  }*/

  //
  if (axis)
    NUGA::axial_rotate(crdR, center, axis, -angle);
  else if (translation)
    for (size_t k = 0; k < crdR.cols(); ++k)
    {
      crdR(0, k) -= translation[0];
      crdR(1, k) -= translation[1];
      crdR(2, k) -= translation[2];
    }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }

#ifdef DEBUG_MATCH_PERIO
  {
    ngon_type ngo(pgsR, false);
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tp::write("D:\\slandier\\DATA\\tmp\\match_perio\\pgsR.tp", crdR, cnto, "NGON");
  }
#endif

  //3. append these faces (geometrically, with their own coordinates)

  std::vector<int> nids;
  K_CONNECT::IdTool::init_inc(nids, m.cnt.PGs.size());

  E_Int shift = m.cnt.PGs.size();

  // append L
  int pt_shift = m.crd.cols();
  pgsL.shift(pt_shift);
  m.crd.pushBack(crdL);
  m.cnt.PGs.append(pgsL);

  for (size_t i = 0; i < oids1.size(); ++i)
    nids[oids1[i]] = i + shift;

  shift = m.cnt.PGs.size();

  // append R
  pt_shift = m.crd.cols();
  pgsR.shift(pt_shift);
  m.crd.pushBack(crdR);
  m.cnt.PGs.append(pgsR);

  for (size_t i = 0; i < oids2.size(); ++i)
    nids[oids2[i]] = i + shift;

  //update face_to_bits
  for (auto& i : face_to_bits)
  {
    int master = i.first;
    auto& bits = i.second;
    for (auto& b : bits) b = nids[b];
  }
}

///
void sync_faces
(NUGA::ph_mesh_t& m, const std::map<int, std::vector<int>>& face_to_bits, double ARTOL)
{
  // 
  std::vector<int> molecPH;
  ngon_unit new_phs;

  // 1. replace faces and keep track of modified PHs
  std::set<E_Int> modifiedPHs;
  for (E_Int i = 0; i < m.cnt.PHs.size(); ++i)
  {
    const E_Int* pPGi = m.cnt.PHs.get_facets_ptr(i);
    E_Int nb_pgs = m.cnt.PHs.stride(i);

    molecPH.clear();

    for (E_Int p = 0; p < nb_pgs; ++p)
    {
      E_Int PGi = *(pPGi + p) - 1;

      auto itBits = face_to_bits.find(PGi);
      if (itBits == face_to_bits.end())
        molecPH.push_back(PGi + 1);
      else
      {
        modifiedPHs.insert(i);
        for (E_Int k = 0; k < itBits->second.size(); ++k)
          molecPH.push_back(itBits->second[k] + 1);
      }
    }

    new_phs.add(molecPH.size(), &molecPH[0]);
  }

  new_phs._type = m.cnt.PHs._type;  // hack to preserve flags (externality)
  new_phs._ancEs = m.cnt.PHs._ancEs;// hack

  m.cnt.PHs = new_phs;
  m.cnt.PHs.updateFacets();
  m.cnt.PHs.remove_duplicated(); //several occurence of the same face in each phs

                                 // 2. merge coincident nodes
  if (ARTOL == 0.) ARTOL = -0.01;

  double TOL = ARTOL;
  if (ARTOL < 0.) // relative
  {
    E_Float Lmin, Lmax;
    E_Int imin, imax;
    ngon_type::edge_length_extrema(m.cnt.PGs, m.crd, Lmin, imin, Lmax, imax);
    TOL = -ARTOL * Lmin;
    assert(TOL > 0.); // no degen in m
  }

  //std::cout << "TOL : " << TOL << std::endl;
  E_Int nb_merges = m.cnt.join_phs(m.crd, TOL);
  //std::cout << "nmerges : " << nb_merges << std::endl;

  //3. close_phs
  std::vector<E_Int> modPHs(ALL(modifiedPHs));
  ngon_type::close_phs(m.cnt, m.crd, &modPHs);

  //4. replace moved master by modified bits : i.e replace original bits by their modified version
  //4.a reverse face_to_bits
  std::map<int, std::vector<int>> face_to_bits_rev, tmp;
  for (auto i : face_to_bits)
  {
    auto bits = i.second;
    int f = i.first;
    for (auto k : bits)
      face_to_bits_rev[k].push_back(f);
  }
  //4.b keep only releavnt
  for (auto i : face_to_bits_rev)
  {
    if (i.second.size() == 1) continue;
    tmp[i.first] = i.second;
  }
  face_to_bits_rev = tmp;

  //4.c replace back
  if (!face_to_bits_rev.empty())
  {
    new_phs.clear();
    for (E_Int i = 0; i < m.cnt.PHs.size(); ++i)
    {
      const E_Int* pPGi = m.cnt.PHs.get_facets_ptr(i);
      E_Int nb_pgs = m.cnt.PHs.stride(i);

      molecPH.clear();

      for (E_Int p = 0; p < nb_pgs; ++p)
      {
        E_Int PGi = *(pPGi + p) - 1;

        auto itBits = face_to_bits_rev.find(PGi);
        if (itBits == face_to_bits_rev.end())
          molecPH.push_back(PGi + 1);
        else
        {
          for (E_Int k = 0; k < itBits->second.size(); ++k)
            molecPH.push_back(itBits->second[k] + 1);
        }
      }

      new_phs.add(molecPH.size(), &molecPH[0]);
    }

    new_phs._type = m.cnt.PHs._type;  // hack to preserve flags (externality)
    new_phs._ancEs = m.cnt.PHs._ancEs;// hack

    m.cnt.PHs = new_phs;
    m.cnt.PHs.updateFacets();

    ngon_type::close_phs(m.cnt, m.crd, &modPHs);
  }

  //5. clean
  std::vector<E_Int> pgnids, phnids;
  m.cnt.remove_unreferenced_pgs(pgnids, phnids);

  // assert closed
}

#endif
