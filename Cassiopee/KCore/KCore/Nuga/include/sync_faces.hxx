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
inline void detect_async_modified_faces(NUGA::ph_mesh_t& vmesh, const E_Float* center, const E_Float* axis, E_Float angle, const E_Float* translation,
  E_Float ARTOL, std::map<E_Int, std::vector<E_Int>>& glob_face_to_bits)
{
  DELAUNAY::Triangulator dt;

  NUGA::pg_smesh_t m;
  std::vector<E_Int> oids, ancestors;

  vmesh.get_boundary(m, oids, ancestors);

  std::map<E_Int, std::vector<E_Int>> loc_face_to_bits;

  glob_face_to_bits.clear();

  // 0. TOL computation
  m.get_nodal_metric2();

  if (ARTOL == 0.) ARTOL = -0.01;

  E_Float TOL = ARTOL;
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
  else if (translation)
  {
    for (E_Int k = 0; k < CsR.cols(); ++k)
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
  std::vector<E_Float> d2s(npgs, IDX_NONE);

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
  //E_Int nmerges = 
  ::merge(ca, TOL, nids);
  //std::cout << "nmerges : " << nmerges << std::endl;

  // check wrong auto_match
  for (E_Int i = 0; i < npgs; ++i)
    assert(nids[i] == i);

#ifdef DEBUG_MATCH_PERIO
  std::vector<E_Int> left, right, remain;
#endif
  for (size_t i = npgs; i < nids.size(); ++i)
  {
    if (nids[i] == E_Int(i))
    {
#ifdef DEBUG_MATCH_PERIO
      remain.push_back(i - npgs);
#endif
      continue;
    }

#ifdef DEBUG_MATCH_PERIO
    left.push_back(i - npgs);
    right.push_back(nids[i]);
#endif

    E_Int nid = nids[i];
    nids[i - npgs] = nid;
    nids[nid] = i - npgs;
  }

  nids.resize(npgs);

#ifdef DEBUG_MATCH_PERIO
  //medith::write("left", m.crd, m.cnt, &left, 0);
  //medith::write("right", m.crd, m.cnt, &right, 0);
  //medith::write("remain", m.crd, m.cnt, &remain, 0);
#endif

  // 5. Filtrate non-perio by boxes
  auto loc = m.get_localizer();
  std::vector<E_Int> cands;

  K_FLD::FloatArray crdR(m.crd);
  if (axis)
    NUGA::axial_rotate(crdR, center, axis, angle);
  else if (translation)
  {
    for (E_Int k = 0; k < crdR.cols(); ++k)
    {
      crdR(0, k) += translation[0];
      crdR(1, k) += translation[1];
      crdR(2, k) += translation[2];
    }
  }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }

  acrd_t acrdR(crdR);
  acrd_t acrd(m.crd);

  for (E_Int i = 0; i < npgs; ++i)
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

    E_Float n[3], nc[3];
    face.normal<K_FLD::FloatArray, 3>(crdR, n);

    E_Float G[3];
    face.iso_barycenter<acrd_t, 3>(acrdR, G);
    E_Float s1 = face.surface<3>(crdR);

    for (size_t c = 0; c< cands.size(); ++c)
    {
      E_Int cand = cands[c];

      auto faceC = m.element(cand);
      faceC.normal<K_FLD::FloatArray, 3>(m.crd, nc);

      // true candidates must be colinear and opposedly oriented, so ps should near -1
      E_Float ps = NUGA::dot<3>(n, nc);
      if (ps > -0.99) continue;

      // check for fully-in case : is face fully inside faceC ?
      bool is_in2{ false };
      E_Float RTOL = 1.e-6;
      for (E_Int k = 0; k < face.nb_nodes(); ++k)
      {
        const E_Float* P = crdR.col(face.node(k));
        // if P is in ae1, ae0 is a piece of ae1
        faceC.fast_is_in_pred<DELAUNAY::Triangulator, 3>(dt, m.crd, P, is_in2, RTOL);
        if (!is_in2) break;
      }

      //reciprocal test : is faceC fully inside face ?
      bool is_in1{ false };
      for (E_Int k = 0; k < faceC.nb_nodes(); ++k)
      {
        const E_Float* P = m.crd.col(faceC.node(k));
        // if P is in ae1, ae0 is a piece of ae1
        face.fast_is_in_pred<DELAUNAY::Triangulator, 3>(dt, crdR, P, is_in1, RTOL);
        if (!is_in1) break;
      }

      if (!is_in1 && !is_in2) continue;

      if (is_in1 && !is_in2)
      {
        /*std::ostringstream o;
        o << i << "_master";
        medith::write(o.str().c_str(), crdR, m.cnt, i);
        o.str("");
        o << i << "_bit_" << cand;
        medith::write(o.str().c_str(), m.crd, m.cnt, cand);*/

        loc_face_to_bits[i].push_back(-(cand + 1)); //pos storing for i to apply +Transfo
        continue;
      }
      else if (is_in2 && !is_in1)
      {
        /*std::ostringstream o;
        o << cand << "_master";
        medith::write(o.str().c_str(), crdR, m.cnt, cand);
        o.str("");
        o << cand << "_bit_" << i;
        medith::write(o.str().c_str(), m.crd, m.cnt, i);*/

        loc_face_to_bits[-(cand + 1)].push_back(i); //neg storing to apply -Transfo
        continue;
      }

      // if mutual inclusion : surface test

      E_Float s2 = faceC.surface<3>(m.crd);

      if (s1 < s2)
      {
        /*std::ostringstream o;
        o << i << "_same";
        medith::write(o.str().c_str(), crdR, m.cnt, i);
        o.str("");
        o << i << "_same2";
        medith::write(o.str().c_str(), m.crd, m.cnt, cand);*/

        loc_face_to_bits[-(cand + 1)].push_back(i);
        continue;
      }
      else
      {
        /*std::ostringstream o;
        o << i << "_same";
        medith::write(o.str().c_str(), crdR, m.cnt, i);
        o.str("");
        o << i << "_same2";
        medith::write(o.str().c_str(), m.crd, m.cnt, cand);*/

        loc_face_to_bits[i].push_back(-(cand + 1));
        continue;
      }
    }
  }

  // LIMITATION : the following is here to avoid mesh corruption
  // It means that input is not processable by this algo (some agglomeration might have occured that breaks the pairing
  // check that surfaces are matching

  std::map<E_Int, std::vector<E_Int>> tmp;

  for (auto i : loc_face_to_bits)
  {
    E_Int master = (i.first < 0) ? -(i.first+1) : i.first;

    auto face = m.element(master);
    E_Float s1 = face.surface<3>(crdR);
    E_Float s2{ 0. };
    for (size_t k = 0; k < i.second.size(); ++k)
    {
      E_Int bitid = (i.second[k] < 0) ? -(i.second[k] + 1) : i.second[k];
      auto fb = m.element(bitid);
      s2 += fb.surface<3>(m.crd);
    }

    E_Float err = ::fabs((s1 - s2) / s2);
    if (err < 1.e-2)
      tmp.insert(std::make_pair(i.first, i.second));
  }

  loc_face_to_bits = tmp;

  // FILTER : ONLY ONE FACE PER PH PER ITER
  std::set<E_Int> involved_phs;

  tmp.clear();
  for (auto ii : loc_face_to_bits)
  {
    E_Int id = ii.first;
    if (id < 0) id = -(id + 1);

    if (!involved_phs.insert(ancestors[id]).second) // already in
      continue;

    tmp.insert(std::make_pair(ii.first, ii.second));
  }

  loc_face_to_bits = tmp;

  //make it reciprocal (ONLY IF bits belong to the same cell)
  std::map<E_Int, std::vector<E_Int>> tmp1;
  for (auto i : loc_face_to_bits)
  {
    E_Int master = i.first;
    auto& bits = i.second;
    tmp1[master] = bits;

    bool samePH = true;
    for (size_t k = 0; (k < bits.size() - 1) && samePH; ++k)
    {
      E_Int bk = bits[k];
      E_Int bkp1 = bits[k + 1];
      bk = (bk < 0) ? -(bk + 1) : bk;
      bkp1 = (bkp1 < 0) ? -(bkp1 + 1) : bkp1;

      samePH = (ancestors[bk] == ancestors[bkp1]);
    }
    if (!samePH) continue;

    for (auto& b : bits)tmp1[b].push_back(master);
  }
  loc_face_to_bits = tmp1;


  // indirection to refer to volume mesh ids
  for (auto i : loc_face_to_bits)
  {
    E_Int lface = i.first;
    E_Int gface = 0;
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

}

///
inline void duplicate_and_move_period_faces
(NUGA::ph_mesh_t& m, const E_Float* center, const E_Float* axis, E_Float angle, const E_Float* translation,
 std::map<E_Int, std::vector<E_Int>>& face_to_bits)
{
  std::set<E_Int> leftF; //apply +Transfo
  std::set<E_Int> rightF;//apply -Transfo

  // 1. update ids (i.e. remove neg vals), separate ids in these two groups to apply appropriate transformation
  std::map<E_Int, std::vector<E_Int>> new_face_to_bits;
  for (auto i : face_to_bits)
  {
    E_Int master = i.first;
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

  std::vector<E_Int> oids1;
  ngon_unit pgsL;
  std::vector<E_Int> lF(ALL(leftF));
  m.cnt.PGs.extract(lF, pgsL, oids1);

  //modify history such moved bits receive target location history (no sense for split master that will vanish anyway)
  //reset first
  pgsL._ancEs.clear();
  pgsL._ancEs.resize(2, pgsL.size(), IDX_NONE);
  K_CONNECT::IdTool::init_inc(pgsL._ancEs, 0, pgsL.size());
  E_Int k = -1;
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
    for (E_Int k = 0; k < crdL.cols(); ++k)
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

  std::vector<E_Int> oids2;
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
  {
    for (E_Int k = 0; k < crdR.cols(); ++k)
    {
      crdR(0, k) -= translation[0];
      crdR(1, k) -= translation[1];
      crdR(2, k) -= translation[2];
    }
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

  std::vector<E_Int> nids;
  K_CONNECT::IdTool::init_inc(nids, m.cnt.PGs.size());

  E_Int shift = m.cnt.PGs.size();

  // append L
  E_Int pt_shift = m.crd.cols();
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
    //E_Int master = i.first;
    auto& bits = i.second;
    for (auto& b : bits) b = nids[b];
  }
}

///
inline void replace_faces(
                  NUGA::ph_mesh_t& m,
                  const std::map<E_Int,
                  std::vector<E_Int>>& face_to_bits,
                  E_Float TOL,
                  std::set<E_Int>& modifiedPHs)
{
  std::vector<E_Int> molecPH;
  ngon_unit new_phs;
  for (E_Int i = 0; i < m.cnt.PHs.size(); ++i)
  {
    const E_Int* pPGi = m.cnt.PHs.get_facets_ptr(i);
    int nb_pgs = m.cnt.PHs.stride(i);

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
        for (size_t k = 0; k < itBits->second.size(); ++k)
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
  {
    std::vector<E_Int> nids, lnids;
    K_CONNECT::IdTool::init_inc(nids, m.crd.cols());

    for (auto i : modifiedPHs)
    {
      K_MESH::Polyhedron<0> ph(m.cnt, i);
      NUGA::haPolyhedron<0> hph;
      hph.set(ph, m.crd);
      lnids.clear();
      hph.join(TOL, lnids);

      //std::ostringstream o;
      //o << "ph_" << i;
      //medith::write(o.str().c_str(), m.crd, m.cnt, i);

      for (size_t k = 0; k < lnids.size(); ++k)
      {
        if (lnids[k] == E_Int(k)) continue;
        E_Int id1 = hph.poids[k] - 1;
        E_Int id2 = hph.poids[lnids[k]] - 1;
        if (id2 < id1) std::swap(id1, id2);
        nids[id2] = id1;
      }
    }

    m.cnt.PGs.change_indices(nids);

    //for (auto i : modifiedPHs)
    //{
    //  std::ostringstream o;
    //  o << "ph1_" << i;
    //  medith::write(o.str().c_str(), m.crd, m.cnt, i);
    //}
  }
  //std::cout << "nmerges : " << nb_merges << std::endl;
}

///
inline bool sync_faces
(NUGA::ph_mesh_t& m, const std::map<E_Int, std::vector<E_Int>>& face_to_bits, E_Float ARTOL)
{
  if (ARTOL == 0.) ARTOL = -0.01;

  E_Float TOL = ARTOL;
  if (ARTOL < 0.) // relative
  {
    E_Float Lmin, Lmax;
    E_Int imin, imax;
    ngon_type::edge_length_extrema(m.cnt.PGs, m.crd, Lmin, imin, Lmax, imax);
    TOL = -ARTOL * Lmin;
    assert(TOL > 0.); // no degen in m
  }

  // 1. replace faces and keep track of modified PHs
  std::set<E_Int> modifiedPHs;
  replace_faces(m, face_to_bits, TOL, modifiedPHs);

  //3. close_phs
  std::vector<E_Int> modPHs(ALL(modifiedPHs));
  bool has_sync = ngon_type::close_phs(m.cnt, m.crd, &modPHs);

  //for (auto i : modifiedPHs)
  //{
  //  std::ostringstream o;
  //  o << "ph2_" << i;
  //  medith::write(o.str().c_str(), m.crd, m.cnt, i);
  //}

  //4. replace moved master by modified bits : i.e replace original bits by their modified version
  //4.a reverse face_to_bits
  std::map<E_Int, std::vector<E_Int>> face_to_bits_rev, tmp;
  for (auto i : face_to_bits)
  {
    auto bits = i.second;
    E_Int f = i.first;
    for (auto k : bits)
      face_to_bits_rev[k].push_back(f);
  }
  //4.b keep only relevant
  for (auto i : face_to_bits_rev)
  {
    if (i.second.size() == 1) continue;
    tmp[i.first] = i.second;
  }
  face_to_bits_rev = tmp;

  //4.c replace back
  
  if (!face_to_bits_rev.empty())
  {
    std::vector<E_Int> molecPH;
    ngon_unit new_phs;

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
          for (size_t k = 0; k < itBits->second.size(); ++k)
            molecPH.push_back(itBits->second[k] + 1);
        }
      }

      new_phs.add(molecPH.size(), &molecPH[0]);
    }

    new_phs._type = m.cnt.PHs._type;  // hack to preserve flags (externality)
    new_phs._ancEs = m.cnt.PHs._ancEs;// hack

    m.cnt.PHs = new_phs;
    m.cnt.PHs.updateFacets();

    {
      std::vector<E_Int> nids, lnids;
      K_CONNECT::IdTool::init_inc(nids, m.crd.cols());
      for (auto i : modPHs)
      {
        K_MESH::Polyhedron<0> ph(m.cnt, i);
        NUGA::haPolyhedron<0> hph;
        hph.set(ph, m.crd);
        lnids.clear();
        hph.join(TOL, lnids);

        //std::ostringstream o;
        //o << "ph3_" << i;
        //medith::write(o.str().c_str(), m.crd, m.cnt, i);

        for (size_t k = 0; k < lnids.size(); ++k)
        {
          if (lnids[k] == E_Int(k)) continue;
          E_Int id1 = hph.poids[k] - 1;
          E_Int id2 = hph.poids[lnids[k]] - 1;
          if (id2 < id1) std::swap(id1, id2);
          nids[id2] = id1;
        }
      }

      m.cnt.PGs.change_indices(nids);

    }

    has_sync |= ngon_type::close_phs(m.cnt, m.crd, &modPHs);
  }

  //5. clean
  std::vector<E_Int> pgnids, phnids;
  m.cnt.remove_unreferenced_pgs(pgnids, phnids);

  if (pgnids.empty())
    K_CONNECT::IdTool::init_inc(pgnids, m.cnt.PGs.size());

  K_CONNECT::IdTool::reverse_indirection(pgnids, m.cnt.PGs._type); //to pass PG history upon exit

  // assert closed

  return has_sync;
}

#endif
