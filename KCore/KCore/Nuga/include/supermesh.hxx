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
void xmatch(const zmesh_t& m0, const zmesh_t& m1, double RTOL, std::vector<E_Int>& anc0, std::vector<E_Int>& anc1, zmesh_t& xm)
{
  using aelt_t = typename zmesh_t::aelt_t;

  auto loc1 = m1.get_localizer();

  m0.get_nodal_metric2();
  m1.get_nodal_metric2();

  E_Int nbcells0 = m0.ncells();
  std::vector<E_Int> cands;

  std::vector<aelt_t> bits;

#ifdef _OPENMP
  E_Int nb_max_threads = omp_get_max_threads();
#else
  E_Int nb_max_threads = 1;
#endif

  std::vector<zmesh_t> xmi(nb_max_threads);
  std::vector<std::vector<E_Int>> anc0i(nb_max_threads);
  std::vector<std::vector<E_Int>> anc1i(nb_max_threads);
  size_t n, k;
  E_Int i,i2;

#pragma omp parallel private(cands, bits, n, k, i, i2)
  {
#ifdef _OPENMP
    int id = omp_get_thread_num();
#else
    int id = 0;
#endif
    

#pragma omp for //if(do_omp)
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

        for (k = 0; k < bits.size(); ++k)
        {
          xmi[id].add(bits[k]);
          anc0i[id].push_back(i);
          anc1i[id].push_back(i2);
        }
      }
    }
  }

  //std::cout << "mx tre ? : " << nb_max_threads << std::endl;
  for (size_t u = 0; u< nb_max_threads; ++u)
  {
    xm.append(xmi[u]);
    //std::cout << "nb cell for rank : " << u << ": " << xmi[u].ncells() << std::endl;
    anc0.insert(anc0.end(), ALL(anc0i[u]));
    anc1.insert(anc1.end(), ALL(anc1i[u]));
  }
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
