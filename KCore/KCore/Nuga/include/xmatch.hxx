/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_XMATCH_HXX
#define NUGA_XMATCH_HXX

//#define XMATCH_DBG

#include "Nuga/include/clipper.hxx"


#ifdef XMATCH_DBG
#include "Nuga/include/medit.hxx"
#endif


namespace NUGA
{

template <typename zmesh_t>
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

}

#endif // NUGA_XMATCH_HXX
