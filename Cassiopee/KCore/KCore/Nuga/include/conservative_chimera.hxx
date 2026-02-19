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

//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef NUGA_CONSERVATIVE_CHIMERA_HXX
#define NUGA_CONSERVATIVE_CHIMERA_HXX

#include "Nuga/include/NGON_BooleanOperator.h"
#include "Nuga/include/Polygon.h"

namespace NUGA
{
  static const E_Int INTERP = 2;
  static const E_Float RESET_VAL = NUGA::FLOAT_MAX;
  static const E_Float COEF_THRESHOLD2 = EPSILON*EPSILON;
  
  template <typename crd_t, typename cnt_t>
  class P1_Conservative_Chimera
  {
    public:
      typedef ngon_t<cnt_t> ngon_type;
      typedef K_FLD::ArrayAccessor<crd_t> acrd_t;
    
    public:
        
  /**
   * Computes conservative chimera donnor coefficients and indices for receiver cells to be interpolated (with cellN == 2)
   * acrdR : accessor for coordinates of the donnor mesh
   * cntR : receiver connectivity
   * acrdD : accessor for coordinates of the receiver mesh
   * cntD : donnor connectivity
   * cellNR : celn nature field for the receiver mesh
   * dindices : donnor indices per receiver cell (ranged with xr)
   * dcoeffs : donnor coefficients per receiver cell (ranged with xr)
   * xr : range delimiter in dindices and dcoefs : the i-th receiver cell has info between xr[i] and xr[i+1]
   */
  
  template <typename TriangulatorType> 
  static E_Int compute_coeffs
  (const crd_t& crdR, const cnt_t& cntR,
   const crd_t& crdD, const cnt_t& cntD,
   const Vector_t<E_Int>& cellNR,
   Vector_t<E_Int>& dindices, Vector_t<E_Float>& dcoeffs,
   Vector_t<E_Int>& xr, Vector_t<E_Int>& roids)
  {
    E_Int err = 0;

    dindices.clear();
    dcoeffs.clear();
    xr.clear();
    xr.push_back(0);
    roids.clear();
    
#ifdef DEBUG_CONS_CHIMERA
    std::cout << "compute_coeffs inputs : " << std::endl;
    std::cout << " -- nb points donnor : " << crdD.getSize() << std::endl;
    std::cout << " -- nb points receiver : " << crdR.getSize() << std::endl;
    std::cout << " -- nb cells donnor : " << cntD[2 + cntD[1]] << std::endl;
    std::cout << " -- nb cells receiver : " << cntR[2 + cntR[1]] << std::endl;
    std::cout << " -- cellN size : " << cellNR.size() << std::endl;
    std::cout << std::endl;
#endif
    
    // Create the operator
    typedef NUGA::NGON_BooleanOperator<crd_t, cnt_t> boolean_t;
    boolean_t oper (crdR, cntR, crdD, cntD, 0., boolean_t::FULL);

    // I. DO THE INTERSECTION
    K_FLD::FloatArray crdS; //named as (crdS, cntS) because the operation is soft, the result is directly the supermes bits
    K_FLD::IntArray cntS;  
    
#ifdef DEBUG_CONS_CHIMERA
    std::cout << "compute_coeffs : step 1 : boolean intersection" << std::endl;
#endif

    err = oper.Intersection(crdS, cntS, boolean_t::SOLID_NONE, boolean_t::PRESERVE_RIGHT);
    if (!oper._ngoper)
    {
#ifdef DEBUG_CONS_CHIMERA
    std::cout << "compute_coeffs : ERROR step 1" << std::endl;
#endif
      return 1;
    }
    
    const ngon_type& ngS = *oper._ngoper; //supermesh
    ngon_type ngD(cntD), ngR(cntR);
    
#ifdef DEBUG_CONS_CHIMERA
    
    K_FLD::IntArray cnto;
    ngS.export_to_array(cnto);
    //medith::write("supermesh.plt", crdS, cnto, "NGON");
    
    std::cout << "compute_coeffs : step 2 : reorient 3D" << std::endl;
#endif
    
    
    
    // II. build neighborhood and orientation on donnor mesh (required to compute the gradient)
    // WARNING : MUST BE DONE BEFORE COMPUTING THE METRICS BECAUSE REORIENTING CHANGES ngD
    
    ngon_unit orientD;
    err = ngon_type::template build_orientation_ngu<TriangulatorType>(crdD, ngD, orientD);
    
    ngon_unit neighborsD;
    ngD.build_ph_neighborhood(neighborsD);
    
    if (err)
    {
#ifdef DEBUG_CONS_CHIMERA
      std::cout << "compute_coeffs : ERROR step 2" << std::endl;
#endif
      return err;
    }
      
    
#ifdef DEBUG_CONS_CHIMERA
    std::cout << "compute_coeffs : step 3 : compute_metrics" << std::endl;
#endif
    
    // III. COMPUTE THE METRICS 
    // volumes
    Vector_t<E_Float> volumeD, volumeR, volumeS; // donnor, receiver, supermesh
    // donnor surface vectors & centroids
    K_FLD::FloatArray surfaceD, centroidD, centroidS;

    ///
    err = compute_metrics<TriangulatorType>(crdD, ngD, neighborsD, crdR, ngR, crdS, ngS, cellNR, 
                                      volumeD, volumeR, volumeS,
                                      centroidD, centroidS, surfaceD);
    
    
   if (err)
    {
#ifdef DEBUG_CONS_CHIMERA
      std::cout << "compute_coeffs : ERROR step 3" << std::endl;
#endif
      return err;
    }
    
#ifdef DEBUG_CONS_CHIMERA
    std::cout << "compute_coeffs : step 4 : receiver partition map" << std::endl;
#endif
    
    // IV. RECEIVER PARTITION MAP
    // count nb of bits covering each receiver cell
    std::map<E_Int, Vector_t<E_Int> > receiver_to_bits;
    E_Int nb_phs = ngS.PHs.size();
    for (E_Int i=0; i < nb_phs; ++i)
    {
      const E_Int& ancR = ngS.PHs._ancEs(0, i);
      
#ifdef DEBUG_CONS_CHIMERA
      if (ancR < 0 || ancR >= cellNR.size())
      {
        std::cout << "ERROR : missing ancestor" << std::endl;
        return 1;
      }
#endif
      
      if (cellNR[ancR] != INTERP) continue;
      
      receiver_to_bits[ancR].push_back(i);
    }
    
#ifdef DEBUG_CONS_CHIMERA
    std::cout << "nb of receivers that have donnors/coeffs  : " << receiver_to_bits.size() << std::endl;
    std::cout << "compute_coeffs : step 5 : coeffs computation" << std::endl;
    //NGON_debug<crd_t, cnt_t>::draw_PHT3(crdD, ngD, 1121);
#endif
    
    // V. COMPUTE THE COEFS
    std::map<E_Int, Vector_t<E_Int> >::const_iterator itR, itb(receiver_to_bits.begin()), ite(receiver_to_bits.end());
    std::map<E_Int, E_Float> molecule;

    E_Float GDkgk[3]; // vector from donnor centroid to bit centroid
    
    //Vector_t<E_Int> th_dinds, th_xr,th_oids; FOR OPENMP
    //Vector_t<E_Float> th_dcoefs;
   
// fixme : parallel not enabled yet because of the variable-stride appending upon exit.. 
//#pragma omp parallel for private (molecule, GDkgk, th_dinds, th_xr,th_oids, th_dcoefs)
    for (itR=itb; itR != ite; ++itR)
    {
      const E_Int& ancR = itR->first;
      molecule.clear();
      
#ifdef DEBUG_CONS_CHIMERA
      //std::cout << "RECEIVER : " << ancR << std::endl;
#endif
      
      const Vector_t<E_Int>& bits = itR->second;
      
      E_Float VtildR=0.; // Covered volume
      E_Int nb_bits(bits.size());
      for (E_Int k=0; k < nb_bits; ++k)
      {
        E_Int bit = bits[k];
        E_Int Dk = ngS.PHs._ancEs(1, bit); //direct current donnor
        
        E_Float vk = volumeS[bit];
        VtildR += vk;

        // 1. Coefficient for direct donnor
        if (molecule.find(Dk) == molecule.end())
          molecule[Dk] = vk;
        else
          molecule[Dk] +=vk;

        // 2. Coefficient for cells adjacent to direct donnor   
        // (reflecting that we take into account for gradient)
        const E_Float* GDk = centroidD.col(Dk);  // centroid of the direct donnor
        const E_Float* gk = centroidS.col(bit); // centroid of the donnor cell's current piece
        NUGA::diff<3>(gk, GDk, GDkgk);

        E_Float coef0 = 0.5 * vk;
        if (volumeD[Dk] != 0.) coef0 /= volumeD[Dk];

#ifdef DEBUG_CONS_CHIMERA
        E_Float cs[] = {0.,0.,0.};
#endif

        const E_Int* pQ  = ngD.PHs.get_facets_ptr(Dk);     // Dk quads
        const E_Int* pDkj = neighborsD.get_facets_ptr(Dk); // Dk neighbors
        const E_Int* pOr = orientD.get_facets_ptr(Dk);     // Dk faces orientations

        E_Int nb_neighs = neighborsD.stride(Dk);
        for (E_Int j=0; j < nb_neighs; ++j)
        {
          E_Int Dkj = *(pDkj + j);

          if (Dkj == IDX_NONE) //Extrapolation policy => modify Dk coeffs
            Dkj = Dk;

          E_Int Qj = *(pQ + j) - 1;  // Quad beteween Dk and Dkj
          E_Int ori = (*(pOr + j) == 1) ? 1. : -1.;

#ifdef DEBUG_CONS_CHIMERA
//          std::cout << "NODES : " ;
//          for (E_Int n=0; n < ngD.PGs.stride(Qj); ++n) std::cout << ngD.PGs.get_facet(Qj, n) << " ";
//          std::cout << std::endl;
//          std::cout << "ori : " << ori << std::endl;          
//          std::cout << "oriented surface : " <<  ori*surfaceD(0,Qj) << " " << ori*surfaceD(1,Qj) << " " << ori*surfaceD(2,Qj) << std::endl;
//          //std::cout << "GDkgk : " <<  GDkgk[0] << " " << GDkgk[1] << " " << GDkgk[2] << std::endl;
//          cs[0] += ori * surfaceD(0,Qj);
//          cs[1] += ori * surfaceD(1,Qj);
//          cs[2] += ori * surfaceD(2,Qj);
#endif

          E_Float coeff = ori *coef0 * NUGA::dot<3>(GDkgk, surfaceD.col(Qj));

          if (molecule.find(Dkj) == molecule.end())
            molecule[Dkj] = coeff;
          else
            molecule[Dkj] +=coeff;
        }

#ifdef DEBUG_CONS_CHIMERA
//        std::cout << cs[0] << " " << cs[1] << " " << cs[2] << std::endl;
//        assert ( (fabs(cs[0]) < EPSILON) && (fabs(cs[1]) < EPSILON) && (fabs(cs[2]) < EPSILON) );
#endif

      }

      // final division by VtildR
      VtildR = 1. / VtildR;
      std::map<E_Int, E_Float>::iterator itC = molecule.begin();
      for (; itC != molecule.end(); ++itC)
      {
        itC->second *= VtildR;
        //std::cout << "DONNOR : " << itC->first << " . coefs : " << itC->second << std::endl;
      }

#ifdef DEBUG_CONS_CHIMERA
//      std::cout << std::endl;
      E_Float coefs_sum = 0.;
#endif

      // info is ready for passing upon exit
      roids.push_back(ancR);
      
      E_Int nb_coefs=0;
      
      for (itC = molecule.begin(); itC != molecule.end(); ++itC)
      {
        
        const E_Float& val = itC->second;
        
        if (val*val < COEF_THRESHOLD2) continue; // UNSIGNIFICANT COEFF

#ifdef DEBUG_CONS_CHIMERA
        if (val > 10./*1. + EPSILON*/)
        {
          std::cout << "ERROR : coef out of range" << std::endl;
          std::cout << "receiver : " << ancR << std::endl;
          std::cout << "donnor " << itC->first << std::endl;
          std::cout << "coefficient " << itC->second << std::endl;
          err = 1;
          return 1;
        }
#endif
        
        dindices.push_back(itC->first);
        dcoeffs.push_back(itC->second);
        
        ++nb_coefs;

#ifdef DEBUG_CONS_CHIMERA
        coefs_sum += val;
#endif
        
      }

      xr.push_back(nb_coefs + xr[xr.size()-1]);
      
#ifdef DEBUG_CONS_CHIMERA
//      E_Float delta = fabs(coefs_sum-1.);
//      E_Int last = xr.size()-1;
//      E_Int nb_donnor =xr[last] - xr[last-1];
//      if (delta >= EPSILON)
//      {
//        std::cout << "partial corvering for cell : " << roids.back() << "." << std::endl;
//        for (E_Int i=xr[last-1]; i < xr[last]; ++i) std::cout << dcoeffs[i] << " ";
//        std::cout << std::endl;
//        for (E_Int i=0; i < nb_donnor; ++i) std::cout << dindices[i] << " ";
//        std::cout << std::endl;
//      }
      //assert (fabs(coefs_sum-1.) < EPSILON) ;
#endif  
      
    }

    return err;
  }
  
  ///
  template <typename TriangulatorType> 
  static E_Int compute_metrics
  (const crd_t& crdD, const ngon_type& ngD, const ngon_unit& neighborsD,
   const crd_t& crdR, const ngon_type& ngR, 
   const crd_t& crdS, const ngon_type& ngS, 
   const Vector_t<E_Int>& cellNR, 
   Vector_t<E_Float>& volumeD,
   Vector_t<E_Float>& volumeR,
   Vector_t<E_Float>& volumeS,
   K_FLD::FloatArray& centroidD,
   K_FLD::FloatArray& centroidS,
   K_FLD::FloatArray& surfaceD)
  {
    TriangulatorType dt;
    E_Int nb_phs = ngS.PHs.size();
    
    volumeD.clear();
    volumeD.resize(ngD.PHs.size(), RESET_VAL);     // donnor
    volumeR.clear();
    volumeR.resize(ngR.PHs.size(), RESET_VAL);     // receiver
    volumeS.clear();
    volumeS.resize(ngS.PHs.size(), RESET_VAL);     // supermesh
    surfaceD.clear();
    surfaceD.resize(3, ngD.PGs.size(), RESET_VAL); // donnor
    centroidD.clear();
    centroidD.resize(3, ngD.PHs.size(), RESET_VAL);// donnor
    centroidS.clear();
    centroidS.resize(3, ngS.PHs.size(), RESET_VAL);// supermesh
    
    E_Int err=0;
    
//#pragma omp parallel for shared(err)/*private (XXX)*/
    for (E_Int i = 0; i < nb_phs ; ++i)
    {
      const E_Int& ancR = ngS.PHs._ancEs(0, i);
      const E_Int& ancD = ngS.PHs._ancEs(1, i);
      
      if (err) continue;
      
      if (ancD  == IDX_NONE)
      {
        std::cout << "faulty donnor ancestor for piece number " << i << " in the supermesh : " <<  ancR << "/" <<  ancD << std::endl;
        err = 1;
        continue;
      }
      if (ancR  == IDX_NONE)
      {
        std::cout << "faulty receiver ancestor for piece number " << i << " in the supermesh : " <<  ancR << "/" <<  ancD << std::endl;
        err = 1;
        continue;
      }
      
      if (cellNR[ancR] != INTERP) continue;
      
      // SUPERMESH CENTROID AND VOLUME : assumed to be centroid-star-shaped as well (because operation is intersection)
      const E_Int* first_pg = ngS.PHs.get_facets_ptr(i);
      E_Int nb_pgs = ngS.PHs.stride(i);
      K_MESH::Polyhedron<STAR_SHAPED>::metrics(dt, crdS, ngS.PGs, first_pg, nb_pgs, volumeS[i], centroidS.col(i));
      
      // DONNOR METRICS : CENTROID, VOLUME AND SURFACE : assumed to be centroid-star-shaped
      // AND NEIGHBORS (for gradients)
      const E_Int* pDkj = neighborsD.get_facets_ptr(ancD);
      E_Int nb_neighs = neighborsD.stride(ancD);
      // so gather in one loop direct donnor (first item) and neighbors
      E_Int DONcur;//current donnor
      for (E_Int j=-1; j<nb_neighs; ++j)
      {
        DONcur = (j==-1) ? ancD : *(pDkj + j); //direct or neighnbor ?
        if (DONcur == IDX_NONE) continue;

        // Loop over current donnor faces
        const E_Int* pF = ngD.PHs.get_facets_ptr(DONcur);
        E_Int nbF = ngD.PHs.stride(DONcur);
        
        for (E_Int k=0; k < nbF; ++k)
        {
          E_Int PGk = *(pF + k) - 1;
          E_Float * pS = surfaceD.col(PGk);
          
          if (*pS != RESET_VAL) continue; //already computed
          
          const E_Int* nodes = ngD.PGs.get_facets_ptr(PGk);
          E_Int nb_nodes = ngD.PGs.stride(PGk);
          K_MESH::Polygon::ndS<crd_t, 3>(crdD, nodes, nb_nodes, 1, pS);
          
        }
        
        if (volumeD[DONcur] != RESET_VAL) continue; //already computed
        // DONNOR VOLUME : Structured inputs => centroid-star-shaped
        K_MESH::Polyhedron<STAR_SHAPED>::metrics(dt, crdD, ngD.PGs, pF, nbF, volumeD[DONcur], centroidD.col(DONcur));
      }
      
      // RECEIVER VOLUME : assumed to be centroid-star-shaped
      if (volumeR[ancR] == RESET_VAL) // not already computed
      {
        // Structured inputs => centroid-star-shaped
        const E_Int* first_pg = ngR.PHs.get_facets_ptr(ancR);
        E_Int nb_pgs = ngR.PHs.stride(ancR);
        E_Float centroid[3];//unused
        K_MESH::Polyhedron<STAR_SHAPED>::metrics(dt, crdR, ngR.PGs, first_pg, nb_pgs, volumeR[ancR], centroid);
      }
    }
    //std::cout << "return value is :" << err << std::endl;
    return err;
  }
  
// fixme : TO VALIDATE WHEN FIXING THE CURVILINEAR ISSUE
  static E_Int chim_conservative_transfer
(const crd_t& dcrd, const K_FLD::IntArray& dcnt, const crd_t& rcrd, const K_FLD::IntArray& rcnt, 
 const K_FLD::FloatArray& fields_don, K_FLD::FloatArray& fields_rec, E_Int field = -1/*-1 means all*/)
{
  //
    
  DELAUNAY::Triangulator dt;

  std::vector<E_Int> rids/*recpetor ids*/, xdon /*delimiter indirection for dids and vcoefs*/, dids, pids;
  std::vector<E_Float> dcoefs;
  K_FLD::FloatArray piece_centroids, pcrd;
  
  ngon_type ngr = rcnt;
  ngon_type ngd = dcnt;
  std::vector<E_Int> cellNR(ngr.PHs.size(), INTERP);
  
  E_Int err = NUGA::P1_Conservative_Chimera<K_FLD::FloatArray, K_FLD::IntArray>::compute_coeffs<DELAUNAY::Triangulator>(rcrd,rcnt , dcrd, dcnt, cellNR, dids, dcoefs, xdon, rids);
  if (err==1)
    return 1;
  
  K_FLD::FloatArray don_centroids;
  //gcc error: expected primary-expression before '>' token
  ngon_type::template centroids<DELAUNAY::Triangulator>(ngd, dcrd, don_centroids);

  E_Int nb_sols = 1;
  E_Int nb_fields = fields_don.rows();
  // if out of sacope value, consider all, only one otherwise
  if (field <= -1 || field >= nb_fields)
  {
    field=0;
    nb_sols = nb_fields;
  }
  
  E_Int nb_rphs = ngr.PHs.size();

  //Redimension in case it is required
  fields_rec.resize(nb_sols, nb_rphs, 0.);

#ifdef FLAG_STEP
  std::cout << "transfer " << nb_sols << " PWL fields on " << nb_rphs << " elements..." << std::endl;
#endif
  
  //E_Float GcGcp[3];
  
  // Now tranfer
  for (E_Int i=0; i < rids.size(); ++i)
  {
    E_Int& reci = rids[i];
    if (reci == 0)
    {
      std::cout << "dealing with " << i << "-th receptor : " << reci << std::endl;
      std::cout << "donnor fieces from " <<  xdon[i] << " to " << xdon[i+1] << std::endl;
    }
    for (E_Int f=0; f < nb_sols; ++f)
    {
      E_Int fld = (field+f)%nb_fields;
      fields_rec(fld, reci) = 0.;
      
      for (E_Int k=xdon[i]; k<xdon[i+1]; ++k)
        fields_rec(fld, reci) += fields_don(fld, dids[k]) * dcoefs[k];
      


    }
  }
  
  return 0;
}
 
};

}

#endif
