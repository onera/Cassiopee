/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_XCELLN_HXX
#define NUGA_XCELLN_HXX

#include "Nuga/include/classifyer.hxx"
#include "Nuga/include/clipper.hxx"

namespace NUGA
{
  ///
  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  class xcellnv : public classifyer<XCELLN_VAL, zmesh_t, bound_mesh_t>
  {
  public:
    using parent_t  = classifyer<XCELLN_VAL, zmesh_t, bound_mesh_t>;
    using wdata_t   = typename parent_t::wdata_t;
    using outdata_t = typename parent_t::outdata_t;
    using elt_t     = typename zmesh_t::elt_t;
    using aelt_t    = typename zmesh_t::aelt_t;

    xcellnv(double RTOL) : parent_t(RTOL){}

    outdata_t __process_X_cells(zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, wdata_t & wdata)
    {
      E_Int ncells = z_mesh.ncells();
      assert(ncells == wdata.size());

      z_mesh.get_nodal_tolerance();

      outdata_t xcelln(ncells, OUT);

      using loc_t = typename zmesh_t::loc_t;

      std::vector<int> cands;
      std::vector<aelt_t> bits, tmpbits; // one clip can produce several bits

      for (E_Int i = 0; i < ncells; ++i)
      {
        color_t const & idata = wdata[i];

        xcelln[i] = double(idata);
        if (xcelln[i] != X) continue;

        bits.clear();
        bits.push_back(z_mesh.aelement(i)); // autonomous element directly stolen by bits (rvalue ref) 

        auto const & bit = bits[0];

        double v0 = bit.extent(); // initial surface/volume
        double vcur(v0), veps(v0*E_EPSILON);
        K_SEARCH::BBox3D bx0;
        bx0.compute(bits[0].m_crd);

        size_t nmask = idata.masks.size();
        for (size_t m = 0; (m < nmask) && (vcur > veps); ++m)
        {
          E_Int mi = idata.masks[m];
          const bound_mesh_t& mask_bit = *mask_bits[mi];
          const loc_t& mask_loc = *(mask_bit.get_localizer());
          
          E_Int idx_start = mask_bit.index_start;
          //
          E_Int nbits = (E_Int)bits.size(); // bits can be appended when split give multiple bits
          for (E_Int b = 0; b < nbits; ++b)
          {
            auto& ae1 = bits[b];

            cands.clear();
            mask_loc.get_candidates(ae1, ae1.m_crd, cands, idx_start);

            bool just_io = cands.empty();
            if (just_io) // try with inital box to catch the bounds and do i/o test
              mask_loc.get_candidates(bx0, cands, idx_start);

            // autonomous cutter front
            bound_mesh_t acut_front(mask_bit, cands);

            //CLIPPING
            tmpbits.clear();
            if (!just_io)
              just_io = !NUGA::CLIP::compute(ae1, acut_front, parent_t::_RTOL, tmpbits); //robust_clip returns true if true clip
 
            // IO : current bit does not intersect front.
            if (just_io)
              if (__classify(ae1, acut_front) == OUT) continue;

            // COMPRESS STRATEGY (with MOVE SEMANTICS) for bits => b can be reduced of 1 to treat the replaced at next iter 
            // if tmpbits is empty (IN) => compress (erase bits if single, put the last instead of current otherwise)
            // else replace the first bit, append the other bits  . 
            __replace_append_and_next_iter(bits, b, tmpbits);
            nbits = bits.size(); //update 
          }

          // accumulated volume
          vcur = 0.;
          for (size_t b = 0; b < nbits; ++b)
            vcur += bits[b].extent();
        }

        xcelln[i] = vcur / v0;
      }

      return std::move(xcelln);
    }

    void __replace_append_and_next_iter(std::vector<aelt_t>& bits, E_Int& b, std::vector<aelt_t>& toadd)
    {
      auto it = toadd.begin();
      if (it != toadd.end())
      {
        bits[b] = std::move(*it);  // replace the first

        bits.reserve(bits.size() + toadd.size() - 1);
        std::move(++it, toadd.end(), std::back_inserter(bits)); // push back remainers
        toadd.clear();                                                          
      }
      else
      {
        if (b < bits.size() - 1) // put the last (if exist) in the current, pop back the last
        {
          bits[b] = std::move(bits.back());
          --b; // to treat it at next iter
        }

        bits.pop_back(); // erase if single, remove the replaced one otherwise
      }
    }

    eClassify __classify(aelt_t const& ae1, bound_mesh_t const& front)
    {
      eClassify res(IN);
      //todo
      return res;
    }
  };
}
#endif