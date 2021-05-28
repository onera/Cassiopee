/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_Q6_HXX
#define NUGA_Q6_HXX

#include "Nuga/include/subdivision.hxx"

namespace NUGA
{
  struct Q6
  {
    template <eSUBDIV_TYPE STYPE>
    static void split(ngon_unit& PGs, const E_Int* nodes, NUGA::eDIR d, E_Int firstChild)
    {
      // set them in PGs
      E_Int* q41 = PGs.get_facets_ptr(firstChild);
      E_Int* q42 = PGs.get_facets_ptr(firstChild+1);

      if (d == Xd)
      {
        q41[0] = nodes[0]; q41[1] = nodes[4]; q41[2] = nodes[5]; q41[3] = nodes[3];
        q42[0] = nodes[4]; q42[1] = nodes[1]; q42[2] = nodes[2]; q42[3] = nodes[5];
      }
      else // d == Y
      {
        q41[0] = nodes[0]; q41[1] = nodes[1]; q41[2] = nodes[4]; q41[3] = nodes[5];
        q42[0] = nodes[5]; q42[1] = nodes[4]; q42[2] = nodes[2]; q42[3] = nodes[3];
      }
    }

    ///
    static void retrieve_ordered_data
    (const ngon_unit& PGs, E_Int PGi, E_Int i0, bool reorient, E_Int* two_childrenPG, E_Int* Q6nodes)
    {
      const E_Int* nodes = PGs.get_facets_ptr(PGi);
      E_Int nb_nodes = PGs.stride(PGi);

      for (int i = 0; i < nb_nodes; i++) Q6nodes[i] = nodes[i];

      const E_Int* pNFils0 = PGs.get_facets_ptr(two_childrenPG[0]);
      const E_Int* pNFils1 = PGs.get_facets_ptr(two_childrenPG[1]);

      if (i0 == 0 && !reorient)
      {
        Q6nodes[0] = nodes[0];
        Q6nodes[1] = nodes[1];
        Q6nodes[2] = nodes[2];
        Q6nodes[3] = nodes[3];
        Q6nodes[4] = pNFils0[1];
        Q6nodes[5] = pNFils0[2];
      }
      else if (i0 == 1 && reorient)//a
      {
        Q6nodes[0] = nodes[1];
        Q6nodes[1] = nodes[0];
        Q6nodes[2] = nodes[3];
        Q6nodes[3] = nodes[2];
        Q6nodes[4] = pNFils0[1];
        Q6nodes[5] = pNFils0[2];

        std::swap(two_childrenPG[0], two_childrenPG[1]);
      }
      else if (i0 == 2 && !reorient)//b
      {
        Q6nodes[0] = nodes[2];
        Q6nodes[1] = nodes[3];
        Q6nodes[2] = nodes[0];
        Q6nodes[3] = nodes[1];
        Q6nodes[4] = pNFils0[2];
        Q6nodes[5] = pNFils0[1];

        std::swap(two_childrenPG[0], two_childrenPG[1]);
      }
      else if (i0 == 3 && reorient)//c
      {
        Q6nodes[0] = nodes[3];
        Q6nodes[1] = nodes[2];
        Q6nodes[2] = nodes[1];
        Q6nodes[3] = nodes[0];
        Q6nodes[4] = pNFils0[2];
        Q6nodes[5] = pNFils0[1];
      }
      else if (i0 == 0 && reorient) //d
      {
        Q6nodes[0] = nodes[0];
        Q6nodes[1] = nodes[3];
        Q6nodes[2] = nodes[2];
        Q6nodes[3] = nodes[1];
        Q6nodes[4] = pNFils0[2];
        Q6nodes[5] = pNFils0[3];
      }
      else if (i0 == 1 && !reorient) //e
      {
        Q6nodes[0] = nodes[3];
        Q6nodes[1] = nodes[0];
        Q6nodes[2] = nodes[1];
        Q6nodes[3] = nodes[2];
        Q6nodes[4] = pNFils0[3];
        Q6nodes[5] = pNFils0[2];

        std::swap(two_childrenPG[0], two_childrenPG[1]);
      }
      else if (i0 == 2 && reorient) //f
      {
        Q6nodes[0] = nodes[2];
        Q6nodes[1] = nodes[1];
        Q6nodes[2] = nodes[0];
        Q6nodes[3] = nodes[3];
        Q6nodes[4] = pNFils0[2];
        Q6nodes[5] = pNFils0[3];

        std::swap(two_childrenPG[0], two_childrenPG[1]);
      }
      else if (i0 == 3 && !reorient) //g
      {
        Q6nodes[0] = nodes[1];
        Q6nodes[1] = nodes[2];
        Q6nodes[2] = nodes[3];
        Q6nodes[3] = nodes[0];
        Q6nodes[4] = pNFils0[2];
        Q6nodes[5] = pNFils0[3];
      }
      
    }
      
  };
}

#endif
