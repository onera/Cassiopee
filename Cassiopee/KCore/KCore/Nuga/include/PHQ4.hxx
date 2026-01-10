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
//Authors : Sam Landier (sam.landier@onera.fr), Jacques Peter (jacques.peter@onera.fr)

#ifndef NUGA_PHQ4_HXX
#define NUGA_PHQ4_HXX

#include "splitting_t.hxx"
#include "Nuga/include/Polyhedron.h"

namespace NUGA
{
  template <>
  class splitting_t<K_MESH::Polyhedron<0>, NUGA::XYZ, 1> : public splitting_base_t
  {

  public:

    E_Int centroidInd1b;

    ///
    template <typename arr_t>
    splitting_t(K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi, E_Int centroidId,
      const K_FLD::IntArray& F2E, const tree<arr_t>& PGtree)
    {
      //calcul et stockage du centroid de PHi dans crd
      E_Int indexstart = 1;
      K_MESH::Polyhedron<UNKNOWN>::iso_barycenter(crd, ng.PGs, ng.PHs.get_facets_ptr(PHi), ng.PHs.stride(PHi), indexstart, crd.col(centroidId));
      centroidInd1b = centroidId + 1;
    }

    ///
    template <typename arr_t>
    void split(ngon_type& ng, E_Int PHi, tree<arr_t>& PHtree, tree<arr_t>& PGtree, K_FLD::IntArray& F2E,
      E_Int firstIntPG, E_Int firstPHchild)
    {
      // 1. dispatcher les quads decoupes (decoupage des faces de PHi a l'etape refine_PGs ) par sommet
      // 2. pour chaque paquet de quad :
      //   a. determiner la ligne polygonale frontiere du paquet
      //   b. "tisser" les quads internes (gestion des doublons ?) avec le centroid et la ligne
      // 3. MAJ de PHtree: PHtree.set_children(PHi, firstPHchild, nbc);
      // 4. MAJ F2E
    
      size_t k, l;
      E_Int j, m, n;

      //== 1 get unodes and number of edges ==============================================================
      //========= build working array for interior faces ====================================  

      std::vector<E_Int> unodes;
      std::vector<E_Int> arity;
      E_Int nb_pgs = ng.PHs.stride(PHi);
      E_Int * first_pg = ng.PHs.get_facets_ptr(PHi);

      K_MESH::Polyhedron<UNKNOWN>::unique_nodes_arity(ng.PGs, first_pg, nb_pgs, unodes, arity);
      E_Int nb_nodesc = unodes.size();

      //== 2 gather all center-edge center-indices about the coarse-level nodes ============================
      //========= build auxiliary data for new (internal) PGi =======================

      std::map <E_Int, std::vector<E_Int>> mapcn2fPG;     // fine level PGs including current coarse node
      std::map <E_Int, std::vector<E_Int>> mapnode;    // corresponding nodes
      E_Int indnode, nb_nodesf, PGi, PGf;
      const E_Int* childPGi;
      E_Int* nodes;
      bool smallfacetoinc;

      for (j = 0; j < nb_nodesc; j++)
      {
        indnode = unodes[j];

        for (k = 0; k < (size_t)nb_pgs; k++)
        {
          PGi = first_pg[k] - 1;
          childPGi = PGtree.children(PGi);

          for (l = 0; l < (size_t)PGtree.nb_children(PGi); l++)
          {
            PGf = childPGi[l];
            nb_nodesf = ng.PGs.stride(PGf);
            nodes = ng.PGs.get_facets_ptr(PGf);
            smallfacetoinc = false;

            for (m = 0; m < nb_nodesf; m++)
            {
              if (nodes[m] == indnode)
              {
                smallfacetoinc = true;
                mapcn2fPG[j].push_back(PGf);
              }
            }
            if (smallfacetoinc == true)
            {
              for (m = 0; m < nb_nodesf; m++)
              {
                if (nodes[m] != indnode)
                  mapnode[j].push_back(nodes[m]);
              }
            }
          }
        }   // loop on children pgs					
      }  // loop on coarse nodes

      //================================================================================================
      //====3 reorganize (the map of) vectors of node indices mapnode ====================
      //======== to get very close to the definition of the fine internal PGs =========================

      std::map <E_Int, std::vector<E_Int>> omapnode;   // reorganized map node
      E_Int threearity, narity;

      for (j = 0; j < nb_nodesc; j++)
      {
        threearity = mapnode[j].size();
        narity = threearity / 3;
        omapnode[j].assign(threearity, 0);

        std::vector<bool> usedinmap(threearity, false);

        // --- check size of  omapnode --------
        omapnode[j][threearity - 1] = mapnode[j][1];
        omapnode[j][0] = mapnode[j][1];
        omapnode[j][1] = mapnode[j][0];
        usedinmap[0] = true;
        usedinmap[1] = true;

        for (n = 0; n < narity - 1; n++)
        {
          for (k = 0; k < (size_t)narity; k++)
          {
            if ((mapnode[j][3 * k] == omapnode[j][1 + 3 * n]) && (usedinmap[3 * k] == false))
            {
              omapnode[j][2 + 3 * n] = mapnode[j][3 * k + 1];
              omapnode[j][3 + 3 * n] = mapnode[j][3 * k + 1];
              omapnode[j][4 + 3 * n] = mapnode[j][3 * k + 2];
              usedinmap[3 * k] = true;
              usedinmap[3 * k + 1] = true;
              usedinmap[3 * k + 2] = true;
            }
            if ((mapnode[j][3 * k + 2] == omapnode[j][1 + 3 * n]) && (usedinmap[3 * k + 2] == false))
            {
              omapnode[j][2 + 3 * n] = mapnode[j][3 * k + 1];
              omapnode[j][3 + 3 * n] = mapnode[j][3 * k + 1];
              omapnode[j][4 + 3 * n] = mapnode[j][3 * k];
              usedinmap[3 * k] = true;
              usedinmap[3 * k + 1] = true;
              usedinmap[3 * k + 2] = true;
            }
          }
        }
      }

      //	
    //====4 loop on coarse-level nodes ===============================================================
    //===== using auxiliary data build  new PGi ======================================================
    // avoid doublon definition of theses PGi by testing second node (middle of edge is in a dedicated list)				
    // no inheritance to be coded check reservation of memory = done in reserve_mem_PHs
    //  mapcn2fPG had only the external faces corresponding to current coarse node > completed with interior faces

      std::vector<E_Int> edgedone;         // to avoid doublon definition of new faces 
      std::map <E_Int, E_Int> midToPGi;       // store indices of middles of edges for which
      E_Int indIntPG = firstIntPG;           // the corresponding new internal face has been created
      bool found = false;

      for (j = 0; j < nb_nodesc; j++)
      {
        indnode = unodes[j];
        E_Int earity = omapnode[j].size() / 3;

        // ==============================================================
        //  (1-based) coarse node number  " << j << " index " << indnode 

        for (k = 0; k < (size_t)earity; k++)
        {
          found = false;
          for (l = 0; l < edgedone.size(); l++)
          {
            if ((omapnode[j])[3 * k + 1] == edgedone[l])
            {
              found = true;
              mapcn2fPG[j].push_back(midToPGi[(omapnode[j])[3 * k + 1]]);
              //std::cout << " just add already created internal face" << midToPGi[(omapnode[j])[3 * k + 1]] << std::endl;
            }
          }
          if (!found)
          {
            E_Int* p = ng.PGs.get_facets_ptr(indIntPG);
            //std::cout << " expected number of nbnodes (hopefully 4) " << ng.PGs.stride(indIntPG) << std::endl;
            mapcn2fPG[j].push_back(indIntPG);
            midToPGi[(omapnode[j])[3 * k + 1]] = indIntPG;
            //std::cout << "   " << (omapnode[j])[3 * k + 1] << " face" << indIntPG << std::endl;
            //   
            indIntPG++;
            p[0] = (omapnode[j])[3 * k];
            p[1] = (omapnode[j])[3 * k + 1];
            p[2] = (omapnode[j])[3 * k + 2];
            p[3] = centroidInd1b;
            //std::cout << " nouveau quadrangle interieur  0-based index " << (indIntPG - 1) << " points 1-based " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
            edgedone.push_back((omapnode[j])[3 * k + 1]);
            //std::cout << " register middle of edge   " << (omapnode[j])[3 * k + 1] << std::endl;
          }
        }
      }
      //std::cout << " ============================================================== " << std::endl;

      //====5 loop on coarse-level nodes ==========================================================
      //======== using auxiliary data build (arity) new PHj ===========================================
      //========== update PHtree ===============================
      //=========  take care of F2E internal faces  =================================

      E_Int indPHChild = firstPHchild;
      STACK_ARRAY(E_Int, nb_nodesc, PHchildren);

      //----------------------------------------------------------------------------------------
      for (j = 0; j < nb_nodesc; j++)        // browse coarse level nodes    
      {
        indnode = unodes[j];
        //std::cout << " ================================================================================= " << std::endl;
        //std::cout << "   create PHi about (1-based) coarse node number " << j << " index " << indnode << std::endl;
        //std::cout << "    arity of the node " << arity[j] << std::endl;
        //for (n = 0; n < mapcn2fPG[j].size(); n++)
          //std::cout << " one 0-based face " << mapcn2fPG[j][n] << std::endl;
        //std::cout << " expected number of faces  " << ng.PHs.stride(indPHChild) << std::endl;

        //------------------------------------------------

        E_Int* p = ng.PHs.get_facets_ptr(indPHChild);

        //std::cout << " pour indPHChild " << indPHChild << " donnees en " << p << std::endl;

        //------------ external faces -----------------------
        // -------- the division kept the orientation -----------------
        //
        for (k = 0; k < (size_t)arity[j]; k++)       // the new PHj has arity[j] old "external" faces (all quadrangles)
                                            // and arity[j] "internal" faces (all quadrangles) 
        {
          p[k] = mapcn2fPG[j][k] + 1;         // mapcn2fPG is 0-based
          if (F2E(0, mapcn2fPG[j][k]) == PHi)             // usefulness of the initialisation
            F2E(0, mapcn2fPG[j][k]) = indPHChild;

          if (F2E(1, mapcn2fPG[j][k]) == PHi)             // usefulness of the initialisation
            F2E(1, mapcn2fPG[j][k]) = indPHChild;
        }
        //------------ internal faces -----------------------
        //--------------------------------------------------
        for (k = arity[j]; k < (size_t)(2 * arity[j]); k++)   // the new PHj has arity[j] old "external" faces (all quadrangles)
                                                 // and arity[j] "internal" faces (all quadrangles) 
        {
          p[k] = mapcn2fPG[j][k] + 1;              // mapcn2fPG est 0-based
          //std::cout << " affectation un indice face =  " << mapcn2fPG[j][k] + 1 << std::endl;

          if (F2E(0, mapcn2fPG[j][k]) == E_IDX_NONE)
            F2E(0, mapcn2fPG[j][k]) = indPHChild;          //0-based 0-based
          else
            F2E(1, mapcn2fPG[j][k]) = indPHChild;          //0-based 0-based
        }
        //--------------------------------------------------
        PHchildren[j] = indPHChild;
        indPHChild++;
      }//coarse nodes

      PHtree.set_children(PHi, PHchildren.get(), nb_nodesc);

      //== 6=== final F2E for internal faces =====================================================
      //========== go for doublon loop on coarse nodes ==========================================

      //E_Int indComFace;
      //E_Bool comInternalFace;
      //
      // Note (Imad): j'ai commenté cette for loop qui n'accomplit rien a priori
      /*
      for (j1 = 0; j1 < nb_nodesc; j1++)        // browse coarse level nodes    
      {
        E_Int indnode1 = unodes[j1];
        for (j2 = 0; j2 < j1; j2++)        // browse coarse level nodes    
        {
          E_Int indnode2 = unodes[j2];
          //std::cout << " couple de noeuds envisage " << indnode1 << " -- " << indnode2 << std::endl;
          //---------------------------------------------------
          //comInternalFace = false;
          // CBX: louche
          for (k1 = 0; k1 < 2 * arity[j1]; k1++)
            for (k2 = 0; k2 < 2 * arity[j2]; k2++)
              if (mapcn2fPG[j1][k1] == mapcn2fPG[j2][k2])
              {
                //comInternalFace = true;
                //indComFace = mapcn2fPG[j2][k2];
                break;
              }
          //-------------------------------------------------
          //if (comInternalFace)
            //std::cout << " face interne et F2E a examiner indice face " << indComFace << std::endl;
        }
      }
      */
    }
  };

  using PHQ4 = splitting_t<K_MESH::Polyhedron<0>, NUGA::XYZ, 1>;
}

#endif
