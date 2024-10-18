/*    
    Copyright 2013-2024 Onera.

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
//Authors : Sam Landier (sam.landier@onera.fr), Imad Hammani (imadeddine.hammani@onera.fr)

#ifndef NUGA_H12_HXX
#define NUGA_H12_HXX

#include "splitting_t.hxx"
#include "Nuga/include/Hexahedron.h"

//#define DEBUG_HIERARCHICAL_MESH

namespace NUGA
{
	template <>
	class splitting_t<K_MESH::Hexahedron, NUGA::Xd, 1> : public splitting_base_t {
	public:
		E_Int FACES[10];
		E_Int nodes[12];

		template <typename arr_t>
		splitting_t(K_FLD::FloatArray& crd, ngon_type& ng, E_Int PHi, E_Int centroidId,
		const K_FLD::IntArray & F2E, const tree<arr_t> & PGtree) :
			splitting_base_t()
		{
			const E_Int* faces = ng.PHs.get_facets_ptr(PHi);
			//E_Int stride = ng.PHs.stride(PHi);

			E_Int PGi = faces[0]-1;

			if (PGtree.nb_children(PGi) != 2) {
				ok_for_split = false;
				return;
			}

			assert(PGtree.nb_children(PGi) == 2);
			const E_Int *pN = ng.PGs.get_facets_ptr(PGi);

			// get global cut direction
			const E_Int *children = PGtree.children(PGi);
			const E_Int *pChild = ng.PGs.get_facets_ptr(children[0]);
			NUGA::eDIR dir;
			if (pN[1] == pChild[1]) dir = Y;
			else dir = Xd;

			// init
			for (E_Int i = 0; i < 4; i++) nodes[i] = pN[i];

			// local direction depending on orientation of bottom
			bool reorient = need_a_reorient(PGi, PHi, true, F2E);
			if (reorient) {
				std::swap(nodes[1], nodes[3]);
				if (dir == Y) dir = Xd;
				else dir = Y;
			}

			if (dir == Xd)
				cut_in_X(ng, PHi, F2E, PGtree);
			else
				cut_in_Y(ng, PHi, F2E, PGtree);
		}

		template <typename arr_t>
		void cut_in_X(ngon_type& ng, E_Int PHi, const K_FLD::IntArray & F2E, const tree<arr_t> & PGtree)
		{
			E_Int* BOT = FACES;        // DIR
			E_Int* TOP = FACES + 2;    // DIR
			E_Int* LEFT = FACES + 4;   // 0
			E_Int* RIGHT = FACES + 5; // 0
			E_Int* FRONT = FACES + 6; // DIR
			E_Int* BACK = FACES + 8;  // DIR

			E_Int* faces = ng.PHs.get_facets_ptr(PHi);
			E_Int local[6];

			if (PGtree.nb_children(faces[0]-1) != 2 || PGtree.nb_children(faces[1]-1) != 2 ||
				PGtree.nb_children(faces[4]-1) != 2 || PGtree.nb_children(faces[5]-1) != 2) {
					ok_for_split = false;
					return;
				}

			// BOTTOM
			E_Int PGi = faces[0] - 1;
			const E_Int *children = PGtree.children(PGi);
			for (E_Int k = 0; k < 2; k++) BOT[k] = children[k];
			bool reorient = need_a_reorient(PGi, PHi, true, F2E);
			E_Int i0 = 0;
			Q6::retrieve_ordered_data(ng.PGs, PGi, 0, reorient, BOT, local);
			assert(local[0] == nodes[0]);
			assert(local[1] == nodes[1]);
			assert(local[2] == nodes[2]);
			assert(local[3] == nodes[3]);
			nodes[8] = local[4];
			nodes[9] = local[5];

			// FRONT
			PGi = faces[4] - 1;
			const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
			children = PGtree.children(PGi);
			E_Int nb_child = PGtree.nb_children(PGi);
			assert(nb_child == 2);
			for (E_Int k = 0; k < nb_child; k++) FRONT[k] = children[k];
			reorient = need_a_reorient(PGi, PHi, true, F2E);
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
			Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, FRONT, local);
			assert(local[0] == nodes[1]);
			assert(local[1] == nodes[0]);
			assert(local[4] == nodes[8]);
			nodes[4] = local[2];
			nodes[5] = local[3];
			nodes[10] = local[5];

			// BACK
			PGi = faces[5] - 1;
			pN = ng.PGs.get_facets_ptr(PGi);
			children = PGtree.children(PGi);
			nb_child = PGtree.nb_children(PGi);
			assert(nb_child == 2);
			for (E_Int k = 0; k < nb_child; k++) BACK[k] = children[k];
			reorient = need_a_reorient(PGi, PHi, false, F2E);
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
			Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, BACK, local);
			assert(local[0] == nodes[2]);
			assert(local[1] == nodes[3]);
			assert(local[4] == nodes[9]);
			nodes[6] = local[3];
			nodes[7] = local[2];
			nodes[11] = local[5];

			// TOP
			PGi = faces[1] - 1;
			pN = ng.PGs.get_facets_ptr(PGi);
			children = PGtree.children(PGi);
			nb_child = PGtree.nb_children(PGi);
			assert(nb_child == 2);
			for (E_Int k = 0; k < nb_child; k++) TOP[k] = children[k];
			reorient = need_a_reorient(PGi, PHi, false, F2E);
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
			Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, TOP, local);
			assert(nodes[4] == local[0]);
			assert(nodes[5] == local[1]);
			assert(nodes[6] == local[2]);
			assert(nodes[7] == local[3]);
			assert(nodes[10] == local[4]);
			assert(nodes[11] == local[5]);

			// LEFT
			PGi = faces[2] - 1;
			pN = ng.PGs.get_facets_ptr(PGi);
			*LEFT = PGi;
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
			reorient = need_a_reorient(PGi, PHi, true, F2E);
			K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
			assert(local[0] == nodes[0]);
			assert(local[1] == nodes[3]);
			assert(local[2] == nodes[7]);
			assert(local[3] == nodes[4]);

			// RIGHT
			PGi = faces[3] - 1;
			pN = ng.PGs.get_facets_ptr(PGi);
			*RIGHT = PGi;
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
			reorient = need_a_reorient(PGi, PHi, false, F2E);
			K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
			assert(local[0] == nodes[1]);
			assert(local[1] == nodes[2]);
			assert(local[2] == nodes[6]);
			assert(local[3] == nodes[5]);
		}

		template <typename arr_t>
		void cut_in_Y(ngon_type& ng, E_Int PHi, const K_FLD::IntArray & F2E, const tree<arr_t> & PGtree)
		{
			E_Int* BOT = FACES;        // DIR
			E_Int* TOP = FACES + 2;    // DIR
			E_Int* LEFT = FACES + 4;   // DIR
			E_Int* RIGHT = FACES + 6; // DIR
			E_Int* FRONT = FACES + 8; // 0
			E_Int* BACK = FACES + 9;  // 0

			E_Int* faces = ng.PHs.get_facets_ptr(PHi);
			E_Int local[6];

			if (PGtree.nb_children(faces[0]-1) != 2 || PGtree.nb_children(faces[1]-1) != 2 ||
				PGtree.nb_children(faces[2]-1) != 2 || PGtree.nb_children(faces[3]-1) != 2) {
					ok_for_split = false;
					return;
				}

			// BOTTOM
			E_Int PGi = faces[0] - 1;
			const E_Int *children = PGtree.children(PGi);
			for (E_Int k = 0; k < 2; k++) BOT[k] = children[k];
			bool reorient = need_a_reorient(PGi, PHi, true, F2E);
			E_Int i0 = 0;
			Q6::retrieve_ordered_data(ng.PGs, PGi, 0, reorient, BOT, local);
			assert(local[0] == nodes[0]);
			assert(local[1] == nodes[1]);
			assert(local[2] == nodes[2]);
			assert(local[3] == nodes[3]);
			nodes[8] = local[4];
			nodes[9] = local[5];

			// LEFT
			PGi = faces[2] - 1;
			const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
			children = PGtree.children(PGi);
			E_Int nb_child = PGtree.nb_children(PGi);
			assert(nb_child == 2);
			for (E_Int k = 0; k < nb_child; k++) LEFT[k] = children[k];
			reorient = need_a_reorient(PGi, PHi, true, F2E);
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
			Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, LEFT, local);
			assert(local[0] == nodes[0]);
			assert(local[1] == nodes[3]);
			assert(local[4] == nodes[9]);
			nodes[7] = local[2];
			nodes[4] = local[3];
			nodes[11] = local[5];

			// RIGHT
			PGi = faces[3] - 1;
			pN = ng.PGs.get_facets_ptr(PGi);
			children = PGtree.children(PGi);
			nb_child = PGtree.nb_children(PGi);
			assert(nb_child == 2);
			for (E_Int k = 0; k < nb_child; k++) RIGHT[k] = children[k];
			reorient = need_a_reorient(PGi, PHi, false, F2E);
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
			Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, RIGHT, local);
			assert(local[0] == nodes[1]);
			assert(local[1] == nodes[2]);
			assert(local[4] == nodes[8]);
			nodes[6] = local[2];
			nodes[5] = local[3];
			nodes[10] = local[5];

			// TOP
			PGi = faces[1] - 1;
			pN = ng.PGs.get_facets_ptr(PGi);
			children = PGtree.children(PGi);
			nb_child = PGtree.nb_children(PGi);
			assert(nb_child == 2);
			for (E_Int k = 0; k < nb_child; k++) TOP[k] = children[k];
			reorient = need_a_reorient(PGi, PHi, false, F2E);
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
			Q6::retrieve_ordered_data(ng.PGs, PGi, i0, reorient, TOP, local);
			assert(local[0] == nodes[4]);
			assert(local[1] == nodes[5]);
			assert(local[2] == nodes[6]);
			assert(local[3] == nodes[7]);
			assert(local[4] == nodes[10]);
			assert(local[5] == nodes[11]);

			// FRONT
			PGi = faces[4] - 1;
			pN = ng.PGs.get_facets_ptr(PGi);
			*FRONT = PGi;
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
			reorient = need_a_reorient(PGi, PHi, true, F2E);
			K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
			assert(local[0] == nodes[1]);
			assert(local[1] == nodes[0]);
			assert(local[2] == nodes[4]);
			assert(local[3] == nodes[5]);

			// BACK
			PGi = faces[5] - 1;
			pN = ng.PGs.get_facets_ptr(PGi);
			*BACK = PGi;
			i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
			reorient = need_a_reorient(PGi, PHi, false, F2E);
			K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
			assert(local[0] == nodes[2]);
			assert(local[1] == nodes[3]);
			assert(local[2] == nodes[7]);
			assert(local[3] == nodes[6]);
		}

		template <typename pg_arr_t, typename ph_arr_t>
		void split(ngon_type& ng, E_Int PHi, tree<ph_arr_t>& PHtree, tree<pg_arr_t>& PGtree,
			K_FLD::IntArray& F2E, E_Int firstIntPG, E_Int firstPHchild)
		{
			const E_Int* faces = ng.PHs.get_facets_ptr(PHi);
			//E_Int stride = ng.PHs.stride(PHi);

			E_Int PGi = faces[0]-1;
			const E_Int *pN = ng.PGs.get_facets_ptr(PGi);

			// get global cut direction
			const E_Int *children = PGtree.children(PGi);
			const E_Int *pChild = ng.PGs.get_facets_ptr(children[0]);
			NUGA::eDIR dir;
			if (pN[1] == pChild[1]) dir = Y;
			else dir = Xd;

			// local direction depending on orientation of bottom
			bool reorient = need_a_reorient(PGi, PHi, true, F2E);
			if (reorient) {
				if (dir == Y) dir = Xd;
				else dir = Y;
			}

			if (dir == Xd)
				split_in_X(ng, PHi, PHtree, PGtree, F2E, firstIntPG, firstPHchild);
			else
				split_in_Y(ng, PHi, PHtree, PGtree, F2E, firstIntPG, firstPHchild);
		}

		template <typename pg_arr_t, typename ph_arr_t>
		void split_in_X(ngon_type& ng, E_Int PHi, tree<ph_arr_t>& PHtree, tree<pg_arr_t>& PGtree,
			K_FLD::IntArray& F2E, E_Int firstIntPG, E_Int firstPHchild)
		{
			// INTERNAL face creation
			E_Int *q4 = ng.PGs.get_facets_ptr(firstIntPG);

			q4[0] = nodes[8];
			q4[1] = nodes[9];
			q4[2] = nodes[11];
			q4[3] = nodes[10];

			// HX8 creation
			E_Int *BOT(FACES), *TOP(FACES + 2), *LEFT(FACES + 4), *RIGHT(FACES + 5),
				*FRONT(FACES + 6), *BACK(FACES + 8);

			E_Int INT[1] = { firstIntPG };

			E_Int* h121 = ng.PHs.get_facets_ptr(firstPHchild);
			E_Int* h122 = ng.PHs.get_facets_ptr(firstPHchild + 1);

			h121[0] = BOT[0] + 1; h121[1] = TOP[0] + 1;
			h121[2] = LEFT[0] + 1; h121[3] = INT[0] + 1;
			h121[4] = FRONT[1] + 1; h121[5] = BACK[1] + 1;

			assert(K_MESH::Polyhedron<0>::is_closed(ng.PGs, h121, 6));

			h122[0] = BOT[1] + 1; h122[1] = TOP[1] + 1;
			h122[2] = INT[0] + 1; h122[3] = RIGHT[0] + 1;
			h122[4] = FRONT[0] + 1; h122[5] = BACK[0] + 1;

			assert(K_MESH::Polyhedron<0>::is_closed(ng.PGs, h122, 6));

			E_Int PHchildren[2];
			for (size_t k = 0; k < 2; ++k) PHchildren[k] = firstPHchild + k;

			// set them in the tree
			PHtree.set_children(PHi, PHchildren, 2);
			update_F2E(ng, PHi, PHchildren, INT, PGtree, F2E);
		}

		template <typename pg_arr_t, typename ph_arr_t>
		void split_in_Y(ngon_type& ng, E_Int PHi, tree<ph_arr_t>& PHtree, tree<pg_arr_t>& PGtree,
			K_FLD::IntArray& F2E, E_Int firstIntPG, E_Int firstPHchild)
		{
			// INTERNAL face creation
			E_Int *q4 = ng.PGs.get_facets_ptr(firstIntPG);

			q4[0] = nodes[8];
			q4[1] = nodes[9];
			q4[2] = nodes[11];
			q4[3] = nodes[10];

			// HX8 creation
			E_Int *BOT(FACES), *TOP(FACES + 2), *LEFT(FACES + 4), *RIGHT(FACES + 6),
				*FRONT(FACES + 8), *BACK(FACES + 9);

			E_Int INT[1] = { firstIntPG };

			E_Int* h121 = ng.PHs.get_facets_ptr(firstPHchild);
			E_Int* h122 = ng.PHs.get_facets_ptr(firstPHchild + 1);

			h121[0] = BOT[0] + 1; h121[1] = TOP[0] + 1;
			h121[2] = LEFT[0] + 1; h121[3] = RIGHT[0] + 1;
			h121[4] = FRONT[0] + 1; h121[5] = INT[0] + 1;

			assert(K_MESH::Polyhedron<0>::is_closed(ng.PGs, h121, 6));

			h122[0] = BOT[1] + 1; h122[1] = TOP[1] + 1;
			h122[2] = LEFT[1] + 1; h122[3] = RIGHT[1] + 1;
			h122[4] = INT[0] + 1; h122[5] = BACK[0] + 1;

			assert(K_MESH::Polyhedron<0>::is_closed(ng.PGs, h122, 6));

			E_Int PHchildren[2];
			for (size_t k = 0; k < 2; ++k) PHchildren[k] = firstPHchild + k;

			// set them in the tree
			PHtree.set_children(PHi, PHchildren, 2);
			update_F2E(ng, PHi, PHchildren, INT, PGtree, F2E);
		}

		template <typename arr_t>
		void update_F2E(const ngon_type& ng, E_Int PHi, const E_Int *children, const E_Int* INT,
			const tree<arr_t>& PGtree, K_FLD::IntArray & F2E)
		{
			//const E_Int &nbc = subdiv_pol<K_MESH::Hexahedron, DIR>::PHNBC;
			__update_outer_F2E(ng, PHi, children, 2/*nbc*/, PGtree, F2E);

			// INTERNAL face
			F2E(0, INT[0]) = *children;
			F2E(1, INT[0]) = *(children + 1);
		}
	};

	using HX12 = splitting_t<K_MESH::Hexahedron, NUGA::Xd, 1>;
}

#endif
