/*
    Copyright 2013-2025 Onera.

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
#include "kdtree.hpp"
#include <algorithm>

namespace K_POST
{
    struct kdtree::node
    {

        node( direction_type dir, E_Int sep_index,
              K_MEMORY::vector_view<E_Int>& indices );

        /**
        * @brief      Test si un noeud de l'arbre contient des sommets
        *
        * @return     Renvoie vrai si le noeud contient des sommets
        */
        bool empty() const { return m_indices.size() == 0; }

        std::pair<E_Int,double> nearest( const point3d& pt, const std::array<K_MEMORY::vector_view<const double>,3>& coords ) const;

        direction_type m_direction; /// Direction de la séparation (en X, en Y, en Z, etc.)
        E_Int m_separation_index;
        K_MEMORY::vector_view<E_Int> m_indices; // Indices des sommets contenus par le noeud
        std::shared_ptr<node> m_left, m_right;
    };
    // ----------------------------------------------------------------------------------
    kdtree::node::node( direction_type dir, E_Int sep_index,
                        K_MEMORY::vector_view<E_Int>& indices) : m_direction(dir), m_separation_index(sep_index),
                                                                 m_indices(indices), 
                                                                 m_left(nullptr), m_right(nullptr)
    {
    }
    // ----------------------------------------------------------------------------------
    std::pair<E_Int,double>
    kdtree::node::nearest(const point3d &pt, const std::array<K_MEMORY::vector_view<const double>, 3> &coords) const
    {
        //std::cout << __PRETTY_FUNCTION__ << std::endl;
        point3d sep_point{ coords[0][m_separation_index], coords[1][m_separation_index], coords[2][m_separation_index]};
        //std::cout << "point de séparation : " << std::string(sep_point) << std::endl;
        double sq_dist_min = square_distance(pt, sep_point);
        //std::cout << "Distance au point de séparation : " << sq_dist_min << std::endl;
        E_Int  min_index = m_separation_index;
        if (sq_dist_min == 0) return {m_separation_index,0.};
        if ((m_left == nullptr) and (m_right == nullptr))
        {
            //std::cout << "Recherche linéaire dans une feuille : " << std::endl;
            // On va rechercher linéairement le noeud le plus proche parmis les noeuds stocké dans ce noeud de l'arbre
            for ( const auto& ind : m_indices )
            {
                point3d vertex{ coords[0][ind], coords[1][ind], coords[2][ind] };
                double sq_dist = square_distance(pt, vertex);
                if (sq_dist < sq_dist_min)
                {
                    sq_dist_min = sq_dist;
                    min_index   = ind;
                }
            }
            //std::cout << "Point trouvé : " << std::string(point3d{coords[0][min_index],coords[1][min_index],coords[2][min_index]}) 
            //          << " à distance² " << sq_dist_min << std::endl;
            return {min_index, sq_dist_min};
        }
        double dx = coords[m_direction][m_separation_index] - pt[m_direction];
        //std::cout << "dx : " << dx << std::endl;
        std::pair<E_Int,double> result_child;
        if (dx>=0)
        {
            if (m_left != nullptr)
            {
                result_child = m_left->nearest(pt, coords);
            }
            else
            {
                result_child = {min_index, sq_dist_min};
            }
        }
        else
        {
            if (m_right != nullptr)
            {
                result_child = m_right->nearest(pt, coords);
            }
            else
            {
                result_child = {min_index, sq_dist_min};                
            }
        }
        if (dx*dx >= result_child.second) return result_child;
        if (result_child.second < sq_dist_min)
        {
            min_index = result_child.first;
            sq_dist_min = result_child.second;
        }
        if (dx >= 0)
        {
            if (m_right != nullptr)
            {
                result_child = m_right->nearest(pt, coords);
            }
            else
            {
                result_child = {min_index, sq_dist_min};                
            }
        }
        else
        {
            if (m_left != nullptr)
            {
                result_child = m_left->nearest(pt, coords);
            }
            else
            {
                result_child = {min_index, sq_dist_min};
            }            
        }
        if ( result_child.second < sq_dist_min)
        {
            //std::cout << "point le plus proche : " << std::string(point3d{coords[0][result_child.first],
            //                                                             coords[1][result_child.first],
            //                                                             coords[2][result_child.first]}) << std::endl;
            //std::cout << "indice: " << result_child.first << " et distance² : " << result_child.second << std::flush << std::endl;
            return result_child;
        }
        //std::cout << "separation_index: " << min_index << " et sq_dist_min : " << sq_dist_min << std::flush << std::endl;
        return {min_index, sq_dist_min};
    }
    // ##################################################################################
    kdtree::kdtree( const K_MEMORY::vector_view<const double>& xCoords, 
                    const K_MEMORY::vector_view<const double>& yCoords,
                    const K_MEMORY::vector_view<const double>& zCoords, unsigned min_nodes ) :
        m_coords{xCoords, yCoords, zCoords},
        m_indices_nodes(xCoords.size())
    {
        assert(xCoords.size() == yCoords.size());
        assert(xCoords.size() == zCoords.size());
        for ( std::vector<E_Int>::size_type i = 0; i < m_indices_nodes.size(); ++i )
            m_indices_nodes[i] = i;
        K_MEMORY::vector_view<E_Int> node_indices(m_indices_nodes.begin(), m_indices_nodes.end());
        m_root = make_tree(0, min_nodes, node_indices);
    }
    // ----------------------------------------------------------------------------------    
    kdtree::kdtree( const std::array<K_MEMORY::vector_view<const double>,3>& coords, unsigned min_nodes ) :
        m_coords(coords),
        m_indices_nodes(coords[0].size())
    {
        for ( std::vector<E_Int>::size_type i = 0; i < m_indices_nodes.size(); ++i )
            m_indices_nodes[i] = i;
        K_MEMORY::vector_view<E_Int> node_indices(m_indices_nodes.begin(), m_indices_nodes.end());
        if (m_indices_nodes.size() < min_nodes)
        {
            m_root = std::make_shared<node>(0, 0, node_indices);
        }
        else m_root = make_tree(0, min_nodes, node_indices);        
    }
    // ----------------------------------------------------------------------------------    
    std::shared_ptr<kdtree::node> 
    kdtree::make_tree(direction_type dir, unsigned min_nodes,
                      K_MEMORY::vector_view<E_Int>& indices)
    {
        if (indices.size() < min_nodes) return nullptr;
        std::size_t n = (indices.size()+1)/2;

        auto cmp_operator = [&,this] ( E_Int ind1, E_Int ind2 ) {
            double val1 = this->m_coords[dir][ind1];
            double val2 = this->m_coords[dir][ind2];
            return val1 < val2;
        };

        std::nth_element(indices.begin(), indices.begin()+n, indices.end(), cmp_operator);
        E_Int separation_index = indices[n];

        std::shared_ptr<node> pt_nd = std::make_shared<node>(dir, separation_index, indices);
        node& nd = *pt_nd;

        dir = (dir+1)%kdtree::dimensions;

        K_MEMORY::vector_view<E_Int> left_indices(indices.begin(), indices.begin()+n);
        nd.m_left = make_tree(dir, min_nodes, left_indices);

        K_MEMORY::vector_view<E_Int> right_indices(indices.begin()+n+1, indices.end());
        nd.m_right = make_tree(dir, min_nodes, right_indices);

        return pt_nd;
    }
    // ----------------------------------------------------------------------------------
    std::pair<E_Int,double> kdtree::nearest(const point3d& pt ) const
    {
        return m_root->nearest(pt, this->m_coords);
    }
}
