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
#include <stdexcept>
#include "Nuga/include/Triangulator.h"
#include "Nuga/include/Polygon.h"

#include "face.hpp"
//#define DEBUG_VERBOSE

namespace K_POST
{
    // Calcul d'une tesselation pour la face. On utilise les sommets de la face + le barycentre
    // pour la tesselation. Un indice valant -1 indique qu'on utilise le barycentre sinon c'est
    // un sommet de la face.
    void face::compute_tessellation() const
    {
        E_Int nb_vertices = this->number_of_vertices();
        if (nb_vertices == 3)
        {
            // La face est déjà un triangle. La tessellation est la face elle-même :
            triangles = std::vector<triangle_indices_type>{{this->indices_vertices[0],
                                                            this->indices_vertices[1],
                                                            this->indices_vertices[2]} };
        }
        else if (nb_vertices == 4)
        {
            // Dans le cas d'un quadrangle, on choisit pour ne pas avoir de cas ambigus,
            // de prendre le barycentre et de faire une triangulation étoilée avec
            // le barycentre pour sommet commun à tous les triangles :
            triangles = std::vector<triangle_indices_type>{ {-1, this->indices_vertices[0], this->indices_vertices[1]},
                                                            {-1, this->indices_vertices[1], this->indices_vertices[2]},
                                                            {-1, this->indices_vertices[2], this->indices_vertices[3]},
                                                            {-1, this->indices_vertices[3], this->indices_vertices[0]},
                                                          };
        }
        else {
            // Dans le cas général, on génère une triangulation de Delaunay qui ne prend en compte que
            // les sommets de la face (sans rajouter de nouveaux sommets).
            assert(nb_vertices > 4); // Sinon, ce n'est pas une face ou alors elle est dégénérée de par sa connectivité !
            DELAUNAY::Triangulator dt;
            K_FLD::FloatArray coord(nb_vertices, 3);
            const auto& zone_coords = this->coordinates_of_zone;
            std::vector<E_Int> nodes(nb_vertices);
            std::vector<E_Int> loc2glob(nb_vertices);
            for ( E_Int i_vert = 0; i_vert < nb_vertices; ++i_vert )
            {
                coord(i_vert, 0) = zone_coords[0][indices_vertices[i_vert]];
                coord(i_vert, 1) = zone_coords[1][indices_vertices[i_vert]];
                coord(i_vert, 2) = zone_coords[2][indices_vertices[i_vert]];
                nodes[i_vert] = i_vert;
            }
            K_FLD::IntArray cT3;
            //E_Int err = 
            K_MESH::Polygon::triangulate<DELAUNAY::Triangulator>(dt, coord, nodes.data(), nb_vertices, 0, cT3);
            triangles.reserve(cT3.cols());
            for ( E_Int icol = 0; icol < cT3.cols(); ++icol)
            {
                E_Int* p_ct = cT3.col(icol);
                triangles.emplace_back(std::array<E_Int,3>{p_ct[0],p_ct[1], p_ct[2]});
            }
        }
    }
    // =====================================================================================
    const point3d& face::get_barycenter() const
    {
        // On utilise le pattern proxy. 
        if (std::isnan(barycenter.x))// Est ce que le barycentre a déjà été calculé ?
        {   // Non, donc on le calcule
            this->barycenter = point3d{0.,0.,0.};
            std::size_t nb_vertices = this->indices_vertices.size();
            for ( E_Int ind_vert : this->indices_vertices )
            {
                this->barycenter.x += this->coordinates_of_zone[0][ind_vert];
                this->barycenter.y += this->coordinates_of_zone[1][ind_vert];
                this->barycenter.z += this->coordinates_of_zone[2][ind_vert];
            }
            this->barycenter.x /= nb_vertices;
            this->barycenter.y /= nb_vertices;
            this->barycenter.z /= nb_vertices;
        }
        return this->barycenter;
    }
    // =====================================================================================
    const vector3d& face::get_normal() const
    {
        if (std::isnan(this->normal.x)) // Est ce que la normale a déjà été calculée ?
        {
            // Non, donc on la calcule :
            E_Int nb_vertices = this->indices_vertices.size();
            point3d v0{this->coordinates_of_zone[0][this->indices_vertices[0]],
                       this->coordinates_of_zone[1][this->indices_vertices[0]],
                       this->coordinates_of_zone[2][this->indices_vertices[0]]};
            point3d v1{this->coordinates_of_zone[0][this->indices_vertices[1]],
                       this->coordinates_of_zone[1][this->indices_vertices[1]],
                       this->coordinates_of_zone[2][this->indices_vertices[1]]};
            point3d vn{this->coordinates_of_zone[0][this->indices_vertices[nb_vertices-1]],
                       this->coordinates_of_zone[1][this->indices_vertices[nb_vertices-1]],
                       this->coordinates_of_zone[2][this->indices_vertices[nb_vertices-1]]};
            vector3d v0v1(v0,v1);
            vector3d v0vn(v0,vn);
            this->normal = (v0v1 ^ v0vn);
        }
        return this->normal;
    }    
    // =====================================================================================
    std::pair<bool,bool> 
    face::is_intersecting_ray( const point3d& origin, const vector3d& direction, bool is_with_perturbation ) const
    {
#if defined(DEBUG_VERBOSE)
        std::cout << "rayon : point de départ <= " << std::string(origin) 
                  << ", direction <= " << std::string(direction) << std::endl;
        std::cout << "coordonnées des sommets de la facette : " << std::endl;
        for ( E_Int ind_vert : this->indices_vertices )
            std::cout << std::string(point3d{this->coordinates_of_zone[0][ind_vert],
                                             this->coordinates_of_zone[1][ind_vert],
                                             this->coordinates_of_zone[2][ind_vert]}) << " ";
        std::cout << std::flush << std::endl;
#endif

        auto comp_sign = [&, this] ( int i0, int i1 ) 
        {
            const point3d Pi0 {this->coordinates_of_zone[0][i0], this->coordinates_of_zone[1][i0], this->coordinates_of_zone[2][i0]};
            const point3d Pi1 {this->coordinates_of_zone[0][i1], this->coordinates_of_zone[1][i1], this->coordinates_of_zone[2][i1]};

            vector3d XPi0(origin, Pi0);
            vector3d XPi1(origin, Pi1);
            vector3d n = (XPi0 ^ XPi1);
            n = (1./abs(n))*n;
            E_Float scal = (direction|n);
            return (scal>1.E-14) - (scal<-1.E-14);// En prenant en compte les erreurs numériques...
        };

        auto signe = comp_sign(this->indices_vertices[this->number_of_vertices()-1], this->indices_vertices[0] );
        bool not_change_sign = true;
        // Si on change de signe, c'est que le rayon n'intersecte pas la face
        for (E_Int i = 0; (i < this->number_of_vertices()-1) && (not_change_sign == true); ++i)
        {
            auto signe2 = comp_sign(this->indices_vertices[i],this->indices_vertices[i+1]);
            // Si signe2 est nul, cela signifie que le rayon passe par cette arête.
            // Dans ce cas, on a une indétermination et on prend pas en compte cette arête
            // dans la détermination de l'intersection (on la considère toujours comme valide). 
            if (signe2 != 0)
            {
                if (signe == 0) signe = signe2;
                if (signe2 != signe) not_change_sign = false;
            }
        }
        if (signe == 0) 
        {
            if (!is_with_perturbation)
            {
                // On va faire une légère perturbation pour voir...
                vector3d n = this->get_normal();
                point3d perturb1_origin = origin + (1.E-6)*n;
                decltype(this->is_intersecting_ray(perturb1_origin, direction)) result;
                try
                {
                    result = this->is_intersecting_ray(perturb1_origin, direction, true);
                } catch(std::underflow_error& err)
                {
                    point3d perturb2_origin = origin + (-1.E-6)*n;
                    result = this->is_intersecting_ray(perturb2_origin, direction, true);
                }
                return result;
            }
            else
                throw std::underflow_error("Indetermination case, coplanar face with direction even with perturbation....");
        }
        return {not_change_sign, signe<0};
    }
    // ============================================================================================
    auto face::compute_intersection( const point3d& origin, const vector3d& direction ) const -> intersection_data
    {
        //  I] Si la face n'est pas encore tessallisée, on calcule la tessellation :
        //     ===================================================================
        if (this->triangles.size() == 0)
        {
            this->compute_tessellation();
        }
        assert(this->triangles.size() > 0);

        // II] On va parcourir tous les triangles composant la face :
        //     ====================================================
        const auto& coordinates = this->coordinates_of_zone;
        point3d v1, v2, v3; // Coordonnees des sommets d'un triangle
        point3d intersection;
        E_Int num_trig = 0;
        for ( const auto& trig_ind : triangles )
        {
            // II.1) Constitution du triangle geometrique :
            // ------------------------------------------
            if (trig_ind[0] == -1)
                v1 = this->get_barycenter();
            else 
                v1 = { coordinates[0][trig_ind[0]], coordinates[1][trig_ind[0]], coordinates[2][trig_ind[0]] };            
            if (trig_ind[1] == -1)
                v2 = this->get_barycenter();
            else
                v2 = { coordinates[0][trig_ind[1]], coordinates[1][trig_ind[1]], coordinates[2][trig_ind[1]] };
            if (trig_ind[2] == -1)
                v3 = this->get_barycenter();
            else
                v3 = { coordinates[0][trig_ind[2]], coordinates[1][trig_ind[2]], coordinates[2][trig_ind[2]] };

#if defined(DEBUG_VERBOSE)
            std::cout << "Test intersection avec triangle : " << std::string(v1) << ", "
                      << std::string(v2) << ", " << std::string(v3) << std::flush << std::endl;
#endif
            // ============================================================================
            // II.2) Détection et calcul éventuel de l'intersection du triangle avec le rayon
            // 
            // Algorithme d'intersection de Möller et Trumbore
            // Pour des explication, voir :
            // --------------------------
            // Möller et Trumbore, « Fast, Minimum Storage Ray-Triangle Intersection », 
            // Journal of Graphics Tools, vol. 2,‎ 1997, p. 21–28 
            // ============================================================================
            constexpr const double epsilon = 1.E-12; // Tolérance d'erreur géométrique
            constexpr const double gepsilon = 1.E-3; // Tolérance d'erreur géométrique
            double nrm2 = std::sqrt((direction|direction));
            vector3d ndir = (1./nrm2)*direction;// Normalisation de la direction du rayon
            vector3d edge1(v1,v2);
            vector3d edge2(v1,v3);
            double nrmedge1 = (edge1|edge1);
            double nrmedge2 = (edge2|edge2);
            vector3d h = (ndir ^ edge2);
            double a = (edge1 | h);
#if defined(DEBUG_VERBOSE)
            std::cerr << "ndir : "<<std::string(ndir) << std::endl;
            std::cerr << "edge1 : "<<std::string(edge1) << std::endl;
            std::cerr << "edge2 : "<<std::string(edge2) << std::endl;
            std::cerr << "h : "<<std::string(h) << std::endl;
            std::cerr << "a : " << a << ", nrmedge1 : " << nrmedge1 << ", nrmedge2 : " << nrmedge2 << "epsilon : " << epsilon << std::endl;
#endif
            if (a*a < epsilon*epsilon*nrmedge1*nrmedge2)
            {
                // Pas d'intersection, on continue avec le prochain triangle
                num_trig += 1;
                continue; 
            }
            double f = 1./a;
            vector3d s(v1,origin);
            double u = f * (s|h);
            if ((u<0.) ||(u>1.))
            {
                // Pas d'intersection, on continue avec le prochain triangle
                num_trig += 1;
                continue;
            }
            vector3d q = (s ^ edge1);
            double v = f * (ndir | q);
#if defined(DEBUG_VERBOSE)
            std::cerr << "u = " << u << " et v = " << v << std::flush << std::endl;
#endif
            if ((v>=-gepsilon) && (u+v <= 1.0+gepsilon))
            {
                // Il y a intersection, on reconstitue à partir des coordonnées
                // barycentriques les coordonnées cartésiennes du point d'intersection.
                double t = f * (edge2 | q);
                intersection = origin + t * ndir;
                break; // Intersection trouvée, on quitte la boucle
            }
#if defined(DEBUG_VERBOSE)
            else
            {
                std::cout << "v : "  << v << " et u+v : " << u+v << std::endl;
                double t = f * (edge2 | q);
                intersection = origin + t * ndir;
                std::cout << "intersection non trouvée sinon ce serait : " << std::string(intersection) << std::endl;
            }
#endif
            num_trig += 1;
        }
        if (num_trig >= (E_Int)triangles.size())
        {
            std::cout << "Warning: streamLine2: no intersection found with face." << std::endl;
            return {intersection,-1};
        }
            //throw std::domain_error("Aucune intersection trouvée avec les triangles composant la face");
        return {intersection, num_trig};
    }
    // ============================================================================================
    std::vector<E_Float> 
    face::compute_interpolated_field( const intersection_data& intersection, const FldArrayF& field ) const
    {
        const auto& coordinates = this->coordinates_of_zone;
        const auto& pos_inter   = intersection.first;
        const auto& triangle_indices = triangles[intersection.second];
        point3d v1, v2, v3; // Coordonnees des sommets d'un triangle
        if (triangle_indices[0] >= 0)
        {
            v1 = { coordinates[0][triangle_indices[0]], coordinates[1][triangle_indices[0]],
                   coordinates[2][triangle_indices[0]] };
        }
        else
            v1 = this->get_barycenter();
        if (triangle_indices[1] >= 0)
        {
            v2 = { coordinates[0][triangle_indices[1]], coordinates[1][triangle_indices[1]],
                   coordinates[2][triangle_indices[1]] };
        }
        else
            v2 = this->get_barycenter();
        if (triangle_indices[2] >= 0)
        {
            v3 = { coordinates[0][triangle_indices[2]], coordinates[1][triangle_indices[2]],
                   coordinates[2][triangle_indices[2]] };
        }
        else
            v3 = this->get_barycenter();
#if defined(DEBUG_VERBOSE)
        std::cout << "Coordonnées du triangle considéré : " << std::string(v1) << ","
                  << std::string(v2) << ", " << std::string(v3) << std::flush << std::endl;
#endif
        vector3d edge1(v1,v2);
        vector3d edge2(v1,v3);
#if defined(DEBUG_VERBOSE)
        std::cout << "edge 1 : " << std::string(edge1) << std::endl;
        std::cout << "edge 2 : " << std::string(edge2) << std::endl;
#endif
        int ind1 = 0, ind2 = 1;
        double det = edge1[ind1]*edge2[ind2] - edge1[ind2]*edge2[ind1];
#if defined(DEBUG_VERBOSE)
        std::cout << "det = " << det << std::flush << std::endl;
#endif
        if (std::abs(det) < 1.E-10)
        {
            ind2 = 2;
            det = edge1[ind1]*edge2[ind2] - edge1[ind2]*edge2[ind1];
#if defined(DEBUG_VERBOSE)
            std::cout << "det = " << det << std::flush << std::endl;
#endif
        }
        if (std::abs(det) < 1.E-10)
        {
            ind1 = 1;
            det = edge1[ind1]*edge2[ind2] - edge1[ind2]*edge2[ind1];
#if defined(DEBUG_VERBOSE)
            std::cout << "det = " << det << std::flush << std::endl;
#endif
        }
        assert(std::abs(det) > 1.E-10);
        double inv_det = 1./det;
        vector3d pmb(v1, pos_inter);
        double s = inv_det * (edge2[ind2] * pmb[ind1] - edge2[ind1] * pmb[ind2]);
        double t = inv_det * (edge1[ind1] * pmb[ind2] - edge1[ind2] * pmb[ind1]);
        assert(s+t >= -1.E-6);
        assert(s+t <= 1+1.E-6);
        // On vérifie bien qu'on a bien la coordonnée barycentrique avec la troisième coordonnée
        assert(std::abs(pmb[3-ind2-ind1] - s * edge1[3-ind2-ind1] - t * edge2[3-ind2-ind1]) < 1.E-6);
        std::vector<E_Float> interpol_field(field.getNfld());
        //std::size_t nb_verts = indices_vertices.size();
        std::vector<E_Float> baryfld(field.getNfld(), 0.);
        if ( (triangle_indices[0] == -1) || (triangle_indices[1] == -1) || (triangle_indices[2] == -1) )
        {
            for ( E_Int ind : this->indices_vertices )
                for ( E_Int f = 0; f < field.getNfld(); ++f )
                    baryfld[f] += field(ind, f+1);
            for ( E_Int f = 0; f < field.getNfld(); ++f )
                baryfld[f] /= this->number_of_vertices();
        }
        // on force les 3 premiers champs (supposes coord) a l'intersection de la face
        interpol_field[0] = pos_inter.x; interpol_field[1] = pos_inter.y; interpol_field[2] = pos_inter.z;
        for ( E_Int f = 3; f < field.getNfld(); ++f )
        {
            double v1_val = (triangle_indices[0] >= 0 ? field(triangle_indices[0],f+1) : baryfld[f]);
            double v2_val = (triangle_indices[1] >= 0 ? field(triangle_indices[1],f+1) : baryfld[f]);
            double v3_val = (triangle_indices[2] >= 0 ? field(triangle_indices[2],f+1) : baryfld[f]);
            double delta1 = v2_val - v1_val;
            double delta2 = v3_val - v1_val;
            interpol_field[f] = v1_val + s * delta1 + t * delta2;
        }
        return interpol_field;
    }
}
