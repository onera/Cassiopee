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
#include <stdexcept>
#include <memory>
#include "structured_data_view.hpp"
#include "unstructured_data_view.hpp"
#include "ngon_data_view.hpp"
#include "Interp/InterpData.h"
#include "Nuga/include/Polygon.h"
#include "stream_line.hpp"
#include "Nuga/include/Triangulator.h"
#include "vector3d.hpp"
#include "point3d.hpp"
using K_POST::vector3d;
using K_POST::point3d;
//#define DEBUG_VERBOSE 

struct K_POST::streamline::Implementation
{
    std::unique_ptr<FldArrayF> pt_field;
    Implementation( const point3d& init_pos, const std::vector<zone_data_view>&  zones, 
                    E_Int max_vertices_for_streamline, bool is_bidirectional );
};

namespace
{
    /**
     * @brief      Retourne le numéro de la zone, le numéro de la cellule et le champs vitesse interpolé
     *             où se situe le point initial init_pos.
     *
     * @param[in]  init_pos                     La position du point initial de la streamline
     * @param[in]  zones                        Les différentes zones
     * @param[in]  max_vertices_for_streamline  Le nombre maximal de sommets pour la streamline
     * @param      streamPt                     Le champs porté par la streamline
     * @param[in]  istream                      La position dans le champs de la valeur initial (utile quand on change de zone)
     * @param[in]  num_zone                     Le numéro de la zone avant changement de zone (-1 : pas de zone courante)
     *
     * @return     Un triplet : numéro de zone, numéro de cellule et vitesse interpolée où se trouve le point init_pos
     */
    std::tuple< E_Int, E_Int, vector3d  >
    init_streamline(const point3d& init_pos, const std::vector<K_POST::zone_data_view>&  zones,
                    E_Int max_vertices_for_streamline, FldArrayF& streamPt, int istream, E_Int num_zone, bool is_with_perturbation=false )
    {
#if defined(DEBUG_VERBOSE)
        std::cout << "Nbre de zones dans lesquels chercher : " << zones.size() << std::endl;        
#endif
        // Recherche du n° de zone et de l'indice élément où se trouve le sommet initial de la streamline :
        decltype(zones.size()) izone = 0;
        E_Int icell = -1;
        for ( const auto& zo : zones )
        {
            if (E_Int(izone) != num_zone)
            {
#if defined(DEBUG_VERBOSE)
                std::cout << "Recherche dans la zone n°" << izone << std::endl;
#endif
                icell = zo.get_interpolation_cell(init_pos);
                if (icell != -1) break;
            }
            ++ izone;
        }
        if (izone == zones.size()) 
        {
            if (!is_with_perturbation)
            {
                // On n'a pas trouvé de domaine adéquat. Le point est sans doute au bord d'un bloc extérieur.
                // On va perturbé le point initial pour trouvé une cellule :
                for ( vector3d perturb : std::array<vector3d,6>{ vector3d{1.E-6,0.,0.}, vector3d{-1.E-6,0.,0.},
                                                                 vector3d{0.,1.E-6,0.}, vector3d{0.,-1.E-6,0.},
                                                                 vector3d{0.,0.,1.E-6}, vector3d{0.,0.,-1.E-6}} )
                {
                    auto pos2 = init_pos + perturb;
                    auto res2 = init_streamline(pos2, zones, max_vertices_for_streamline, streamPt, istream, num_zone, true);
                    if (std::get<0>(res2) != -1) return res2;
                }
            }
            //throw std::domain_error("Wrong starting point : no zones contain this point.");
            return  std::tuple< E_Int, E_Int, vector3d >( -1, -1, {0.,0.,0.});
        }
        //if (izone > zones.size()) throw std::domain_error("Wrong starting point : no zones contain this point.");

#if defined(DEBUG_VERBOSE)
        std::cout << "bloc ou se trouve le point initial = " << izone << std::endl;
        std::cout << "Cellule contenant le point initial : " << icell << std::flush << std::endl;
#endif
        E_Int num_blk = izone;
        const K_POST::zone_data_view& zone = zones[izone];

        // Interpolation du champs de la zone sur le sommet de départ de la streamline :
        E_Int nfld = zone.get_fields()->getNfld();
        if (istream ==0) streamPt.malloc(max_vertices_for_streamline, nfld);
        zone.compute_interpolated_field(init_pos, icell, istream, streamPt);
        auto npos_vel = zone.get_position_of_velocity();
        vector3d velocity_pt{streamPt(istream, npos_vel[0]),streamPt(istream, npos_vel[1]),streamPt(istream, npos_vel[2])};

        // Initialise les coordonnées du premier sommet de la streamline aux coordonnées du point initial :
        auto pos_crds = zones[num_blk].get_position_of_coordinates();
        streamPt(istream,pos_crds[0]) = init_pos.x; streamPt(istream,pos_crds[1]) = init_pos.y; streamPt(istream,pos_crds[2]) = init_pos.z;
        //return {num_blk, icell, velocity_pt};
        return std::tuple< E_Int, E_Int, vector3d >(num_blk, icell, velocity_pt);
    }

    void build_streamline(const point3d& init_pos, 
                          const std::vector<K_POST::zone_data_view>&  zones, 
                          E_Int max_vertices_for_streamline,
                          FldArrayF& streamPt, bool is_reversed = false )
    {
        E_Int ind_cell = -1;
        E_Int num_blk  = -1;
        int   istream;
        vector3d velocity;
        point3d  cur_point = init_pos;
        for (istream = 0; istream < max_vertices_for_streamline; ++istream )
        {
            if (ind_cell == -1)
            {
#               if defined(DEBUG_VERBOSE)
                std::cout << "Recherche d'un nouveau bloc contenant le point " << std::string(cur_point)
                          << " pour le " << istream << " sommet de streamline" << std::endl;
#               endif
                vector3d old_velocity = velocity;

                std::tie(num_blk, ind_cell, velocity) = 
                        init_streamline(cur_point, zones, max_vertices_for_streamline, streamPt, istream, num_blk);
#               if defined(DEBUG_VERBOSE)
                std::cout << "bloc trouvé : " << num_blk << " indice cellule trouvée : " << ind_cell << std::endl;
#               endif
                if (ind_cell == -1)
                {
                    // On perturbe légèrement pour trouver une cellule :
                    point3d newcur_point=cur_point + 1.E-6*old_velocity;
                    std::tie(num_blk, ind_cell, velocity) = 
                        init_streamline(newcur_point, zones, max_vertices_for_streamline, streamPt, istream, num_blk);
#                   if defined(DEBUG_VERBOSE)
                    std::cout << "bloc trouvé : " << num_blk << " indice cellule trouvée : " << ind_cell << std::endl;
#                   endif
                }
                if (ind_cell == -1) break; // On n'a pas trouvé de zones correspondantes...
                continue;
            }
#if defined(DEBUG_VERBOSE)
            E_Int facette_in = -1;
#endif
            E_Int facette_out = -1;
            auto facettes_candidates = zones[num_blk].get_faces_of_element(ind_cell, num_blk);

#if defined(DEBUG_VERBOSE)
            std::cout << "position initiale : " << std::string(cur_point) << std::endl;
            std::cout << "Vitesse initiale  : " << std::string(velocity) << std::endl;
            std::cout << "Facettes candidates : " << std::endl;
            for (const auto& f : facettes_candidates )
            {
                const auto& zone_coords = f.coordinates_of_zone;
                std::cout << "\tindice des coordonnées : ";
                for ( auto fi : f.indices_vertices ) std::cout << fi << " ";
                std::cout << std::flush << std::endl;
                std::cout << "\tCoordonnées des sommets ";
                for ( size_t ivert = 0; ivert < f.indices_vertices.size(); ivert++ ) 
                    std::cout << std::string(point3d{zone_coords[0][f.indices_vertices[ivert]],
                                                     zone_coords[1][f.indices_vertices[ivert]],
                                                     zone_coords[2][f.indices_vertices[ivert]]}) << " ";
                std::cout << std::endl;
            }
#endif
            if ((velocity|velocity) < 1.E-14)
            {
#if defined(DEBUG_VERBOSE)
                std::cout << "Warning: streamLine2: null speed detected. End Streamline computation ..." << std::flush << std::endl;
#endif
                //istream -= 1;
                break;
            }
            //
            if (is_reversed) velocity = -velocity;
            //
            bool is_intersecting, is_entering;
            E_Int ind_facette = 0;
            for ( const auto& facette : facettes_candidates )
            {
                std::tie(is_intersecting, is_entering) = facette.is_intersecting_ray(cur_point, velocity);
                if (is_intersecting)
                {
                    if (is_entering) 
                    { 
#if defined(DEBUG_VERBOSE)
                        facette_in = ind_facette;
#endif
                    }
                    else facette_out= ind_facette;
                }
                ind_facette ++;
            }
#if defined(DEBUG_VERBOSE)
            std::cout << "Facette sortante n°" << facette_out << " trouvée pour intersection rayon : " << std::flush << std::endl;
            auto& out_facette = facettes_candidates[facette_out];
            std::cout << "\tindice des coordonnées : ";
            for ( auto fi : out_facette.indices_vertices ) std::cout << fi << " ";
            std::cout  << std::flush << std::endl;
            std::cout << "Facette rentrante n°" << facette_in << " trouvée pour intersection rayon ! " << std::flush << std::endl;
            auto& in_facette = facettes_candidates[facette_in];
            std::cout << "\tindice des coordonnées : ";
            for ( auto fi : in_facette.indices_vertices ) std::cout << fi << " ";
            std::cout << std::flush << std::endl;
#endif
            // Dans un premier temps, on recherche en aval. On cherchera en amont uniquement si is_bidirectional est vrai
            if (facette_out == -1) break;

            auto& facette = facettes_candidates[facette_out];
#if defined(DEBUG_VERBOSE)
            std::cout << "facette sélectionnée : " << std::endl;
            std::cout << "\tindice des coordonnées : ";
            for ( auto fi : facette.indices_vertices ) std::cout << fi <<  " ";
            std::cout  << std::flush << std::endl;
#endif
            decltype(facette.compute_intersection(cur_point, velocity)) intersect_data;
/*#if !defined(DEBUG_VERBOSE)            
            try
            {
#endif*/
                intersect_data = facette.compute_intersection(cur_point, velocity);
#if defined(DEBUG_VERBOSE)
                //std::cout << ".";
                //std::cout << "Nouvelle position: " << std::string(intersect_data.first) << std::endl;
                //std::cout << "No triangle d'intersection: " << intersect_data.second << "." << std::endl;
#endif
/*#if !defined(DEBUG_VERBOSE)            
            }
            catch(std::domain_error& err)
            {
                std::cerr << "Warning: streamLine2: geometric intersection problem: " << err.what() << ". On arête prématurément le calcul de cette streamline." << std::endl;
                //break;
            }
#endif*/
            if (intersect_data.second >= 0)
            {
            const FldArrayF& field = *zones[num_blk].get_fields();
            std::vector<E_Float> interpfld = facette.compute_interpolated_field(intersect_data, field);
            for ( decltype(interpfld.size()) f = 0; f < interpfld.size(); ++f )
                streamPt(istream, f+1) = interpfld[f];
            auto npos_vel = zones[num_blk].get_position_of_velocity();
            velocity = { streamPt(istream, npos_vel[0]), streamPt(istream, npos_vel[1]), streamPt(istream, npos_vel[2]) };
            ind_cell  = facette.indices_cells.second;
            cur_point = intersect_data.first;
            }
        } // End for (istream)
        E_Int nfld = streamPt.getNfld();
//        std::cout << "Reallocation..." << istream << "," << nfld << std::endl << std::flush ;
        streamPt.reAllocMatSeq(istream,nfld);
//        std::cout << "Fin reallocation..." << std::flush << std::endl;
    }

}

namespace K_POST
{
    streamline::Implementation::Implementation(const point3d& init_pos, const std::vector<zone_data_view>&  zones, 
                                               E_Int max_vertices_for_streamline, bool is_bidirectional)
    {
        this->pt_field = std::unique_ptr<FldArrayF>(new FldArrayF);//std::make_unique<FldArrayF>();
        if (zones.size() == 0) return;
        FldArrayF& streamPt = *this->pt_field;
        build_streamline(init_pos, zones, max_vertices_for_streamline, streamPt, false);
        if (is_bidirectional)
        {
            const K_POST::zone_data_view& zone = zones[0];

            // Interpolation du champs de la zone sur le sommet de depart de la streamline :
            E_Int nfld = zone.get_fields()->getNfld();

            FldArrayF streamPtRev;
            streamPtRev.malloc(max_vertices_for_streamline, nfld);
            if (streamPt.getSize()==0) return;
            for ( E_Int f = 0; f < nfld; ++f )
                streamPtRev(0,f+1) = streamPt(0,f+1);
            build_streamline(init_pos, zones, max_vertices_for_streamline, streamPtRev, true);
            // Fusion des deux streamlines pour en faire une seule
            std::unique_ptr<FldArrayF> pt_globStreamLine= std::unique_ptr<FldArrayF>(new FldArrayF);
            FldArrayF& globStreamLine = *pt_globStreamLine;
            globStreamLine.malloc(streamPt.getSize()+streamPtRev.getSize()-1, nfld);
            for ( E_Int f = 0; f < nfld; ++f )
            {
                for ( E_Int i = 0; i < streamPt.getSize(); ++i )
                {
                    globStreamLine(i + streamPtRev.getSize()-1, f+1) = streamPt(i,f+1);
                }
                for ( E_Int i = 1; i < streamPtRev.getSize(); ++i )
                {
                    globStreamLine(streamPtRev.getSize()-1-i,f+1) = streamPtRev(i, f+1);
                }
            }
            this->pt_field = std::move(pt_globStreamLine);
        }
    }
    // ==============================================================================================================
    streamline::streamline( const point3d& init_pos, const std::vector<zone_data_view>&  zones, 
                            E_Int max_vertices_for_streamline, bool is_bidirectional ) :
        ptr_impl(new Implementation(init_pos, zones, max_vertices_for_streamline, is_bidirectional))
    {
    }// Fin du constructeur
    //_______________________________ Destructeur _____________________________
    streamline::~streamline()
    {
        delete ptr_impl;
    }
    // -------------------------------------------------------------------------------------------------------------
    const FldArrayF& streamline::field() const
    {
        return *this->ptr_impl->pt_field;
    }
    // -------------------------------------------------------------------------------------------------------------
    FldArrayF& streamline::field()
    {
        return *this->ptr_impl->pt_field;
    }
}
