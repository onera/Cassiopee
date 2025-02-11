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
#include <cassert>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <list>
#include "triangulated_polyhedron.hpp"
#include "Nuga/include/Polygon.h"
#include "Nuga/include/Triangulator.h"
#include "ngon_data_view_p.hpp"
#include "vector3d.hpp"
#include "linear_algebra.hpp"
#include "rotational.hpp"
#include "volume.hpp"

//#define DEBUG_VERBOSE
namespace K_POST
{
    //# #####################################################################################################
    //_ _                                    Mise en oeuvre de la classe privée                             _
    //# #####################################################################################################

    //_ ______________________________________ Constructeur _________________________________________________
    ngon_data_view::Implementation::Implementation( E_Int nb_faces, const const_view_type<E_Int>& face2verts, 
                                                    E_Int nb_elts , const const_view_type<E_Int>& elt2Faces,
                                                    const fields_type& fields, const coordinates_npos& pos_coords, 
                                                    const coordinates_npos& pos_velocity) : 
        zone_data_view::Implementation( zone_data_view::NGON, pos_coords, fields, pos_velocity)
    {
        // On reconstruit la connectivité à partir de la connectivité CGNS qui n'est pas top top.
        // On construit le tableau d'indirection pour les faces vers les sommets :
        std::vector<E_Int>(nb_faces+1).swap(m_beg_face2verts);
        m_beg_face2verts[0] = 0;
        E_Int index = 0;
        for ( E_Int iface = 0; iface < nb_faces; ++iface)
        {
            E_Int nb_verts = face2verts[index];
            index += nb_verts+1;
            m_beg_face2verts[iface+1] = m_beg_face2verts[iface] + nb_verts;
        }
        // Puis la connectivité face vers les sommets :
        m_face2verts.reserve(m_beg_face2verts[nb_faces]);
        index = 0;
        for ( E_Int iface = 0; iface < nb_faces; ++iface)
        {
            E_Int nb_verts = face2verts[index];
            ++ index;
            for (E_Int ivert = 0; ivert < nb_verts; ++ivert, ++index)
            {
                m_face2verts.push_back(face2verts[index]-1);
            }
        }
#       if defined(DEBUG_VERBOSE)
        std::cout << "connectivité face vers sommets : \n";
        for (E_Int iface = 0; iface < nb_faces; ++iface )
        {
            std::cout << "\tface n°" << iface << " : ";
            for ( E_Int ivert = m_beg_face2verts[iface]; ivert < m_beg_face2verts[iface+1]; ++ivert )
            {
                std::cout << m_face2verts[ivert] << " ";
            }
            std::cout << std::endl;
        }
#       endif
        // Reconstruction de la connectivité éléments vers faces :
        // D'abord le tableau d'indirection pour les éléments vers les faces :
        std::vector<E_Int>(nb_elts+1).swap(m_beg_elt2faces);
        m_beg_elt2faces[0] = 0;
        index = 0;
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt)
        {
            E_Int nfaces = elt2Faces[index];
            index += nfaces + 1;
            m_beg_elt2faces[ielt+1] = m_beg_elt2faces[ielt] + nfaces;
        }
        // Puis la connectivité éléments vers faces :
        m_elt2faces.reserve(m_beg_elt2faces[nb_elts]);
        index = 0;
        for ( E_Int ielt = 0; ielt < nb_elts; ++ielt )
        {
            E_Int nfaces = elt2Faces[index];
            ++ index;
            for ( E_Int iface = 0; iface < nfaces; ++iface, ++index )
                m_elt2faces.push_back(elt2Faces[index]-1);
        }
#       if defined(DEBUG_VERBOSE)
        std::cout << "connectivité élément vers faces : " << std::endl;
        for ( E_Int ielt = 0; ielt < nb_elts; ++ielt )
        {
            std::cout << "\nélément n°" << ielt << " : ";
            for ( E_Int iface = m_beg_elt2faces[ielt]; iface < m_beg_elt2faces[ielt+1]; ++iface )
            {
                std::cout << m_elt2faces[iface] << " ";
            } 
            std::cout << std::endl;
        }
#       endif
        // Pour le reste de la connectivité, il va falloir calculer :
        // la connectivité éléments vers sommets :
        // En premier : calcul du tableau d'indirection
        std::vector<E_Int>(nb_elts+1).swap(m_beg_elt2verts);
        m_beg_elt2verts[0] = 0;
        for ( E_Int ielt = 0; ielt < nb_elts; ++ielt )
        {
            std::unordered_set<E_Int> set_of_vertex_indices;
            for (E_Int iface = m_beg_elt2faces[ielt]; iface < m_beg_elt2faces[ielt+1]; ++iface )
            {
                E_Int index_face = m_elt2faces[iface];
                set_of_vertex_indices.insert(m_face2verts.begin()+m_beg_face2verts[index_face  ],
                                             m_face2verts.begin()+m_beg_face2verts[index_face+1]);
            }
            m_beg_elt2verts[ielt+1] = m_beg_elt2verts[ielt] + set_of_vertex_indices.size();
        }
        // Puis calcul du tableau de connectivité élément vers sommets :
        std::vector<E_Int>(m_beg_elt2verts[nb_elts]).swap(m_elt2verts);
#       pragma omp parallel for        
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt )
        {
            std::unordered_set<E_Int> set_of_vertex_indices;
            for (E_Int iface = m_beg_elt2faces[ielt]; iface < m_beg_elt2faces[ielt+1]; ++iface )
            {
                E_Int index_face = m_elt2faces[iface];
                set_of_vertex_indices.insert(m_face2verts.begin()+m_beg_face2verts[index_face  ],
                                             m_face2verts.begin()+m_beg_face2verts[index_face+1]);
            }
            E_Int ipos = 0;
            for ( E_Int index_vertex : set_of_vertex_indices )
            {
                m_elt2verts[m_beg_elt2verts[ielt]+ipos] = index_vertex;
                ++ipos;
            }
        }
#       if defined(DEBUG_VERBOSE)
        auto coords = this->getCoordinates();
        std::cout << "Connectivité élément vers sommets : " << std::endl;
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt )
        {
            std::cout << "Elément n°" << ielt << " : ";
            for ( E_Int ivert = m_beg_elt2verts[ielt]; ivert < m_beg_elt2verts[ielt+1]; ++ivert )
                std::cout << m_elt2verts[ivert] << " ";
            std::cout << std::endl;
            // On va calculer l'aire des tétraèdres pour vérifier qu'ils sont non nuls à l'aide de la formule de Tartaglia (tétraèdre quelconque)
            if (m_beg_elt2verts[ielt+1]-m_beg_elt2verts[ielt] == 4)
            {
                E_Int beg_ind = m_beg_elt2verts[ielt];
                E_Int ind1    = m_elt2verts[beg_ind+0];
                E_Int ind2    = m_elt2verts[beg_ind+1];
                E_Int ind3    = m_elt2verts[beg_ind+2];
                E_Int ind4    = m_elt2verts[beg_ind+3];

                point3d p1{coords[0][ind1],coords[1][ind1],coords[2][ind1]};
                point3d p2{coords[0][ind2],coords[1][ind2],coords[2][ind2]};
                point3d p3{coords[0][ind3],coords[1][ind3],coords[2][ind3]};
                point3d p4{coords[0][ind4],coords[1][ind4],coords[2][ind4]};

                vector3d v1(p1,p2);
                double d12 = norm(v1);
                vector3d v2(p1,p3);
                double d13 = norm(v2);
                vector3d v3(p1,p4);
                double d14 = norm(v3);
                vector3d v4(p2,p3);
                double d23 = norm(v4);
                vector3d v5(p3,p4);
                double d34 = norm(v5);
                vector3d v6(p2,p4);
                double d24 = norm(v6);

                // CAlcul du carré du volume
                double V2 = (d12*(-d12*d34 +d13*(-d23 + d24 + d34) + d14*(d23 -d24 +d34) + d34*(d23 + d24 - d34)) 
                             + d13* (-d13*d24 +d14*d23 +d14*d24 -d14*d34 +d23*d24 -d24*d24 +d24*d34) 
                             + d14* (-d14*d23 - d23*d23 +d23*d24 +d23*d34) - d23*d24*d34)/144.;

                if ((V2<0) || (std::sqrt(V2) < 1.E-14))
                { // A ce point, on considère le volume nul, non ?
                    if (V2 > 0)
                        std::cerr << "L'élément n°" << ielt << " a un volume de " << std::sqrt(V2) << " ce qui est trop faible !" << std::endl;
                    else 
                        std::cerr << "L'élément n°" << ielt << " a un volume au carré négatif ????" << std::endl;
                    std::cerr << "Ces sommets sont : " << std::string(p1) << ", " << std::string(p2) << ", "
                              << std::string(p3) << ", " << std::string(p4) << std::endl;
                    std::cerr << "Indices des points considérés : " << ind1 << ", " << ind2 << ", " << ind3 << ", " << ind4 << std::endl;
                    std::cerr << "Indice des Faces composant l'élément (donnée brute du maillage) : ";
                    for ( E_Int iface = m_beg_elt2faces[ielt]; iface < m_beg_elt2faces[ielt+1]; ++iface )
                    {
                        std::cerr << m_elt2faces[iface] << " ";
                    } 
                    std::cerr << std::endl;
                    std::cerr << "Coordonnées des sommets des faces (pour voir si c'est bien consistant) : " << std::endl;
                    for ( E_Int iface = m_beg_elt2faces[ielt]; iface < m_beg_elt2faces[ielt+1]; ++iface )
                    {
                        E_Int ind_face = m_elt2faces[iface];
                        std::cerr << "Face " << ind_face << " : [ ";
                        for ( E_Int ivert = m_beg_face2verts[ind_face]; ivert < m_beg_face2verts[ind_face+1]; ivert++)
                        {
                            E_Int ind_vert = m_face2verts[ivert];
                            std::cerr << std::string(point3d{coords[0][ind_vert],coords[1][ind_vert],coords[2][ind_vert]}) << " ";
                        }
                        std::cerr << "]" << std::endl;
                    }                    
                }
            }
        }
#       endif
        // Calcul de la connectivité sommet vers éléments.
        // Comme toujours, on va commencer par le tableau d'indirection, parque qu'il le vaut bien :-p
        E_Int nb_verts = fields->getSize();
        std::vector<E_Int>(nb_verts+1).swap(m_beg_vert2elts);
        std::vector<E_Int> nb_elts_per_vert(nb_verts,0);
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt)
        {
            for (E_Int ivert = m_beg_elt2verts[ielt]; ivert < m_beg_elt2verts[ielt+1]; ++ivert)
                nb_elts_per_vert[m_elt2verts[ivert]] += 1;
        }        
        std::vector<E_Int>(nb_verts+1).swap(m_beg_vert2elts);
        m_beg_vert2elts[0] = 0;
        for (E_Int ivert = 0; ivert < nb_verts; ++ivert)
        {
            m_beg_vert2elts[ivert+1] = m_beg_vert2elts[ivert] + nb_elts_per_vert[ivert];
        }
        m_vert2elts.resize(m_beg_vert2elts[nb_verts]);
        // Et calcul de la connectivité sommet vers éléments :
        nb_elts_per_vert.assign(nb_elts_per_vert.size(), 0);
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt)
        {
            for (E_Int ivert = m_beg_elt2verts[ielt]; ivert < m_beg_elt2verts[ielt+1]; ++ivert)
            {
                E_Int index_vert = m_elt2verts[ivert];
                m_vert2elts[m_beg_vert2elts[index_vert]+nb_elts_per_vert[index_vert]] = ielt;
                nb_elts_per_vert[index_vert] += 1;
            }
        }
#       if defined(DEBUG_VERBOSE)
        std::cout << "Connectivité sommet vers éléments : " << std::endl;
        for (E_Int ivert = 0; ivert < nb_verts; ++ivert )
        {
            std::cout << "Sommet n°" << ivert << " : ";
            for (E_Int ielt = m_beg_vert2elts[ivert]; ielt < m_beg_vert2elts[ivert+1]; ielt++)
                std::cout << m_vert2elts[ielt] << " ";
            std::cout << std::endl;
        }
#       endif
        // Dernière connectivité, mais curieusement pas la moindre : celle des faces vers les éléments
        // La difficulté ici est surtout de trouver si la définition de la face donnée en entrée
        // est une face directe ou indirecte (sortante ou rentrante pour la normale de la face) pour
        // l'élément considéré. Si cette face est directe, on donne le numéro de l'élément (commençant pour ce tableau
        // exceptionnellement à 1), et sinon si la face est indirecte, on donne l'opposé du numéro de l'élément (donc un
        // indice négatif).
        // 
        // Pour trouver si la face est directe ou indirecte, je propose l'algorithme suivant :
        // Puisqu'on parcourt la connectivité élément vers faces pour construire face vers éléments,
        // on peut considérer toutes les faces d'un élément. On prend la première face définie pour
        // l'élément, on calcule sa normale puis on détecte le nombre d'intersection de cette normale
        // avec les autres faces (grâce à l'algorithme rapide mis en oeuvre dans la structure face).
        // Si cette intersection est impaire, c'est que la face est dans le sens indirect (la normale rentre
        // dans le n-gon), sinon c'est que la face est dans le sens direct.
        // Pour les autres faces, on regarde les arêtes communes avec notre face (si il y en a),
        // et on regarde si la numérotation des sommets est dans le même ordre, ce qui veut dire que le sens de cette face
        // est inversée par rapport à la face de référence, ou dans l'ordre inverse, ce qui veut dire que le sens de cette
        // face est la même que la face de référence. Si la face n'a pas d'arête commune avec notre face de référence, on la
        // met dans un pool de faces non visitées. On prend ensuite la face visitée suivante et on teste de nouveau
        // quelles sont les faces non visitées ayant une arête commune avec cette nouvelle face de référence et détecter
        // de la même manière leurs orientations. Puis on itère de nouveau sur les faces visitées par rapport aux faces
        // non visitées jusqu'à ce qu'il n'y ait plus de faces non visitées. Normalement, puisqu'on peut raisonnablement
        // penser qu'un n-gon est connexe (mais pas forcément convexe), et que le maillage est conforme, toutes
        // les faces auront ainsi été visitées et détectées comme directes ou indirectes à la fin des itérations
        // pour l'élément donné. Il ne reste plus qu'à remplir correctement le tableau de connectivité face2elts
        // avec cette nouvelle information.
        std::vector<std::pair<E_Int,E_Int>>(nb_faces,{0,0}).swap(m_face2elts);
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt) // Pour chaque élément :
        {
            E_Int nfaces = m_beg_elt2faces[ielt+1] - m_beg_elt2faces[ielt];
            std::vector<face> faces; faces.reserve(nfaces);
            std::vector<bool> is_direct(nfaces);
            for (E_Int iface = m_beg_elt2faces[ielt]; iface < m_beg_elt2faces[ielt+1]; ++iface)
            {
                E_Int index_face = m_elt2faces[iface];
                std::vector<E_Int> indCoords(m_face2verts.begin()+m_beg_face2verts[index_face  ],
                                             m_face2verts.begin()+m_beg_face2verts[index_face+1]);
                faces.emplace_back( indCoords, std::pair<E_Int,E_Int>{-1,-1}, this->getCoordinates() );
            }
            E_Int ref_face = 0;
            bool found_ref_face = false;
            if (nfaces > 4)
            {
                // Calcul de la bouding box alignée:
                double xmin, ymin, zmin;
                auto coords = this->getCoordinates();
                E_Int nb_verts = m_beg_elt2verts[ielt+1] - m_beg_elt2verts[ielt];
                E_Int beg_verts = m_beg_elt2verts[ielt];
                xmin = coords[0][m_elt2verts[beg_verts]];
                ymin = coords[1][m_elt2verts[beg_verts]];
                zmin = coords[2][m_elt2verts[beg_verts]];
                for (E_Int ivert = 1; ivert < nb_verts; ++ivert )
                {
                    E_Int ind_vert = m_elt2verts[beg_verts+ivert];
                    xmin = std::min(coords[0][ind_vert],xmin);
                    ymin = std::min(coords[1][ind_vert],xmin);
                    zmin = std::min(coords[2][ind_vert],xmin);
                }
                while (!found_ref_face)
                {
                    const point3d &    origin = faces[ref_face].get_barycenter();
                    const vector3d& direction = faces[ref_face].get_normal();
                    vector3d normal = (1./abs(direction))*direction;
                    point3d  o = origin + (-1.E-8)*normal;
                    point3d  r{xmin-1.E-3,o.y, o.z};
                    vector3d d{o.x-r.x,0.,0.};
                    E_Int nb_intersects = 0;
                    bool is_intersecting, is_entering;
                    try { // Normalement, dans 99.99% des cas, tout devrait bien se passer
                        for (E_Int iface = 0; iface < nfaces; ++iface )
                        {                            
                            std::tie(is_intersecting,is_entering) = faces[iface].is_intersecting_ray(r, d);
                            if (is_intersecting)
                            {
                                const point3d &    o = faces[iface].get_barycenter();
                                const vector3d&    df= faces[iface].get_normal();
                                vector3d n = (1./abs(df))*df;
                                vector3d ofr(o, r);
                                // On va calculer le point d'intersection :
                                // Plan "moyen" de la face : (o,df)
                                // Rayon : (r, d) où d = vector3d(ro)
                                double nxdx = n.x*d.x;
                                if (std::abs(nxdx) < 1.E-14) throw std::underflow_error("Face parallèle à la direction");
                                double t = -(ofr|n)/nxdx;
                                if ((0<=t) && (t<=1)) nb_intersects += 1;  
                            } 
                        }
                        is_direct[ref_face] = (nb_intersects%2==1);// Si impair, point à l'intérieur, donc face bien orientée
                        found_ref_face = true;
                    }
                    catch(std::underflow_error& err) // Mais ça arrive qu'on croise des arêtes ou des sommets...
                    { 
                        // On essaie du coup dans la direction y :
                        r = point3d{o.x,ymin-1.E-3,o.z};
                        d = vector3d{0.,o.y-r.y,0.};
                        try { 
                            for (E_Int iface = 0; iface < nfaces; ++iface )
                            {
                                std::tie(is_intersecting,is_entering) = faces[iface].is_intersecting_ray(r, d);
                                if (is_intersecting)
                                {
                                    const point3d &    o = faces[iface].get_barycenter();
                                    const vector3d&    df= faces[iface].get_normal();
                                    vector3d n = (1./abs(df))*df;
                                    vector3d ofr(o, r);
                                    // On va calculer le point d'intersection :
                                    // Plan "moyen" de la face : (o,df)
                                    // Rayon : (r, d) où d = vector3d(ro)
                                    double nydy = n.y*d.y;
                                    if (std::abs(nydy) < 1.E-14) throw std::underflow_error("Face parallèle à la direction");
                                    double t = -(ofr|n)/nydy;
                                    if ((0<=t) && (t<=1)) nb_intersects += 1;  
                                } 
                            }
                            is_direct[ref_face] = (nb_intersects%2==1);// Si impair, point à l'intérieur, donc face bien orientée
                            found_ref_face = true;
                        }
                        catch(std::underflow_error& err)
                        {
                            // Dans la direction z ?
                            r = point3d{o.x,o.y,zmin-1.E-3};
                            d = vector3d{0.,0.,o.z-r.z};
                            try {
                                for (E_Int iface = 0; iface < nfaces; ++iface )
                                {
                                    std::tie(is_intersecting,is_entering) = faces[iface].is_intersecting_ray(r, d);
                                    if (is_intersecting)
                                    {
                                        const point3d &    o = faces[iface].get_barycenter();
                                        const vector3d&    df= faces[iface].get_normal();
                                        vector3d n = (1./abs(df))*df;
                                        vector3d ofr(o, r);
                                        // On va calculer le point d'intersection :
                                        // Plan "moyen" de la face : (o,df)
                                        // Rayon : (r, d) où d = vector3d(ro)
                                        double nzdz = n.z*d.z;
                                        if (std::abs(nzdz) < 1.E-14) throw std::underflow_error("Face parallèle à la direction");
                                        double t = -(ofr|n)/nzdz;
                                        if ((0<=t) && (t<=1)) nb_intersects += 1;  
                                    } 
                                }
                                is_direct[ref_face] = (nb_intersects%2==1);// Si impair, point à l'intérieur, donc face bien orientée
                                found_ref_face = true;
                            }
                            catch(std::underflow_error& err)
                            {
                                // Pour la face courante, on est dans une configuration indéterminée, donc je prend une autre face
                                // comme face de référence.
                                ref_face ++;
                                assert(ref_face < nfaces);
                            }// catch
                        }// catch
                    }// catch
                }// while (!found_ref_face)
            }// nfaces > 4...
            else 
            {   // Pour les tétraèdres, on va faire un truc numériquement plus stable ?
                // Calcul barycentre du barycentre légèrement perturbé par la normale :
                //std::cerr << "Testing with barycenters." << std::endl;
                auto coords = this->getCoordinates();
                const point3d &    origin = faces[ref_face].get_barycenter();
                const vector3d& direction = faces[ref_face].get_normal();
                // On récupère les sommets du tétraèdre et on calcul le barycentre du tétraèdre :
                assert(m_beg_elt2verts[ielt+1]-m_beg_elt2verts[ielt] == 4);
                E_Int ivert = m_beg_elt2verts[ielt];
                std::array<point3d,4> sommets;
                point3d barycentre{0.,0.,0.};
                for ( E_Int i = 0; i < 4; ++i)
                {
                    E_Int ind_vert = m_elt2verts[ivert+i];
                    sommets[i] = point3d{coords[0][ind_vert],coords[1][ind_vert],coords[2][ind_vert]};
                    barycentre.x += 0.25*sommets[i].x;
                    barycentre.y += 0.25*sommets[i].y;
                    barycentre.z += 0.25*sommets[i].z;
                }
                vector3d sortante(barycentre, origin);
                if (std::abs((sortante|direction)) < 1.E-16) std::cout << "Warning: too much small dot product..." << std::endl;
                if ((sortante|direction) > 0)
                    is_direct[ref_face] = true;
                else
                    is_direct[ref_face] = false;
            }

            /**
             * @brief      Fonction lambda retournant 0 si par d'arêtes communes (donc cas encore indéterminé), 
             !                                        1 si même orientation
             !                                       -1 si orientation opposée
             *
             * @param[in]  f1    La première face
             * @param[in]  f2    La seconde face
             *
             * @return     L'orientation de la facette f2 (1 ou -1) ou cas indéterminé (0)
             */
            auto same_orientation = [] (const face& f1, const face& f2)
            {
                E_Int pos_first_common_vertex_f1 = -1, pos_first_common_vertex_f2 = -1;
                E_Int nb_verts1 = f1.indices_vertices.size();
                E_Int nb_verts2 = f2.indices_vertices.size();

                for (E_Int ipos1 = 0; (ipos1 < nb_verts1) && (pos_first_common_vertex_f1==-1); ++ipos1 )
                {
                    for (E_Int ipos2 = 0; (ipos2 < nb_verts2) && (pos_first_common_vertex_f2==-1); ++ipos2)
                    {
                        if (f1.indices_vertices[ipos1] == f2.indices_vertices[ipos2])
                        {
                            pos_first_common_vertex_f1 = ipos1;
                            pos_first_common_vertex_f2 = ipos2;
                        }
                    }
                }
                if (pos_first_common_vertex_f1 == -1)
                {
                    assert(pos_first_common_vertex_f2 == -1);
                    return 0;
                }
                assert(pos_first_common_vertex_f2 != -1);
                if (pos_first_common_vertex_f1 == nb_verts1-1) // Dans ce cas, qu'un seul sommet commun, cas indéterminé
                    return 0;
                E_Int next_index_f1 = f1.indices_vertices[pos_first_common_vertex_f1+1];
                if (next_index_f1 == f2.indices_vertices[(pos_first_common_vertex_f2+1)%nb_verts2])
                    return -1;
                if (next_index_f1 == f2.indices_vertices[(pos_first_common_vertex_f2-1+nb_verts2)%nb_verts2])
                    return +1;
                E_Int prev_index_f1 = f1.indices_vertices[(pos_first_common_vertex_f1-1+nb_verts1)%nb_verts1];
                if (prev_index_f1 == f2.indices_vertices[(pos_first_common_vertex_f2+1)%nb_verts2])
                    return +1;
                if (prev_index_f1 == f2.indices_vertices[(pos_first_common_vertex_f2-1+nb_verts2)%nb_verts2])
                    return -1;
                // Un seul sommet commun, on retourne une indétermination :
                return 0;
            };
            // Création d'une liste de face dont on a déjà trouvé l'orientation
            std::list<E_Int> is_visited; is_visited.push_back(ref_face);
            auto iter_visited = is_visited.begin();
            // Création d'un pool de faces dont l'orientation est encore inconnue
            std::list<E_Int> pool;
            for (E_Int iface = 0; iface < nfaces; ++iface)
                if (iface != ref_face)
                    pool.push_back(iface);

            while (!pool.empty()) // Tant qu'il y a encore des candidats à tester pour leur orientation :
            {
                for (auto iter = pool.begin(); iter != pool.end(); )
                {
                    E_Int res = same_orientation(faces[*iter_visited], faces[*iter]);
                    if (res == 1) // Cas même orientation :
                    {
                        is_direct[*iter] = is_direct[*iter_visited];
                        is_visited.push_back(*iter);
                        iter = pool.erase(iter);

                    }
                    if (res == -1) // Cas orientation opposée :
                    {
                        is_direct[*iter] = !is_direct[*iter_visited];
                        is_visited.push_back(*iter);
                        iter = pool.erase(iter);
                    }
                    if (res == 0) ++ iter; // Cas indéterminé, on passe à la prochaine face
                }// Fin du test des faces encore dans le pool
                if (!pool.empty()) // Si il y a encore des faces dans le pool
                {   // On prend la prochaine face déjà visitée comme face de référence
                    iter_visited ++;
                    assert(iter_visited != is_visited.end());
                }
            }// Tant que le pool n'est pas vide
            // On remplit la connectivité face2elts avec l'orientation trouvée :
            for (E_Int iface = 0; iface < nfaces; ++iface )
            {
                E_Int index_face = m_elt2faces[m_beg_elt2faces[ielt]+iface];
                if (m_face2elts[index_face].first == 0)
                {
                    if (is_direct[iface]) m_face2elts[index_face].first =   ielt+1;
                    else                  m_face2elts[index_face].first = -(ielt+1);
                }
                else
                {
                    assert(m_face2elts[index_face].second == 0);
                    if (is_direct[iface]) m_face2elts[index_face].second =   ielt+1;
                    else                  m_face2elts[index_face].second = -(ielt+1);

                }
            }
        }// End for (ielt = ) => Fin de la boucle sur les éléments
#       if defined(DEBUG_VERBOSE)
        std::cout << "Connectivité face vers éléments : (négatif = indirect ) : " << std::endl;
        for ( E_Int iface = 0; iface < nb_faces; ++iface )
        {
            std::cout << "\tface n°" << iface << " : {";
            std::cout << m_face2elts[iface].first << ", " << m_face2elts[iface].second << "}" << std::endl;
        }
#       endif
    }
    //_ _______________________________ Retourne les faces constituant un élément donné __________________
    std::vector<face> ngon_data_view::Implementation::get_faces_of_element( E_Int number, E_Int no_zone ) const
    {
        E_Int beg_e2f = m_beg_elt2faces[number  ];
        E_Int end_e2f = m_beg_elt2faces[number+1];
        E_Int nb_faces_per_elt = end_e2f - beg_e2f;
        std::vector<face> faces; faces.reserve(nb_faces_per_elt);
        for (E_Int iface = 0; iface < nb_faces_per_elt; ++iface)
        {
            E_Int index_face = m_elt2faces[beg_e2f+iface];
            E_Int beg_vertices = m_beg_face2verts[index_face  ];
            E_Int end_vertices = m_beg_face2verts[index_face+1];
            const auto& f2e = m_face2elts[index_face];
            E_Int afirst = std::abs(f2e.first)-1;
            E_Int cur_elt= ( number == afirst? f2e.first : f2e.second );
            E_Int op_elt = ( number == afirst ? f2e.second : f2e.first );
            op_elt = (op_elt == 0 ? -1 : op_elt > 0 ? op_elt - 1 : - 1 - op_elt);
            std::vector<E_Int> coords;
            if (cur_elt > 0) // OK, la facette est bien dans le sens direct :
                coords.assign(m_face2verts.begin()+beg_vertices, m_face2verts.begin()+end_vertices);
            else
            {   // Ah non, face stockée avec normale rentrante pour cet élément. Donc on va
                // recréer la face en parcourant ses sommets à l'envers :
                coords.reserve(end_vertices-beg_vertices);
                for (E_Int ivert = end_vertices-1; ivert >= beg_vertices; --ivert)
                    coords.push_back(m_face2verts[ivert]);
            }
            faces.emplace_back(coords, std::pair<E_Int,E_Int>{number,op_elt}, this->getCoordinates());
        }
        return faces;
    }
    //_ ____________________________ Test si un élément donné contient un point donné ____________________
    bool 
    ngon_data_view::Implementation::is_containing(E_Int ind_elt, const point3d& pt) const
    {
        const auto& crds = this->getCoordinates();
        std::array<std::vector<double>,3> coords; // Nombre de sommet dépend du type d'élément (avec barycentre pour certains)
        E_Int nb_verts = this->m_beg_elt2verts[ind_elt+1] - this->m_beg_elt2verts[ind_elt];
        // On va calculer le nombre de barycentres :
        // Remarque :
        // =========
        // A la lecture de l'interface de la triangulation de Delaunay (pour le cas où le nombre de sommets de la facette
        // est plus grand que 4), il semble qu'aucun nouveau point est généré. Donc nous n'avons qu'à compter que le
        // nombre de barycentres créés pour le cas des quadrangles.
        E_Int nb_barycenters = 0;
        E_Int nb_faces = m_beg_elt2faces[ind_elt+1] - m_beg_elt2faces[ind_elt];
        for (E_Int iface = 0; iface < nb_faces; ++iface)
        {
            E_Int index_face = m_elt2faces[m_beg_elt2faces[ind_elt] + iface];
            E_Int nb_vertices_for_face = m_beg_face2verts[index_face+1] - m_beg_face2verts[index_face];
            if (nb_vertices_for_face == 4) nb_barycenters ++;
        } 
        coords[0].reserve(nb_verts+nb_barycenters); coords[1].reserve(nb_verts+nb_barycenters); 
        coords[2].reserve(nb_verts+nb_barycenters);
        std::unordered_map<E_Int, E_Int> glob2loc;
        // On extrait tous les points de la cellule :
        double xb = 0, yb = 0, zb = 0;
        for (E_Int ivert = 0; ivert < nb_verts; ++ivert)
        {
            E_Int ind_vert = m_elt2verts[this->m_beg_elt2verts[ind_elt]+ivert];
            glob2loc[ind_vert] = ivert;
            xb += crds[0][ind_vert];
            yb += crds[1][ind_vert];
            zb += crds[2][ind_vert];
            coords[0].push_back(crds[0][ind_vert]);
            coords[1].push_back(crds[1][ind_vert]);
            coords[2].push_back(crds[2][ind_vert]);
        }
        point3d bary(xb/nb_verts,yb/nb_verts, zb/nb_verts);
        // On doit utiliser des triangles en tessalisant les faces de l'élément :
        using triangle_type = triangulated_polyhedron::triangle_type;
        // Par contre, il me semble difficile de connaître le nombre de triangle d'avance...
        // Je vais calculer un nombre minimal de triangles qu'on doit avoir :
        E_Int nb_triangles = 0;
        for (E_Int iface = 0; iface < nb_faces; ++iface)
        {
            E_Int index_face = m_elt2faces[m_beg_elt2faces[ind_elt] + iface];
            E_Int nb_vertices_for_face = m_beg_face2verts[index_face+1] - m_beg_face2verts[index_face];
            if (nb_vertices_for_face == 3) nb_triangles += 1;
            else if (nb_vertices_for_face == 4) nb_triangles += 4;
            else 
            {
                assert(nb_vertices_for_face > 4);
                // On va prendre comme nombre de triangles : nb sommets - 2
                nb_triangles += nb_vertices_for_face-2;
            }
        }
        std::vector<triangle_type> faces; faces.reserve(nb_triangles);
        // On parcourt les faces du n-gon :
        for ( E_Int iface = 0; iface < nb_faces; ++iface)
        {
            E_Int index_face = m_elt2faces[m_beg_elt2faces[ind_elt] + iface];
            E_Int nb_vertices_for_face = m_beg_face2verts[index_face+1] - m_beg_face2verts[index_face];
            if (nb_vertices_for_face == 3)
            {   // C'est un triangle, on l'inclu directement :
                // En faisant attention à l'orientation !!!!!
                auto f2e = m_face2elts[index_face];
                E_Int ind1, ind2, ind3;
                if ((f2e.first == ind_elt+1) ||(f2e.second == ind_elt+1))
                { // Sens direct
                    ind1 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+0]];
                    ind2 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+1]];
                    ind3 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+2]];
                }
                else
                {   // Sens indirect, on permute les indices
                    assert((-f2e.first == ind_elt+1) ||(-f2e.second == ind_elt+1));
                    ind1 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+0]];
                    ind2 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+2]];
                    ind3 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+1]];                    
                }
                point3d p1{coords[0][ind1], coords[1][ind1],coords[2][ind1]};
                point3d p2{coords[0][ind2], coords[1][ind2],coords[2][ind2]};
                point3d p3{coords[0][ind3], coords[1][ind3],coords[2][ind3]};
                point3d barf = point3d{(p1.x+p2.x+p3.x)/3.,(p1.y+p2.y+p3.y)/3.,(p1.z+p2.z+p3.z)/3.};
                vector3d ob(bary,barf);
                vector3d no = vector3d{p1,p2} ^ vector3d{p1,p3};
                if ((no|ob) < 0)
                {
                    std::cerr << "Facette mal orientée !";
                    std::cerr << "ind1 : " << ind1 << " ind2 : " << ind2 << ", ind3 : " << ind3 << std::endl;

                    std::cerr << "Sommets du tetrahèdre : " << std::endl;
                    for (E_Int ivert = 0; ivert < nb_verts; ++ivert)
                        printf("{%13.9g,%13.9g,%13.9g}",coords[0][ivert],coords[1][ivert],coords[2][ivert]);
                    std::cerr << "Sommets orientés de la face : " << std::endl;
                    std::cerr << std::string(p1) << ", " << std::string(p2) << ", " << std::string(p3) << std::endl;
                    std::cerr << "no : " << std::string(no) << std::endl;
                    std::cerr << "ob : " << std::string(ob) << std::endl;
                    std::cerr << "<no|ob> = " << (no|ob) << std::endl;
                }
                faces.emplace_back(triangle_type{ind1, ind2, ind3}); // On rajoute directement la facette
            }
            else if (nb_vertices_for_face == 4)
            {   // C'est un quadrangle. On va devoir rajouter le barycentre pour
                // tessaliser sans être ambigu entre éléments.
                // Attention à l'orientation :
                auto f2e = m_face2elts[index_face];
                E_Int ind1 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+0]];
                E_Int ind2 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+1]];
                E_Int ind3 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+2]];
                E_Int ind4 = glob2loc[m_face2verts[m_beg_face2verts[index_face]+3]];
                E_Int ind5 = coords[0].size(); // Index du point barycentre mis à la fin de coords
                // Calcul d'un barycentre et rajout d'un point à coords :
                coords[0].push_back(0.25*(coords[0][ind1]+coords[0][ind2]+coords[0][ind3]+coords[0][ind4]));
                coords[1].push_back(0.25*(coords[1][ind1]+coords[1][ind2]+coords[1][ind3]+coords[1][ind4]));
                coords[2].push_back(0.25*(coords[2][ind1]+coords[2][ind2]+coords[2][ind3]+coords[2][ind4]));
                if ((f2e.first == ind_elt+1) ||(f2e.second == ind_elt+1))
                { // Sens direct
                    faces.emplace_back(triangle_type{ind1, ind2, ind5});
                    faces.emplace_back(triangle_type{ind2, ind3, ind5});
                    faces.emplace_back(triangle_type{ind3, ind4, ind5});
                    faces.emplace_back(triangle_type{ind4, ind1, ind5});
                }
                else
                {
                    // Sens indirect, on permute des indices dans les triangles :
                    faces.emplace_back(triangle_type{ind2, ind1, ind5});
                    faces.emplace_back(triangle_type{ind3, ind2, ind5});
                    faces.emplace_back(triangle_type{ind4, ind3, ind5});
                    faces.emplace_back(triangle_type{ind1, ind4, ind5});                    
                }
            }
            else
            {
                assert(nb_vertices_for_face > 4);
                DELAUNAY::Triangulator dt;
                K_FLD::FloatArray coord(nb_verts, 3);
                std::vector<E_Int> nodes(nb_verts);
                for ( E_Int i_vert = 0; i_vert < nb_vertices_for_face; ++i_vert )
                {
                    E_Int loc_ind_vert = glob2loc[m_face2verts[m_beg_face2verts[index_face]+i_vert]];
                    coord(i_vert, 0) = coords[0][loc_ind_vert];
                    coord(i_vert, 1) = coords[0][loc_ind_vert];
                    coord(i_vert, 2) = coords[0][loc_ind_vert];

                    nodes[i_vert] = i_vert;
                 }
                K_FLD::IntArray cT3;
                //E_Int err = 
                K_MESH::Polygon::triangulate<DELAUNAY::Triangulator>(dt, coord, nodes.data(), nb_vertices_for_face, 0, cT3);

                auto f2e = m_face2elts[index_face];
                if ((f2e.first == ind_elt+1) ||(f2e.second == ind_elt+1))
                { // Sens direct
                    for ( E_Int icol = 0; icol < cT3.cols(); ++icol)
                    {
                        E_Int* p_ct = cT3.col(icol);
                        faces.emplace_back(triangle_type{p_ct[0],p_ct[1], p_ct[2]});
                    }
                }
                else
                {
                    for ( E_Int icol = 0; icol < cT3.cols(); ++icol)
                    {
                        E_Int* p_ct = cT3.col(icol);
                        faces.emplace_back(triangle_type{p_ct[1],p_ct[0], p_ct[2]});
                    }                    
                }
            }
        }
        // Construction du polyèdre :
        triangulated_polyhedron polyedre( faces, {
                                     K_MEMORY::vector_view<const double>(coords[0].begin(),coords[0].end()),
                                     K_MEMORY::vector_view<const double>(coords[1].begin(),coords[1].end()),
                                     K_MEMORY::vector_view<const double>(coords[2].begin(),coords[2].end())
                                                 } );
        bool is_inside;
        try {
            is_inside = polyedre.is_containing(pt);
        } catch(std::invalid_argument& err)
        {
            // On est sur la frontière de l'élément :
            //std::cerr << "Warning: interpolated point is on interface. Possibility to have two points in same location in the stream line"
            //          << std::flush << std::endl;
            is_inside = true; // Dans ce cas, on considère qu'on est à l'intérieur (on prend l'élément comme un fermé topologique)
        }
        return is_inside;

    }
    //_ ___________________________ Retourne la cellule contenant un point donné _________________________
    E_Int 
    ngon_data_view::Implementation::get_interpolation_cell( const point3d& point ) const
    {
        if (not this->aabbox.contains(point)) return -1;
        E_Int ind_nearest_vertex;
        double dist_nearest_vertex;
        std::tie(ind_nearest_vertex, dist_nearest_vertex) = this->tree.nearest(point);
        // On va chercher l'appartenance du point aux cellules contenant le point le plus proche :
        for (E_Int ielt = m_beg_vert2elts[ind_nearest_vertex]; ielt < m_beg_vert2elts[ind_nearest_vertex+1]; ++ielt )
        {
            E_Int elt = this->m_vert2elts[ielt];
            bool found = this->is_containing(elt, point);
            if (found) return elt;
        }
        // Et sinon on cherche chez les voisins des voisins pour une couverture complète de localisation :
        for (E_Int ielt = m_beg_vert2elts[ind_nearest_vertex]; ielt < m_beg_vert2elts[ind_nearest_vertex+1]; ++ielt )
        {
            E_Int elt = this->m_vert2elts[ielt];
            for ( E_Int ivert = m_beg_elt2verts[elt]; ivert < m_beg_elt2verts[elt+1]; ++ivert )
            {
                E_Int index_vert = m_elt2verts[ivert];
                for ( E_Int ielt2 = m_beg_vert2elts[index_vert]; ielt2 < m_beg_vert2elts[index_vert+1]; ielt2++ )
                {
                    E_Int elt2 = this->m_vert2elts[ielt2];
                    bool found = this->is_containing(elt2, point);
                    if (found) return elt2;
                }
            }
        }
        return -1;
    }
    //_ ________________________ Interpole un champ sur un point donné ___________________________________
    void 
    ngon_data_view::Implementation::compute_interpolated_field( const point3d& pt, E_Int ind_cell, 
                                                                E_Int ipos, FldArrayF& interpolatedField ) const
    {
        // Ne connaissant le nombre de sommets qu'à l'exécution, l'idée ici est de faire une simple
        // interpolation linéaire en prenant le tétraèdre formé par les quatre sommets les plus proches
        // du point d'interpolation.
        E_Int nb_verts_per_elt = this->m_beg_elt2verts[ind_cell+1] - this->m_beg_elt2verts[ind_cell];
        const auto& crds = this->getCoordinates();
        std::vector<point3d> coords;// Nombre de sommet dépend du type d'élément (avec barycentre pour certains)
        coords.reserve(nb_verts_per_elt);
        FldArrayF* fld = this->fields;
        E_Int nfld = fld->getNfld();
        std::vector<std::vector<double>> values(nfld);
        for (E_Int ifld = 0; ifld < nfld; ++ifld)
        {
            values[ifld].reserve(nb_verts_per_elt);
        }
        // On extrait tous les points de la cellule selon l'ordre CGNS :
        for (E_Int ivert = 0; ivert < nb_verts_per_elt; ++ivert)
        {
            E_Int ind_vert = this->m_elt2verts[this->m_beg_elt2verts[ind_cell]+ivert];
            coords.emplace_back(crds[0][ind_vert], crds[1][ind_vert], crds[2][ind_vert]);
            for (E_Int ifld = 0; ifld < nfld; ++ifld)
            {
                values[ifld].push_back((*fld)(ind_vert,ifld+1));
            }
        }

        std::vector<E_Int> vertices_index; vertices_index.reserve(nb_verts_per_elt);
        for ( E_Int ivert = 0; ivert < nb_verts_per_elt; ++ivert) vertices_index.push_back(ivert);
        std::sort(vertices_index.begin(), vertices_index.end(), [&](E_Int i, E_Int j) {
            double dip = square_distance(coords[i], pt);
            double djp = square_distance(coords[j], pt);
            return dip < djp;
        });
        // Interpolation sur le tetraèdre T(pᵢ₀, pᵢ₁, pᵢ₂, pᵢ₃)
        //                                        →  →  →
        // On choisit comme repère barycentrique (e₁,e₂,e₃)
        //       →                         →                          →
        // avec  e₁ le vecteur (pᵢ₀, pᵢ₁), e₂ le vecteur (pᵢ₀,pᵢ₂) et e₃ le vecteur (pᵢ₀,pᵢ₃)
        // 
        // On calcule les coordonnées barycentriques (𝛼,𝛽,𝛾) du point p où on doit interpoler : :
        //        →      →      →
        // p₀ + 𝛼.e₁ + 𝛽.e₂ + 𝛾.e₃ = p avec 𝛼 + 𝛽 + 𝛾 ≤ 1, 𝛼 ≥ 0, 𝛽 ≥ 0, 𝛾 ≥ 0
        // 
        // ce qui revient à inverser un système linéaire de dimension trois.
        // Pour calculer le champs interpolé, il suffit alors de calculer :
        // 
        // f(p) = (1-𝛼-𝛽-𝛾).f(pᵢ₀) + 𝛼.f(pᵢ₁) + 𝛽.f(pᵢ₂) + 𝛾.f(pᵢ₃)
        // 
        vector3d e1(coords[vertices_index[0]],coords[vertices_index[1]]);
        vector3d e2(coords[vertices_index[0]],coords[vertices_index[2]]);
        // Pour le dernier sommet, on prend le plus éloigné pour s'assurer d'avoir une matrice inversible 
        // (bien qu'elle puisse être non inversible quand même, faudra réfléchir dans ce cas à choisir le troisième
        // sommet jusqu'à avoir un repère en 3D...)
        vector3d e3(coords[vertices_index[0]],coords[vertices_index[nb_verts_per_elt-1]]);
        vector3d pp0(coords[vertices_index[0]], pt);

        matrix_3x3_type A{std::array<double,3>{e1.x,e2.x,e3.x},
                                              {e1.y,e2.y,e3.y},
                                              {e1.z,e2.z,e3.z}
                         };
        auto barycrds = inverse_linear_system(A, pp0);
        double alpha = barycrds[0];
        double beta  = barycrds[1];
        double gamma = barycrds[2];
        double umabg = 1. - alpha - beta - gamma;
        for ( E_Int ifld = 0; ifld < nfld; ++ifld)
        {
            interpolatedField(ipos,ifld+1) = umabg * values[ifld][vertices_index[0]] + 
                                             alpha * values[ifld][vertices_index[1]] + 
                                             beta  * values[ifld][vertices_index[2]] +
                                             gamma * values[ifld][vertices_index[3]];
        }
    }
    //_ ___________________________ Retourne les indices des sommets de l'élément ________________________
    std::vector<E_Int> ngon_data_view::Implementation::get_indices_of_vertices(E_Int icell) const
    {
        // Ne connaissant le nombre de sommets qu'à l'exécution, l'idée ici est de faire une simple
        // interpolation linéaire en prenant le tétraèdre formé par les quatre sommets les plus proches
        // du point d'interpolation.
        E_Int nb_verts_per_elt = this->m_beg_elt2verts[icell+1] - this->m_beg_elt2verts[icell];
        std::vector<E_Int> ind_verts;
        ind_verts.reserve(nb_verts_per_elt);
        for (E_Int ivert = 0; ivert < nb_verts_per_elt; ++ivert)
            ind_verts.push_back(this->m_elt2verts[this->m_beg_elt2verts[icell]+ivert]);
        return ind_verts;
    }
    //_ ___________________________ Calcul le rotationel pour un élément donné ___________________________
    vector3d ngon_data_view::Implementation::compute_rotational_in_cell(E_Int ind_cell) const
    {
        auto faces = this->get_faces_of_element(ind_cell, 0);
        return compute_rotational_on_ngon(faces, 
                                          {
                                            this->getField(this->pos_velocity[0]),
                                            this->getField(this->pos_velocity[1]),
                                            this->getField(this->pos_velocity[2])
                                          } );
    }
    //_ ___________________________ Calcul le volume d'un élément donné __________________________________
    double ngon_data_view::Implementation::compute_volume_of_cell(E_Int ind_cell) const
    {
        auto faces = this->get_faces_of_element(ind_cell, 0);
        return compute_volume_ngon(faces);
    }
    //# ##################################################################################################
    //_ _                                 Mise en oeuvre de la classe publique                           _
    //# ##################################################################################################

    //_ _________________________________________ Constructeur ___________________________________________
    ngon_data_view::ngon_data_view(E_Int nb_faces, const const_view_type<E_Int>& face2verts, 
                                   E_Int nb_elts,  const const_view_type<E_Int>& elt2faces,  const fields_type& fields,
                                   const coordinates_npos& pos_coords, const coordinates_npos& pos_velocity) :
        zone_data_view(new Implementation(nb_faces, face2verts, nb_elts, elt2faces, fields, 
                                          pos_coords, pos_velocity))
    {
    }
    //_ __________________________ Retourne le nombre d'éléments dans le maillage _______________________
    E_Int ngon_data_view::number_of_elements() const
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::NGON);
        auto& impl = static_cast<ngon_data_view::Implementation&>(*implementation);
        return impl.m_beg_elt2faces.size() - 1;
    }
    //_ __________________________ Retourne la connectivité face vers sommets __________________________
    auto ngon_data_view::get_face_to_vertices() const 
                                -> std::pair<const std::vector<E_Int>&, const std::vector<E_Int>&>
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::NGON);
        auto& impl = static_cast<ngon_data_view::Implementation&>(*implementation);
        return {impl.m_beg_face2verts, impl.m_face2verts};
    }
    //_ _________________________ Retourne la connectivité face vers éléments _________________________
    auto ngon_data_view::get_face_to_elements() const -> const std::vector<std::pair<E_Int,E_Int>>&
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::NGON);
        auto& impl = static_cast<ngon_data_view::Implementation&>(*implementation);
        return impl.m_face2elts;   
    }
    //_ ______________________________ Retourne la connectivité élément vers faces ___________________
    auto ngon_data_view::get_element_to_faces() const 
                                -> std::pair<const std::vector<E_Int>&,const std::vector<E_Int>&>
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::NGON);
        auto& impl = static_cast<ngon_data_view::Implementation&>(*implementation);
        return {impl.m_beg_elt2faces, impl.m_elt2faces};        
    }
}
