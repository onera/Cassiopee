/*
    Copyright 2013-2020 Onera.

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
#include <algorithm>
#include <set>
#include "vector3d.hpp"
#include "linear_algebra.hpp"
#include "triangulated_polyhedron.hpp"
#include "unstructured_data_view_p.hpp"
#include "connectivity.hpp"

namespace K_POST
{
    constexpr const std::array<unsigned char,unstructured_data_view::Implementation::number_of_element_types> 
    unstructured_data_view::Implementation::number_of_vertices_per_element;

    constexpr const std::array<unsigned char,unstructured_data_view::Implementation::number_of_element_types> 
    unstructured_data_view::Implementation::number_of_vertices_for_polyhedron;

    constexpr const std::array<unsigned char,unstructured_data_view::Implementation::number_of_element_types> 
    unstructured_data_view::Implementation::number_of_faces_per_element;

    constexpr const std::array<unsigned char,unstructured_data_view::Implementation::number_of_element_types> 
    unstructured_data_view::Implementation::total_nb_of_vertices_for_faces_per_element;

    constexpr const std::array<unsigned char,unstructured_data_view::Implementation::number_of_element_types> 
    unstructured_data_view::Implementation::number_of_triangles_per_element;

    const std::array<std::vector<unsigned char>,unstructured_data_view::Implementation::number_of_element_types> 
    unstructured_data_view::Implementation::number_of_vertices_per_face_per_element =
    { // Nombre de sommets par face pour chaque élément à l'ordre 1
        std::vector<unsigned char>{3,3,3,3},
                                  {4,3,3,3,3},
                                  {4,4,4,3,3},
                                  {4,4,4,4,4,4}
    };

    const std::array<std::vector<std::vector<unsigned char>>,unstructured_data_view::Implementation::number_of_element_types>
    unstructured_data_view::Implementation::vertices_per_face_per_element =
    {
        // Pour le tetraèdre :
        std::vector<std::vector<unsigned char>> 
        {
            std::vector<unsigned char>{0,2,1},
                                      {0,1,3},
                                      {1,2,3},
                                      {2,0,3}
        },
        // Pour la pyramide :
        {
            std::vector<unsigned char>{0,3,2,1},
                                      {0,1,4},
                                      {1,2,4},
                                      {2,3,4},
                                      {3,0,4}
        },
        // Pour le pentaèdre :
        {
            std::vector<unsigned char>{0,1,4,3},
                                      {1,2,5,4},
                                      {2,0,3,5},
                                      {0,2,1},
                                      {3,4,5}
        },
        // Pour l'hexaèdre
        {
            std::vector<unsigned char>{0,3,2,1},
                                      {0,1,5,4},
                                      {1,2,6,5},
                                      {2,3,7,6},
                                      {0,4,7,3},
                                      {4,5,6,7}
        }
    };
    // ########################################## Private part #################################### 
    void unstructured_data_view::Implementation::compute_faces_connectivity()
    {
        /**
         * L'algorithme utilisé ici est basé sur l'algorithme formel présenté dans l'article suivant :
         * 
         *             A. Logg (2008), 'Efficient Representation of Computational Meshes'
         * 
         * posté sur arXiv le quatorze mai 2012.
         * 
         */


        E_Int nb_verts = m_beg_vert2elts.size()-1;
        E_Int nb_elts  = m_elt2verts->getSize();
        E_Int nb_faces_per_elt = number_of_faces_per_element[m_type_of_element];
        //std::cout << "nb_faces_per_elt : " << nb_faces_per_elt << std::endl;
        const FldArrayI& e2v = *m_elt2verts;
        m_elt2faces.resize(nb_elts * nb_faces_per_elt);

        // Foncteur pour détecter si deux faces f1 et f2 sont identiques. 
        // On suppose le maillage orienté, donc deux faces identiques, provenant de
        // deux éléments du maillage différents présentent une orientation différente
        // (car dans les deux cas, normale sortante).
        // Supposons que pour chaque face crée ou testée, on fait en sorte que le premier sommet
        // de la face soit celui d'indice minimum. Dans ce cas, nous n'avons qu'à tester que :
        // f1[0] == f2[0], f1[1] == f2[nb_verts-1] et f1[2] == f2[nb_verts-2] !
        // on suppose aussi le maillage conforme, c'est à dire qu'une face d'un élément est commune
        // à un autre élément en entier, et pas seulement en partie !
        // Notons que le ordonnancement des sommets dans chaque face (en commençant par l'indice le plus petit),
        // n'a pas d'incidence sur le reste du code hormis la constitution géométrique des faces (qui ne demande
        // qu'à ce que les faces soient orientées vers l'extérieur)
        // puisqu'on se base principalement sur la connectivité element to vertices.
        auto same_face = [] (const E_Int nb_verts, const E_Int* indices_face1, const E_Int* indices_face2)
        {
            if ( (indices_face1[0] == indices_face2[0]) &&
                 (indices_face1[1] == indices_face2[nb_verts-1]) &&
                 (indices_face1[2] == indices_face2[nb_verts-2]) )
                return true;
            return false;
        };

        // Taille initiale prise au max...
        m_beg_face2verts.reserve(nb_elts*nb_faces_per_elt);
        m_face2verts.reserve(nb_elts * total_nb_of_vertices_for_faces_per_element[m_type_of_element]);

        m_beg_face2verts.push_back(0);
        for ( E_Int ielt = 0; ielt < nb_elts; ++ielt )
        {
            // On constitue les faces candidates de notre ième élément :
            std::vector<std::vector<E_Int>> candidates_faces(nb_faces_per_elt);
#           pragma omp parallel for
            for ( E_Int iface = 0; iface < vertices_per_face_per_element[m_type_of_element].size(); ++iface )
            {
                // Indices locaux suivant la norme CGNS des sommets constituant les faces :
                const auto& loc_faces = vertices_per_face_per_element[m_type_of_element][iface];
                E_Int nb_verts_for_face = loc_faces.size(); // Nombre de sommet dans une face
                std::vector<E_Int> face(nb_verts_for_face);
                candidates_faces[iface].reserve(nb_verts_for_face);
                for ( E_Int ivert = 0; ivert < nb_verts_for_face; ++ivert )
                    face[ivert] = e2v(ielt,loc_faces[ivert]+1)-1;// Indices globaux des sommets de la iface ème face
                // Recherche du plus petit indice dans face :
                E_Int imin = 0;
                for ( E_Int ivert = 1; ivert < nb_verts_for_face; ++ivert )
                    imin = (face[imin] < face[ivert] ? imin : ivert );
                // On commence chaque face par l'indice de sommet le plus petit en conservant l'orientation :
                for ( E_Int ivert = 0; ivert < nb_verts_for_face; ++ivert )
                    candidates_faces[iface].push_back(face[(imin+ivert)%nb_verts_for_face]);
            }
            // On va chercher pour chaque face si elle a déjà été définie par les éléments précédents
            if (ielt>0) // Sauf si c'est le premier élément :-)
            {
                // On va chercher les éléments voisins de l'élément courant. Ce sont les seuls qui peuvent avoir une face commune
                // avec notre élément.
                E_Int nb_elts_for_verts = 0;
                for (E_Int ivert = 0; ivert < number_of_vertices_per_element[m_type_of_element]; ++ivert)
                {
                    // Attention, m_elt2verts commence ses indices à un !
                    E_Int ind_vert = e2v(ielt,ivert+1)-1;
                    assert(ind_vert>=0);
                    assert(ind_vert<nb_verts);
                    nb_elts_for_verts += m_beg_vert2elts[ind_vert+1] - m_beg_vert2elts[ind_vert];
                }
                std::vector<E_Int> cells; cells.reserve(nb_elts_for_verts);
                for (E_Int ivert = 0; ivert < number_of_vertices_per_element[m_type_of_element]; ++ivert)
                {
                    E_Int ind_vert = e2v(ielt,ivert+1)-1;
                    for ( E_Int icell = m_beg_vert2elts[ind_vert]; icell < m_beg_vert2elts[ind_vert+1]; ++icell)
                        cells.push_back(m_vert2elts[icell]);
                }
                std::set<E_Int> unique_cells(cells.begin(), cells.end()); // On fait en sorte de ne pas tester plusieurs fois la même cellule voisine
                for ( size_t iface = 0; iface < candidates_faces.size(); ++iface )
                {
                    const auto& face = candidates_faces[iface]; // Pour chaque face candidate comme nouvelle face :
                    bool found_face = false;
                    // On cherche si la face de la cellule a déjà été prise en compte par une cellule voisine :
                    for ( E_Int jelt : unique_cells )
                    {
                        for (E_Int jface = 0; jface < nb_faces_per_elt; ++jface )
                        {
                            E_Int index_jface = m_elt2faces[nb_faces_per_elt*jelt + jface];
                            E_Int nb_vert_for_jface = m_beg_face2verts[index_jface+1] - m_beg_face2verts[index_jface];
                            if (nb_vert_for_jface != face.size() ) continue;// Pas même nombre de sommet, faces différentes !
                            if (same_face(nb_vert_for_jface, face.data(), m_face2verts.data() + m_beg_face2verts[index_jface]))
                            {   // On a trouvé la face déjà stockée. On met à jour l'index de la face dans m_elt2faces
                                // pour l'élément courant :
                                m_elt2faces[nb_faces_per_elt*ielt + iface] = index_jface;
                                found_face = true;
                                break; // On sort de la boucle sur les faces du jième élément
                            }
                        }
                        if (found_face) break;// OK, on a trouvé la face déjà définie, on sort de la boucle sur les éléments précédents
                    }
                    if (found_face == false) // Non, on n'a pas trouvé cette face déjà définie
                    {
                        // On la rajoute :
                        m_elt2faces[nb_faces_per_elt*ielt + iface] = m_beg_face2verts.size()-1;
                        m_beg_face2verts.push_back(m_beg_face2verts.back()+face.size());
                        m_face2verts.insert(m_face2verts.end(), face.begin(),face.end());
                    }
                }// Fin for (E_Int iface ... ) : On passe à la face suivante de l'élément en train d'être traité
            }// Fin if (ielt> 0)
            else
            {
                // C'est le premier élément visité, donc toutes les faces sont nouvelles !
                for( E_Int iface = 0; iface < candidates_faces.size(); ++iface )
                {
                    const auto& face = candidates_faces[iface]; // Pour chaque face candidate comme nouvelle face :
                    m_elt2faces[iface] = m_beg_face2verts.size()-1;
                    m_beg_face2verts.push_back(m_beg_face2verts.back()+face.size());
                    m_face2verts.insert(m_face2verts.end(), face.begin(),face.end());
                }
            }
        } // For (E_int ielt ... ) : on passe à l'élément suivant pour traiter ses faces
        m_elt2faces.shrink_to_fit();
        m_beg_face2verts.shrink_to_fit();
        m_face2verts.shrink_to_fit();
        /*std::cout << "elt2verts : " << std::endl;
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt)
        {
            std::cout << "\telt n°" << ielt << " : ";
            for ( E_Int iv = 0; iv < m_elt2verts->getNfld(); ++iv )
                std::cout << (*m_elt2verts)(ielt,iv+1)-1 << " ";
            std::cout << std::endl;
        }
        std::cout << "elt2faces : " << std::endl;
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt)
        {
            std::cout << "elt n°" << ielt << " : ";
            for ( E_Int iv = 0; iv < number_of_vertices_per_element[m_type_of_element]; ++iv)
            {
                std::cout << m_elt2faces[ielt*number_of_vertices_per_element[m_type_of_element]+iv] << " ";
            }
            std::cout << std::endl;
        }
        for ( E_Int ef : m_elt2faces)
            std::cout << ef << " ";
        std::cout << "\nface2verts : " << std::endl;
        for ( E_Int pt_f  = 0; pt_f < m_beg_face2verts.size()-1; ++pt_f )
        {
            std::cout << "\tface n°" << pt_f << " : ";
            for ( E_Int iface = m_beg_face2verts[pt_f]; iface < m_beg_face2verts[pt_f+1]; ++iface )
            {
                std::cout << m_face2verts[iface] << " ";
            }
            std::cout << std::endl;
        }*/
        // On calcule m_face2elts en sachant que first pour element correspondant à face directe
        // et second à élément avec face indirecte
        // Au vue de la construction ci--dessus, le premier élément faisant référence à une face donnée aura forcément
        // la normale à la face sortante (ce que l'on veut), tandis que le second élément y faisant référence
        // aura de son point de vue la normale à la face rentrante. Donc en parcourant les éléments, il suffit de remplir
        // la première fois pour une face le champs first avec le n° de l'élément (qui aura la face sortante) et la seconde
        // fois on remplit le champs second avec le n° de l'élément (qui aura la face rentrante)
        E_Int nb_faces = m_beg_face2verts.size()-1;
        std::vector<std::pair<E_Int,E_Int>>(nb_faces,{-1,-1}).swap(m_face2elts);
        for (E_Int ielt = 0; ielt < nb_elts; ++ielt)
        {
            for (E_Int ieltface = 0; ieltface < nb_faces_per_elt; ++ieltface )
            {
                E_Int index_face = m_elt2faces[ielt*nb_faces_per_elt+ieltface];
                if (m_face2elts[index_face].first == -1)
                    m_face2elts[index_face].first = ielt;
                else
                {
                    assert(m_face2elts[index_face].second == -1);
                    m_face2elts[index_face].second = ielt;
                }
            }
        }
        /*std::cout << "face 2elts : ";
        for ( const auto& fe : m_face2elts )
        {
            std::cout << "{" << fe.first << ", " << fe.second << "}";
        }
        std::cout << std::endl;*/
    }
    // ########################################### Public part ####################################
    unstructured_data_view::Implementation::Implementation( const char* str_elt_type, const pt_connect_type elt2verts,
                                                            const fields_type& fields, const coordinates_npos& pos_coords, 
                                                            const coordinates_npos& pos_velocity) : 
        zone_data_view::Implementation( zone_data_view::ELEMENT_TYPE, pos_coords, fields, pos_velocity),
        m_type_of_element(str_type_of_element_to_element_type(str_elt_type)),
        m_elt2verts(elt2verts)
    {
        //std::cout << "Type d'élément détecté : " << m_type_of_element << std::endl;
        //std::cout << "nb faces par élément : " << int(number_of_faces_per_element[m_type_of_element]) << std::flush <<  std::endl;
        std::tie(m_beg_vert2elts, m_vert2elts) = compute_vertex_to_elements( fields->getSize(), *elt2verts); 
        this->compute_faces_connectivity();
    }
    // ---------------------------------------------------------------------------------------------------
    std::vector<face> unstructured_data_view::Implementation::get_faces_of_element( E_Int number, E_Int no_zone ) const
    {
        E_Int nb_faces_per_elt = number_of_faces_per_element[m_type_of_element];
        std::vector<face> faces; faces.reserve(nb_faces_per_elt);
        for (E_Int iface = 0; iface < nb_faces_per_elt; ++iface)
        {
            E_Int index_face = m_elt2faces[nb_faces_per_elt*number+iface];
            E_Int beg_vertices = m_beg_face2verts[index_face  ];
            E_Int end_vertices = m_beg_face2verts[index_face+1];
            const auto& f2e = m_face2elts[index_face];
            bool is_direct = (f2e.first == number ? true : false);
            std::vector<E_Int> coords;
            if (is_direct) // Ok, face stockée avec normale sortante pour cet élément :
                std::vector<E_Int>(m_face2verts.begin()+beg_vertices, m_face2verts.begin()+end_vertices).swap(coords);
            else
            {   // Ah non, face stockée avec normale rentrante pour cet élément. Donc on va
                // recréer la face en parcourant ses sommets à l'envers :
                coords.reserve(end_vertices-beg_vertices);
                for (E_Int ivert = end_vertices-1; ivert >= beg_vertices; --ivert)
                    coords.push_back(m_face2verts[ivert]);
            }
            E_Int op_elt = ( is_direct ? f2e.second : f2e.first );
            faces.emplace_back(coords, std::pair<E_Int,E_Int>{number,op_elt}, this->getCoordinates());
        }
        return faces;
    }
    // ---------------------------------------------------------------------------------------------------
    bool 
    unstructured_data_view::Implementation::is_containing(E_Int ind_elt, const point3d& pt) const
    {
        const auto& crds = this->getCoordinates();
        std::array<std::vector<double>,3> coords; // Nombre de sommet dépend du type d'élément (avec barycentre pour certains)
        coords[0].reserve(number_of_vertices_for_polyhedron[m_type_of_element]);
        coords[1].reserve(number_of_vertices_for_polyhedron[m_type_of_element]);
        coords[2].reserve(number_of_vertices_for_polyhedron[m_type_of_element]);
        E_Int nb_verts_per_elt = number_of_vertices_per_element[m_type_of_element];
        // On extrait tous les points de la cellule selon l'ordre CGNS :
        for (E_Int ivert = 0; ivert < nb_verts_per_elt; ++ivert)
        {
            E_Int ind_vert = (*m_elt2verts)(ind_elt,ivert+1)-1;
            coords[0].push_back(crds[0][ind_vert]);
            coords[1].push_back(crds[1][ind_vert]);
            coords[2].push_back(crds[2][ind_vert]);
        }
        // On doit utiliser des triangles en tessalisant les faces de l'élément :
        using triangle_type = triangulated_polyhedron::triangle_type;
        std::vector<triangle_type> faces; faces.reserve(number_of_triangles_per_element[m_type_of_element]);
        // On parcourt les faces dans l'ordre CGNS :
        for ( E_Int iface = 0; iface < number_of_faces_per_element[m_type_of_element]; ++iface)
        {
            if (number_of_vertices_per_face_per_element[m_type_of_element][iface] == 3)
            {  // C'est un triangle, pas besoin de rajouter un barycentre :
                assert(vertices_per_face_per_element[m_type_of_element][iface].size() == 3);
                E_Int ind1 = vertices_per_face_per_element[m_type_of_element][iface][0];
                E_Int ind2 = vertices_per_face_per_element[m_type_of_element][iface][1];
                E_Int ind3 = vertices_per_face_per_element[m_type_of_element][iface][2];
                faces.emplace_back(triangle_type{ind1, ind2, ind3}); // On rajoute directement la facette
            }
            else
            {   // C'est un quadrangle. On va devoir rajouter le barycentre pour
                // tessaliser sans être ambigu entre éléments.
                assert(number_of_vertices_per_face_per_element[m_type_of_element][iface] == 4);
                E_Int ind1 = vertices_per_face_per_element[m_type_of_element][iface][0];
                E_Int ind2 = vertices_per_face_per_element[m_type_of_element][iface][1];
                E_Int ind3 = vertices_per_face_per_element[m_type_of_element][iface][2];
                E_Int ind4 = vertices_per_face_per_element[m_type_of_element][iface][3];
                E_Int ind5 = coords[0].size(); // Index du point barycentre mis à la fin de coords
                // Calcul d'un barycentre et rajout d'un point à coords :
                coords[0].push_back(0.25*(coords[0][ind1]+coords[0][ind2]+coords[0][ind3]+coords[0][ind4]));
                coords[1].push_back(0.25*(coords[1][ind1]+coords[1][ind2]+coords[1][ind3]+coords[1][ind4]));
                coords[2].push_back(0.25*(coords[2][ind1]+coords[2][ind2]+coords[2][ind3]+coords[2][ind4]));
                faces.emplace_back(triangle_type{ind1, ind2, ind5});
                faces.emplace_back(triangle_type{ind2, ind3, ind5});
                faces.emplace_back(triangle_type{ind3, ind4, ind5});
                faces.emplace_back(triangle_type{ind4, ind1, ind5});
            }
        }
        // Construction du polyèdre :
        /*std::cout << "point à localiser : " << std::string(pt) << std::endl;
        std::cout << "polyèdre généré : " << std::endl;
        for ( const auto& face : faces )
        {
            std::cout << "{" << std::string(point3d{coords[0][face[0]],coords[1][face[0]],coords[2][face[0]]})
                      << "," << std::string(point3d{coords[0][face[1]],coords[1][face[1]],coords[2][face[1]]})
                      << "," << std::string(point3d{coords[0][face[2]],coords[1][face[2]],coords[2][face[2]]})
                      << "}" << std::endl;
        }*/

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
            std::cerr << "Warning : interpolated point is on interface. Possibility to have two points in same location in the stream line"
                      << std::flush << std::endl;
            is_inside = true; // Dans ce cas, on considère qu'on est à l'intérieur (on prend l'élément comme un fermé topologique)
        }
        //std::cout << "is inside " << std::boolalpha << is_inside << std::endl;
        return is_inside;
    }
    // ---------------------------------------------------------------------------------------------------
    E_Int
    unstructured_data_view::Implementation::get_interpolation_cell( const point3d& point ) const
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
            for ( E_Int ivert = 0; ivert <  number_of_vertices_per_element[m_type_of_element]; ++ivert )
            {
                E_Int index_vert = (*this->m_elt2verts)(elt,ivert+1)-1;
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
    // ---------------------------------------------------------------------------------------------------
    void 
    unstructured_data_view::Implementation::compute_interpolated_field( const point3d& pt, E_Int ind_cell, 
                                                                        E_Int ipos, FldArrayF& interpolatedField ) const
    {
        E_Int nb_verts_per_elt = number_of_vertices_per_element[m_type_of_element];
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
            E_Int ind_vert = (*m_elt2verts)(ind_cell,ivert+1)-1;
            coords.emplace_back(crds[0][ind_vert], crds[1][ind_vert], crds[2][ind_vert]);
            for (E_Int ifld = 0; ifld < nfld; ++ifld)
            {
                values[ifld].push_back((*fld)(ind_vert,ifld+1));
            }
        }
        // L'interpolation va dépendre ici du type d'élément :
        if (this->m_type_of_element == tetraedre)
        {
            // L'interpolation la plus simple :
            // Soit le tetraèdre T(p₀, p₁, p₂, p₃)    →  →  →
            // On choisit comme repère barycentrique (e₁,e₂,e₃)
            //       →                       →                        →
            // avec  e₁ le vecteur (p₀, p₁), e₂ le vecteur (p₀,p₂) et e₃ le vecteur (p₀,p₃)
            // 
            // On calcule les coordonnées barycentriques (𝛼,𝛽,𝛾) du point p où on doit interpoler : :
            //        →      →      →
            // p₀ + 𝛼.e₁ + 𝛽.e₂ + 𝛾.e₃ = p avec 𝛼 + 𝛽 + 𝛾 ≤ 1, 𝛼 ≥ 0, 𝛽 ≥ 0, 𝛾 ≥ 0
            // 
            // ce qui revient à inverser un système linéaire de dimension trois.
            // Pour calculer le champs interpolé, il suffit alors de calculer :
            // 
            // f(p) = (1-𝛼-𝛽-𝛾).f(p₀) + 𝛼.f(p₁) + 𝛽.f(p₂) + 𝛾.f(p₃)
            // 
            vector3d e1(coords[0],coords[1]);
            vector3d e2(coords[0],coords[2]);
            vector3d e3(coords[0],coords[3]);
            vector3d pp0(coords[0], pt);

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
                interpolatedField(ipos,ifld+1) = umabg * values[ifld][0] + alpha * values[ifld][1] + beta * values[ifld][2] +
                                                 gamma * values[ifld][3];
            }
        }
        if (this->m_type_of_element == pyramide)
        {
            // Soit la pyramide P(p₀, p₁, p₂, p₃, p₄) →  →  →
            // On choisit comme repère barycentrique (e₁,e₂,e₃)
            //      →                      →                        →
            // avec e₁ le vecteur (p₀,p₁), e₂ le vecteur (p₀,p₃) et e₃ le vecteur (p₀, p₄)
            // ( en suivant la norme CGNS de numérotation des sommets )
            // 
            // On va chercher les polynômes d'interpolation de chaque champs à l'aide de polynômes de la forme :
            // 
            // f(𝛼,𝛽,𝛾) = a₀ + a₁.𝛼 + a₂.𝛽 + a₃.𝛾 + a₄.𝛼𝛽
            // 
            // On calcule les coordonnées barycentriques (𝛼₂,𝛽₂,𝛾₂) du point p₂ tel que :
            //         →       →       →
            // p₀ + 𝛼₂.e₁ + 𝛽₂.e₂ + 𝛾₂.e₃ = p₂
            // 
            // en résolvant un système linéaire de dimension trois ainsi que les coordonnées barycentriques (𝛼ₚ, 𝛽ₚ, 𝛾ₚ)
            //  du point p où on interpole :
            //         →       →      →
            // p₀ + 𝛼ₚ.e₁ + 𝛽ₚ.e₂ + 𝛾ₚ.e₃ = p
            // 
            // Puisqu'on connaît les valeurs du champs aux sommets de l'élément, c'est à dire qu'on connait
            // en particuliers les valeurs de f(0,0,0), f(1,0,0), f(0,1,0) et f(0,0,1), on en déduit que :
            // 
            // a₀ = f(0,0,0)       → valeur du champs au point p₀
            // a₁ = f(1,0,0) - a₀  → valeur du champs au point p₁
            // a₂ = f(0,1,0) - a₀  → valeur du champs au point p₃
            // a₃ = f(0,0,1) - a₀  → valeur du champs au point p₄
            // 
            // On connaît également f(𝛼₂,𝛽₂,𝛾₂), la valeur du champs au point p₂, ce qui nous donne a₄ :
            // 
            // a₄ = (f(𝛼₂,𝛽₂,𝛾₂) - a₀ - 𝛼₂.a₁ - 𝛽₂.a₂ - 𝛾₂.a₃)/(𝛼₂.𝛽₂)
            // 
            // Il s'en suit que la valeur du champs au point p vaut :
            // 
            // f(𝛼ₚ, 𝛽ₚ, 𝛾ₚ) = a₀ + a₁.𝛼ₚ + a₂.𝛽ₚ + a₃.𝛾ₚ + a₄.𝛼ₚ𝛽ₚ
            // 
            vector3d e1(coords[0],coords[1]);
            vector3d e2(coords[0],coords[3]);
            vector3d e3(coords[0],coords[4]);
            vector3d p0p2(coords[0], coords[2]);
            vector3d p0p(coords[0], pt);

            matrix_3x3_type A{std::array<double,3>{e1.x,e2.x,e3.x},
                                                  {e1.y,e2.y,e3.y},
                                                  {e1.z,e2.z,e3.z}
                             };
            matrix_3x3_type invA = inverse(A);
            auto p2_barycrds = invA*p0p2;
            double alpha_2 = p2_barycrds[0];
            double beta_2  = p2_barycrds[1];
            double gamma_2 = p2_barycrds[2];
            auto p_barycrds  = invA*p0p;
            double alpha_p = p_barycrds[0];
            double beta_p  = p_barycrds[1];
            double gamma_p = p_barycrds[2];
            for ( E_Int ifld = 0; ifld < nfld; ++ifld)
            {
                double a0 = values[ifld][0];
                double a1 = values[ifld][1] - a0;
                double a2 = values[ifld][3] - a0;
                double a3 = values[ifld][4] - a0;                
                double a4 = (values[ifld][2] - a0 - alpha_2*a1 -beta_2*a2 - gamma_2*a3)/(alpha_2*beta_2);
                interpolatedField(ipos,ifld+1) = a0 + a1 * alpha_p + a2 * beta_p + a3 * gamma_p + a4 * alpha_p * beta_p;
            }
        }
        if (this->m_type_of_element == pentaedre)
        {
            // Soit le pentaedre P(p₀, p₁, p₂, p₃, p₄, p₅)
            //                                     →  →  →
            // On choisit le repère barycentrique (e₁,e₂,e₃) avec
            //      →                      →                        →
            // avec e₁ le vecteur (p₀,p₁), e₂ le vecteur (p₀,p₂) et e₃ le vecteur (p₀, p₃)
            // ( en suivant la norme CGNS de numérotation des sommets )
            // 
            // On va chercher les polynômes d'interpolation de chaque champs à l'aide de polynômes de la forme :
            // 
            // f(𝛼,𝛽,𝛾) = a₀ + a₁.𝛼 + a₂.𝛽 + a₃.𝛾 + a₄.𝛼𝛾 + a₅.𝛽𝛾
            // 
            // On calcule les coordonnées barycentriques (𝛼ᵢ,𝛽ᵢ,𝛾ᵢ) (i=4 ou 5) du point pᵢ tel que :
            //         →       →       →
            // p₀ + 𝛼ᵢ.e₁ + 𝛽ᵢ.e₂ + 𝛾ᵢ.e₃ = pᵢ
            // 
            // en résolvant un système linéaire de dimension trois ainsi que les coordonnées barycentriques (𝛼ₚ, 𝛽ₚ, 𝛾ₚ)
            //  du point p où on interpole :
            //         →       →      →
            // p₀ + 𝛼ₚ.e₁ + 𝛽ₚ.e₂ + 𝛾ₚ.e₃ = p
            // 
            // Puisqu'on connaît les valeurs du champs aux sommets de l'élément, c'est à dire qu'on connait
            // en particuliers les valeurs de f(0,0,0), f(1,0,0), f(0,1,0) et f(0,0,1), on en déduit que :
            // 
            // a₀ = f(0,0,0)       → valeur du champs au point p₀
            // a₁ = f(1,0,0) - a₀  → valeur du champs au point p₁
            // a₂ = f(0,1,0) - a₀  → valeur du champs au point p₂
            // a₃ = f(0,0,1) - a₀  → valeur du champs au point p₃
            //
            // Connaissant également la valeur du champs aux points p₄ et p₅, on peut déduire les valeurs de a₄ et a₅
            // en résolvant le système linéaire suivant :
            // 
            // ⎛ 𝛼₄𝛾₄  𝛽₄𝛾₄⎞ ⎛ a₄⎞   ⎛ f(𝛼₄,𝛽₄,𝛾₄) - a₀ - a₁.𝛼₄ - a₂.𝛽₄ - a₃.𝛾₄⎞
            // ⎝ 𝛼₅𝛾₅  𝛽₅𝛾₅⎠ ⎝ a₅⎠ = ⎝ f(𝛼₅,𝛽₅,𝛾₅) - a₀ - a₁.𝛼₅ - a₂.𝛽₅ - a₃.𝛾₅⎠
            // 
            // Il ne reste plus qu'à interpoler le champs au point p :
            // 
            // f(𝛼ₚ, 𝛽ₚ, 𝛾ₚ) = a₀ + a₁.𝛼ₚ + a₂.𝛽ₚ + a₃.𝛾ₚ + a₄.𝛼ₚ𝛾ₚ + a₅.𝛽ₚ𝛾ₚ
            // 
            vector3d e1(coords[0],coords[1]);
            vector3d e2(coords[0],coords[2]);
            vector3d e3(coords[0],coords[3]);
            vector3d p0p4(coords[0], coords[4]);
            vector3d p0p5(coords[0], coords[5]);
            vector3d p0p(coords[0], pt);

            matrix_3x3_type A{std::array<double,3>{e1.x,e2.x,e3.x},
                                                  {e1.y,e2.y,e3.y},
                                                  {e1.z,e2.z,e3.z}
                             };
            matrix_3x3_type invA = inverse(A);
            auto p4_barycrds = invA*p0p4;
            double alpha_4 = p4_barycrds[0];
            double beta_4  = p4_barycrds[1];
            double gamma_4 = p4_barycrds[2];
            auto p5_barycrds = invA*p0p5;
            double alpha_5 = p5_barycrds[0];
            double beta_5  = p5_barycrds[1];
            double gamma_5 = p5_barycrds[2];
            auto p_barycrds  = invA*p0p;
            double alpha_p = p_barycrds[0];
            double beta_p  = p_barycrds[1];
            double gamma_p = p_barycrds[2];            
            for ( E_Int ifld = 0; ifld < nfld; ++ifld)
            {
                double a0 = values[ifld][0];
                double a1 = values[ifld][1] - a0;
                double a2 = values[ifld][2] - a0;
                double a3 = values[ifld][3] - a0;
                matrix_2x2_type A2{ std::array<double,2>{alpha_4*gamma_4, beta_4*gamma_4},
                                                        {alpha_5*gamma_5, beta_5*gamma_5}
                                  };
                auto invA2 = inverse(A2);
                vector2d b2{ values[ifld][4] - a0 - a1*alpha_4 - a2*beta_4 - a3*gamma_4,
                             values[ifld][5] - a0 - a1*alpha_5 - a2*beta_5 - a3*gamma_5 };
                vector2d res = invA2 * b2;
                double a4 = res[0];
                double a5 = res[1];
                interpolatedField(ipos,ifld+1) = a0 + a1 * alpha_p + a2 * beta_p + a3 * gamma_p + a4 * alpha_p * gamma_p +
                                                 a5 * beta_p * gamma_p;
            }
        }
        if (this->m_type_of_element == hexaedre)
        {
            // Soit un hexaèdre représenté par huit sommets : H : {p₀, p₁, p₂, p₃, p₄, p₅, p₆, p₇}
            //         →                       →                        →
            // On note e₁ le vecteur (p₀, p₁), e₂ le vecteur (p₀,p₃) et e₃ le vecteur (p₀,p₄)
            // ( en suivant la norme CGNS de numérotation des sommets )
            // 
            // On va utiliser une interpolation trilinéaire en
            // prenant pour polynôme d'interpolation f({𝛼,𝛽,𝛾}) = a₀ + a₁.𝛼 + a₂.𝛽 + a₃.𝛾 + a₄.𝛼𝛽 + a₅.𝛼𝛾 + a₆.𝛽𝛾 + a₇.𝛼𝛽𝛾
            // avec {𝛼,𝛽,𝛾} ∈ [0;1]³           →      →      →
            // Il est clair que :  p₀ = p₀ + 0.e₁ + 0.e₂ + 0.e₃ soit {𝛼,𝛽,𝛾} = {0,0,0}
            //                                 →      →      →
            //                     p₁ = p₀ + 1.e₁ + 0.e₂ + 0.e₃ soit {𝛼,𝛽,𝛾} = {1,0,0}
            //                                 →      →      →
            //                     p₃ = p₀ + 0.e₁ + 1.e₂ + 0.e₃ soit {𝛼,𝛽,𝛾} = {0,1,0}
            //                                 →      →      →
            //                     p₄ = p₀ + 0.e₁ + 0.e₂ + 1.e₃ soit {𝛼,𝛽,𝛾} = {0,0,1}
            //                                                              →   →     →
            // Pour trouver les autres sommets de l'hexahèdre par rapport à e₁, e₂ et e₃, il faut résoudre le système
            // linéaire suivant (pour le sommet pᵢ, i ∈ {2,5,6,7}) :
            //   →    →    →   ⎛𝛼⎞
            //  (e₁ | e₂ | e₃) ⎜𝛽⎟ = pᵢ - p₀
            //                 ⎝𝛾⎠                                                                        →   →   →
            // On obtient alors un triplet (𝛼ᵢ,𝛽ᵢ,𝛾ᵢ) permettant de représenter le sommet dans le repère e₁, e₂, e₃
            // Puisqu'on connaît la valeur du champs aux sommets de l'élément, c'est à dire que
            //  f(pᵢ), i ∈ {1,2,3,4,5,6,7}  sont connus.
            //  D'après le choix du polynôme d'interpolation qu'on a fait plus haut :
            //      f(p₀) = a₀
            //      f(p₁) = a₀ + a₁ ⇒ a₁ = f(p₁) - a₀
            //      f(p₃) = a₀ + a₂ ⇒ a₂ = f(p₃) - a₀
            //      f(p₄) = a₀ + a₃ ⇒ a₃ = f(p₄) - a₀
            //  et pour les autres sommets, on aura donc :
            //      f(pᵢ) = a₀ + a₁.𝛼ᵢ + a₂.𝛽ᵢ + a₃.𝛾ᵢ + a₄.𝛼ᵢ𝛽ᵢ + a₅.𝛼ᵢ𝛾ᵢ + a₆.𝛽ᵢ𝛾ᵢ + a₇.𝛼ᵢ𝛽ᵢ𝛾ᵢ
            //  ce qui nous donne un système linéaire de dimension quatre à résoudre :
            //      a₄.𝛼ᵢ𝛽ᵢ + a₅.𝛼ᵢ𝛾ᵢ + a₆.𝛽ᵢ𝛾ᵢ + a₇.𝛼ᵢ𝛽ᵢ𝛾ᵢ = f(pᵢ) - a₀ - a₁.𝛼ᵢ - a₂.𝛽ᵢ - a₃.𝛾ᵢ
            //   avec i ∈ {2,5,6,7}.
            //   Il faudra ensuite déterminer les coordonnées (𝛼ₚ,𝛽ₚ,𝛾ₚ) de notre point à interpoler, puis appliquer la fonction
            //   polynomiale ainsi calculée.
            //   
            vector3d e1(coords[0], coords[1]);
            vector3d e2(coords[0], coords[3]);
            vector3d e3(coords[0], coords[4]);
            matrix_3x3_type A{std::array<double,3>{ e1.x ,e2.x, e3.x},
                                                  { e1.y ,e2.y, e3.y},
                                                  { e1.z ,e2.z, e3.z}
                             };
            matrix_3x3_type invA = inverse(A);
            vector3d bar2 = invA*vector3d(coords[0],coords[2]);
            vector3d bar5 = invA*vector3d(coords[0],coords[5]);
            vector3d bar6 = invA*vector3d(coords[0],coords[6]);
            vector3d bar7 = invA*vector3d(coords[0],coords[7]);

            vector3d bary_pt= invA*vector3d(coords[0], pt);

            //      a₄.𝛼ᵢ𝛽ᵢ + a₅.𝛼ᵢ𝛾ᵢ + a₆.𝛽ᵢ𝛾ᵢ + a₇.𝛼ᵢ𝛽ᵢ𝛾ᵢ = f(pᵢ) - a₀ - a₁.𝛼ᵢ - a₂.𝛽ᵢ - a₃.𝛾ᵢ
            matrix_4x4_type B{ std::array<double,4>{bar2.x*bar2.y, bar2.x*bar2.z, bar2.y*bar2.z, bar2.x*bar2.y*bar2.z},
                                                   {bar5.x*bar5.y, bar5.x*bar5.z, bar5.y*bar5.z, bar5.x*bar5.y*bar5.z},
                                                   {bar6.x*bar6.y, bar6.x*bar6.z, bar6.y*bar6.z, bar6.x*bar6.y*bar6.z},
                                                   {bar7.x*bar7.y, bar7.x*bar7.z, bar7.y*bar7.z, bar7.x*bar7.y*bar7.z}
                             };
            auto LUB = factorize(B);
            for ( E_Int ifld = 0; ifld < nfld; ++ifld)
            {
                double a0 = values[ifld][0];        // a₀ = f(p₀)
                double a1 = values[ifld][1] - a0;   // a₁ = f(p₁) - a₀
                double a2 = values[ifld][3] - a0;   // a₂ = f(p₃) - a₀
                double a3 = values[ifld][4] - a0;   // a₃ = f(p₄) - a₀
                vector4d b{
                    values[ifld][2] - a0 -a1*bar2.x -a2*bar2.y - a3*bar2.z, 
                    values[ifld][5] - a0 -a1*bar5.x -a2*bar5.y - a3*bar5.z, 
                    values[ifld][6] - a0 -a1*bar6.x -a2*bar6.y - a3*bar6.z, 
                    values[ifld][7] - a0 -a1*bar7.x -a2*bar7.y - a3*bar7.z
                          };
                auto x = inverse_linear_system(LUB, b);
                double a4 = x[0];
                double a5 = x[1];
                double a6 = x[2];
                double a7 = x[3];
                interpolatedField(ipos,ifld+1) = a0 + a1*bary_pt.x + a2*bary_pt.y + a3*bary_pt.z +
                                                 a4*bary_pt.x*bary_pt.y + a5*bary_pt.x*bary_pt.z + a6*bary_pt.y*bary_pt.z +
                                                 a7*bary_pt.x*bary_pt.y*bary_pt.z;
            }
        }
    }
    // ===================================================================================================
    unstructured_data_view::unstructured_data_view(const char* eltType, const pt_connect_type elt2verts,
                                                   const fields_type& fields,
                                                   const coordinates_npos& pos_coords, const coordinates_npos& pos_velocity) :
        zone_data_view(new Implementation(eltType, elt2verts, fields, pos_coords, pos_velocity))
    {
    }
    // ===================================================================================================
    E_Int unstructured_data_view::number_of_elements() const
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::ELEMENT_TYPE);
        auto& impl = static_cast<unstructured_data_view::Implementation&>(*implementation);
        assert(impl.m_elt2faces.size() % impl.number_of_faces_per_element[impl.m_type_of_element] == 0);
        return impl.m_elt2faces.size() / impl.number_of_faces_per_element[impl.m_type_of_element];
    }
    // ---------------------------------------------------------------------------------------------------
    auto unstructured_data_view::get_face_to_vertices() const 
                                -> std::pair<const std::vector<E_Int>&, const std::vector<E_Int>&>
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::ELEMENT_TYPE);
        auto& impl = static_cast<unstructured_data_view::Implementation&>(*implementation);
        return {impl.m_beg_face2verts, impl.m_face2verts};
    }
    // ---------------------------------------------------------------------------------------------------
    auto unstructured_data_view::get_face_to_elements() const 
                                                        -> const std::vector<std::pair<E_Int,E_Int>>&
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::ELEMENT_TYPE);
        auto& impl = static_cast<unstructured_data_view::Implementation&>(*implementation);
        return impl.m_face2elts;        
    }
    // ---------------------------------------------------------------------------------------------------
    auto unstructured_data_view::get_element_to_faces() const -> const std::vector<E_Int>&
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::ELEMENT_TYPE);
        auto& impl = static_cast<unstructured_data_view::Implementation&>(*implementation);
        return impl.m_elt2faces;        
    }
}