#include "ribbon_stream_line.hpp"
#include "point3d.hpp"
#include "vector3d.hpp"
#include "volume.hpp"
using K_POST::point3d;
using K_POST::vector3d;
//#define DEBUG_VERBOSE
//# ####################################################################
//# # Définition de la structure privée d'implémentation de la classe  #
//# ####################################################################
struct K_POST::ribbon_streamline::Implementation {
    FldArrayF fields;
    FldArrayI connectivity;
    Implementation( const point3d &init_pos, 
                    const std::vector<zone_data_view>& zones,
                    E_Int max_vertices_for_streamline, 
                    E_Float width, bool is_bidirectional );
};


//# ###########################################################################
//# #                  Fonctions de mise en oeuvre privées                    #
//# ###########################################################################
namespace
{
    //_ _____________________ Initialisation de la ribbon stream ______________
    //@brief     Retourne n° de zone, n° de cellule, champ interpolé du pt init

    //@param[in]  init_pos                     La position du point initial  de
    //-                                         la streamline                  
    //@param[in]  zones                        Les différentes zones           
    //@param[in]  max_vertices_for_streamline  Le  nombre  maximal  de  sommets
    //-                                        pour la streamline              
    //@param[in]  width                        Largeur du ruban                
    //@param      streamPt                     Le champs porté par  streamline 
    //@param[in]  istream                      La position dans le champs de la
    //-                                        valeur initial (utile  quand  on
    //-                                        change de zone)                 
    //@param[in]  num_zone                     Numéro de la zone avt changement
    //-                                        de zone (-1 : pas de zone)      
    //@return     Triplet : n° de zone, n° de cellule et vitesse interpolée  où
    //-           se trouve le point init_pos et le vecteur porteur du ruban   
    std::tuple<E_Int, E_Int, vector3d, vector3d>
    init_ribbon_streamline( const point3d &init_pos, 
                            const std::vector<K_POST::zone_data_view> &zones,
                            E_Int max_vertices_for_streamline, E_Float width, 
                            FldArrayF &streamPt, int istream, E_Int num_zone, 
                            bool is_with_perturbation = false )
    {
#if defined( DEBUG_VERBOSE )
        std::cout << "Nbre de zones dans lesquels chercher : " << zones.size() << std::endl;
#endif
        // Recherche du n° de zone et de l'indice élément où se trouve le sommet 
        // initial de la ribbon stream :
        decltype( zones.size() ) izone = 0;
        E_Int icell = -1;
        for ( const auto &zo : zones ) {
            if ( E_Int( izone ) != num_zone ) {
#if defined( DEBUG_VERBOSE )
                std::cout << "Recherche dans la zone n°" << izone << std::endl;
#endif
                icell = zo.get_interpolation_cell( init_pos );
                if ( icell != -1 ) break;
            }
            ++izone;
        }
        if ( izone == zones.size() ) {
            if ( !is_with_perturbation ) {
                // On n'a pas trouvé de domaine adéquat. Le point est sans doute au bord d'un bloc extérieur.
                // On va perturbé le point initial pour trouvé une cellule :
                for ( vector3d perturb : std::array<vector3d, 6>{ vector3d{ 1.E-6, 0., 0. }, vector3d{ -1.E-6, 0., 0. },
                                                                  vector3d{ 0., 1.E-6, 0. }, vector3d{ 0., -1.E-6, 0. },
                                                                  vector3d{ 0., 0., 1.E-6 }, vector3d{ 0., 0., -1.E-6 } } ) {
                    auto pos2 = init_pos + perturb;
                    auto res2 = init_ribbon_streamline( pos2, zones, max_vertices_for_streamline, width, streamPt, istream, num_zone, true );
                    if ( std::get<0>( res2 ) != -1 ) return res2;
                }
            }
            //throw std::domain_error("Wrong starting point : no zones contain this point.");
            return std::tuple<E_Int, E_Int, vector3d, vector3d>( -1, -1, { 0., 0., 0. }, {0.,0.,0.} );
        }
        //if (izone > zones.size()) throw std::domain_error("Wrong starting point : no zones contain this point.");

#if defined( DEBUG_VERBOSE )
        std::cout << "bloc ou se trouve le point initial = " << izone << std::endl;
        std::cout << "Cellule contenant le point initial : " << icell << std::flush << std::endl;
#endif
        E_Int num_blk = izone;
        const K_POST::zone_data_view &zone = zones[ izone ];

        // Interpolation du champs de la zone sur le sommet de départ de la streamline :
        E_Int nfld = zone.get_fields()->getNfld();
        // Deux fois dans l'allocation car on a un ribbon stream !!
        if ( istream == 0 ) streamPt.malloc( 2*max_vertices_for_streamline, nfld );
        zone.compute_interpolated_field( init_pos, icell, istream, streamPt );
        for ( int ifld = 0; ifld < streamPt.getNfld(); ++ifld )
        {
            //std::cout << "interpolated field : " << streamPt(istream, ifld+1) << " ";
            streamPt(istream+1,ifld+1) = streamPt(istream  ,ifld+1); 
        }
        auto npos_vel = zone.get_position_of_velocity();
        vector3d velocity_pt{ streamPt( istream, npos_vel[ 0 ] ), streamPt( istream, npos_vel[ 1 ] ), streamPt( istream, npos_vel[ 2 ] ) };

        // Initialise les coordonnées du premier sommet de la streamline aux coordonnées du point initial :
        auto pos_crds = zones[ num_blk ].get_position_of_coordinates();
        streamPt( istream, pos_crds[ 0 ] ) = init_pos.x;
        streamPt( istream, pos_crds[ 1 ] ) = init_pos.y;
        streamPt( istream, pos_crds[ 2 ] ) = init_pos.z;
        // On rajoute un point orthogonal à la vitesse (pris arbitrairement dans le plan orthogonal au vecteur vitesse)
        // si c'est les deux premiers points du ribbon :
#if defined( DEBUG_VERBOSE )
        std::cout << "Calcul du vecteur orthogonal initial pour la ribbon" << std::endl;
#endif
        vector3d dir_ortho;
        if (istream == 0)
        {
            dir_ortho = vector3d{velocity_pt[1],-velocity_pt[0],0.};
            E_Float nrm_dir_ortho = abs(dir_ortho);
            if (nrm_dir_ortho < 1.E-14)
            {
                dir_ortho = vector3d{0.,-velocity_pt[2],velocity_pt[1]};
                nrm_dir_ortho = abs(dir_ortho);
                assert(nrm_dir_ortho < 1.E-14);
            }
            dir_ortho = (width/nrm_dir_ortho)*dir_ortho;
            streamPt( istream+1, pos_crds[0] ) = streamPt(istream, pos_crds[0]) + dir_ortho[0];
            streamPt( istream+1, pos_crds[1] ) = streamPt(istream, pos_crds[1]) + dir_ortho[1];
            streamPt( istream+1, pos_crds[2] ) = streamPt(istream, pos_crds[2]) + dir_ortho[2];
        }
        else
        {
#if defined( DEBUG_VERBOSE )
        std::cout << "Copie du vecteur orthogonal" << std::endl;
#endif
            //std::cout << "Doublement des points par changement de frontière" << std::endl;
            // Le point est doublé à la frontière entre les deux ensembles, on 
            streamPt( istream+1, pos_crds[0] ) = streamPt( istream-1, pos_crds[0]);
            streamPt( istream+1, pos_crds[1] ) = streamPt( istream-1, pos_crds[1]);
            streamPt( istream+1, pos_crds[2] ) = streamPt( istream-1, pos_crds[2]);
            dir_ortho = vector3d(streamPt( istream+1, pos_crds[0] )-streamPt( istream, pos_crds[0] ),
                                 streamPt( istream+1, pos_crds[1] )-streamPt( istream, pos_crds[1] ),
                                 streamPt( istream+1, pos_crds[2] )-streamPt( istream, pos_crds[2] ));
        }
        //return {num_blk, icell, velocity_pt};
        return std::tuple<E_Int, E_Int, vector3d, vector3d>( num_blk, icell, velocity_pt, dir_ortho );
    }

    //_ __________________ Construit la ribbon stream _________________________
    void build_ribbon_streamline( const point3d &init_pos,
                                 const std::vector<K_POST::zone_data_view> &zones,
                                 E_Int max_vertices_for_streamline,
                                 FldArrayF &streamPt, E_Float ribbon_width = 1., bool is_reversed = false )
    {
        E_Int ind_cell = -1;
        E_Int num_blk = -1;
        int istream;
        vector3d velocity, dir_ortho;
        point3d cur_point = init_pos;
        for ( istream = 0; istream < max_vertices_for_streamline; ++istream ) {
            //? On est dans un cas où on ne connaît pas la cellule d'appartenance du point courant
            //? de la ribbon stream. On va donc effectuer une recherche de la cellule ( et de la zone
            //? correspondante )
            if ( ind_cell == -1 ) {
#if defined( DEBUG_VERBOSE )
                std::cout << "Recherche d'un nouveau bloc contenant le point " << std::string( cur_point )
                          << " pour le " << istream << " sommet de streamline" << std::endl;
#endif
                std::tie( num_blk, ind_cell, velocity, dir_ortho ) =
                    init_ribbon_streamline( cur_point, zones, max_vertices_for_streamline, ribbon_width, 
                                           streamPt, 2*istream, num_blk );
#if defined( DEBUG_VERBOSE )
                std::cout << "bloc trouvé : " << num_blk << " indice cellule trouvée : " << ind_cell << std::endl;
#endif
                if ( ind_cell == -1 ) break;  // On n'a pas trouvé de zones correspondantes...
                continue;
            }
#if defined( DEBUG_VERBOSE )
            E_Int facette_in = -1;
#endif
            E_Int facette_out = -1;
            auto facettes_candidates = zones[ num_blk ].get_faces_of_element( ind_cell, num_blk );

#if defined( DEBUG_VERBOSE )
            std::cout << "position initiale : " << std::string( cur_point ) << std::endl;
            std::cout << "Vitesse initiale  : " << std::string( velocity ) << std::endl;
            std::cout << "Facettes candidates : " << std::endl;
            for ( const auto &f : facettes_candidates ) {
                const auto &zone_coords = f.coordinates_of_zone;
                std::cout << "\tindice des coordonnées : ";
                for ( auto fi : f.indices_vertices ) std::cout << fi << " ";
                std::cout << std::flush << std::endl;
                std::cout << "\tCoordonnées des sommets ";
                for ( E_Int ivert = 0; ivert < f.indices_vertices.size(); ivert++ )
                    std::cout << std::string( point3d{ zone_coords[ 0 ][ f.indices_vertices[ ivert ] ],
                                                       zone_coords[ 1 ][ f.indices_vertices[ ivert ] ],
                                                       zone_coords[ 2 ][ f.indices_vertices[ ivert ] ] } )
                              << " ";
                std::cout << std::endl;
            }
#endif
            if ( ( velocity | velocity ) < 1.E-14 ) {
#if defined( DEBUG_VERBOSE )
                std::cout << "Warning: vitesse nulle. fini calcul Streamline ..." << std::flush << std::endl;
#endif
                //istream -= 1;
                break;
            }
            //
            if ( is_reversed ) velocity = -velocity;
            //
            bool is_intersecting, is_entering;
            E_Int ind_facette = 0;
            for ( const auto &facette : facettes_candidates ) {
                std::tie( is_intersecting, is_entering ) = facette.is_intersecting_ray( cur_point, velocity );
                if ( is_intersecting ) {
                    if ( is_entering )
                    {
#if defined( DEBUG_VERBOSE )
                        facette_in = ind_facette;
#endif
                    }
                    else
                        facette_out = ind_facette;
                }

                ind_facette++;
            }
#if defined( DEBUG_VERBOSE )
            std::cout << "Facette sortante n°" << facette_out << " trouvée pour intersection rayon : " << std::flush << std::endl;
            auto &out_facette = facettes_candidates[ facette_out ];
            std::cout << "\tindice des coordonnées : ";
            for ( auto fi : out_facette.indices_vertices ) std::cout << fi << " ";
            std::cout << std::flush << std::endl;
            std::cout << "Facette rentrante n°" << facette_in << " trouvée pour intersection rayon ! " << std::flush << std::endl;
            auto &in_facette = facettes_candidates[ facette_in ];
            std::cout << "\tindice des coordonnées : ";
            for ( auto fi : in_facette.indices_vertices ) std::cout << fi << " ";
            std::cout << std::flush << std::endl;
#endif
            // Dans un premier temps, on recherche en aval. On cherchera en amont uniquement si is_bidirectional est vrai
            auto &facette = facettes_candidates[ facette_out ];
#if defined( DEBUG_VERBOSE )
            std::cout << "facette sélectionnée : " << std::endl;
            std::cout << "\tindice des coordonnées : ";
            for ( auto fi : facette.indices_vertices ) std::cout << fi << " ";
            std::cout << std::flush << std::endl;
#endif
            decltype( facette.compute_intersection( cur_point, velocity ) ) intersect_data;
#if !defined( DEBUG_VERBOSE )
            try {
#endif
                intersect_data = facette.compute_intersection( cur_point, velocity );
#if defined( DEBUG_VERBOSE )
                std::cout << "Nouvelle position : " << std::string( intersect_data.first ) << std::endl;
                std::cout << "No triangle d'intersection : " << intersect_data.second << "." << std::endl;
#endif
#if !defined( DEBUG_VERBOSE )
            } catch ( std::domain_error &err ) {
                std::cerr << "Warning: probleme intersection géométrique sur facette : " << err.what() << ". On arête prématurément le calcul de cette streamline." << std::endl;
                break;
            }
#endif
            const FldArrayF &field = *zones[ num_blk ].get_fields();
            std::vector<E_Float> interpfld = facette.compute_interpolated_field( intersect_data, field );
            auto pos_crds = zones[ num_blk ].get_position_of_coordinates();
            auto pos_velo = zones[ num_blk ].get_position_of_velocity();
            //? Quelles sont les cellules de la facette ?
            E_Int icell1 = facette.indices_cells.first;
            E_Int icell2 = facette.indices_cells.second;
            icell2 = (icell2 == -1 ? icell1 : icell2);
#if defined( DEBUG_VERBOSE )
            std::cout << "cell 1 : " << icell1 << ", cell 2 : " << icell2 << std::endl;
#endif
            //?  On lit les indices des sommets des deux cellules :
            std::vector<E_Int> ind_verts_cell1 = zones[ num_blk ].get_indices_of_vertices(icell1); 
#if defined( DEBUG_VERBOSE )
            std::cout << "Indices des sommets associés aux deux cellules. Première cellule : ";
            for ( const auto& iv : ind_verts_cell1)
                std::cout << iv << " ";
            std::cout << std::flush << std::endl;
#endif
            std::vector<E_Int> ind_verts_cell2 = zones[ num_blk ].get_indices_of_vertices(icell2); 
#if defined( DEBUG_VERBOSE )
            std::cout  << "Indices deuxième cellule : ";
            for ( const auto& iv : ind_verts_cell2)
                std::cout << iv << " ";
            std::cout << std::flush << std::endl;
#endif
            //? Calcul du rotationnel et du volume pour les deux cellules :
#if defined( DEBUG_VERBOSE )
            std::cout << "Calcul rotationel dans première cellule : " << std::flush;
#endif            
            vector3d rot1 = zones[ num_blk ].compute_rotational_in_cell(icell1);
#if defined( DEBUG_VERBOSE )
            std::cout << "rot1 : " << rot1 << std::endl;
            std::cout << "Calcul rotationel dans deuxième cellule : " << std::flush;
#endif
            vector3d rot2 = zones[ num_blk ].compute_rotational_in_cell(icell2);
#if defined( DEBUG_VERBOSE )
            std::cout << "rot2 : " << rot2 << std::flush << std::endl;
            std::cout << "Calcul Volume dans première cellule : " << std::flush;
#endif
            E_Float vol1 = zones[ num_blk ].compute_volume_of_cell(icell1);
#if defined( DEBUG_VERBOSE )
            std::cout << vol1 << std::flush << std::endl;
            std::cout << "Calcul volume dans deuxième cellule : "; 
#endif
            E_Float vol2 = zones[ num_blk ].compute_volume_of_cell(icell2);
#if defined( DEBUG_VERBOSE )
            std::cout << vol2 << std::flush << std::endl;
#endif
            //std::cout << "vol 1 : " << vol1 << ", vol 2 : " << vol2 << " ";
            vector3d rot = (1./(vol1+vol2))*(vol1*rot1 + vol2*rot2);
            //? Calcul de l'angle de rotation du vecteur normal en fonction  du
            //? rotationel :
            //$                                                →               
            //$  dθ        → →     →                      →    u               
            //$ ──── = 0.5(ω.s) où ω est le rotationel et s = ───              
            //$  dt                                            →               
            //$                                               ∥u∥              
            //? avec u le vecteur vitesse
            //? Calcul du temps de  parcours  de  la  particule  de  l'ancienne 
            //? position à la nouvelle en fonction des deux points pᵢ et pᵢ₊₁ de
            //? la streamline et de la vitesse moyenne u de la particule :
            // $                  →                                            
            // $  t = d(pᵢ,pᵢ₊₁)/∥u∥                                           
            // $                                                               ⁣
            //@ Voir article Fast Algorithms for Visualizing  Fluild Motion  in
            //@              steady flow on unstructured grids                 
            //@              Authors : S.K. Uengn K. Sikorski and Kwan-Liu Ma  
            //@              Sixth IEEE visualization conference               
            //@              1995                                              
            point3d pi{streamPt( 2*(istream-1), pos_crds[0]),
                       streamPt( 2*(istream-1), pos_crds[1]),
                       streamPt( 2*(istream-1), pos_crds[2])};
            point3d pip1{interpfld[pos_crds[0]-1],
                         interpfld[pos_crds[1]-1],
                         interpfld[pos_crds[2]-1]};
            E_Float d = distance(pi,pip1);
#if defined( DEBUG_VERBOSE )
            std::cout << "istream : " << istream << "pi : " << pi << ", pi+1 : " << pip1 << " ";
            std::cout << "d : " << d << " " << std::flush;
#endif
            vector3d ui{streamPt( 2*(istream-1), pos_velo[0]),
                        streamPt( 2*(istream-1), pos_velo[1]),
                        streamPt( 2*(istream-1), pos_velo[2])};
            vector3d uip1{interpfld[pos_velo[0]-1],
                          interpfld[pos_velo[1]-1],
                          interpfld[pos_velo[2]-1] };
#if defined( DEBUG_VERBOSE )
            std::cout << "u_i : " << std::string(ui) << ", et u_ip1 : " << std::string(uip1) << std::flush << std::endl;
#endif
            vector3d u = 0.5*(ui + uip1);
            //std::cout << "u : " << u << ", rot : " << rot << std::endl;
            E_Float nrmu = abs(u);
            u = (1./nrmu)*u;
            //std::cout << "||u||/||u|| ! : " << abs(u) << " ";
            E_Float t = d / nrmu;
            //std::cout << " t <= " << t << " ";
            // ? Calcul de l'angle en fonction de t :
            E_Float theta = 0.5*(rot|u)*t;
            //std::cout << "theta : " << theta << std::endl
             //         << "------------------------------------------------------------------------" << std::endl;
            // ? Orthogonalisation de dir_ortho par rapport à u :
            dir_ortho = dir_ortho - (dir_ortho|u)*u;
            // ? Renormalisation pour conserver la largeur de la bande
            dir_ortho = (ribbon_width/abs(dir_ortho))*dir_ortho;
            //std::cout << "<d|u> = " << (dir_ortho|u) << std::endl;
            // ? Puis rotation selon la divergence
            dir_ortho = rotate(u, theta, dir_ortho);
            //std::cout << "<d|u> = " << (dir_ortho|u) << std::endl;
            for ( decltype( interpfld.size() ) f = 0; f < interpfld.size(); ++f )
            {
                streamPt( 2*istream  , f + 1 ) = interpfld[ f ];
                streamPt( 2*istream+1, f + 1 ) = interpfld[ f ];
            }
            streamPt(2*istream+1, pos_crds[0] ) = interpfld[ pos_crds[0]-1 ] + dir_ortho.x;
            streamPt(2*istream+1, pos_crds[1] ) = interpfld[ pos_crds[1]-1 ] + dir_ortho.y;
            streamPt(2*istream+1, pos_crds[2] ) = interpfld[ pos_crds[2]-1 ] + dir_ortho.z;
            auto npos_vel = zones[ num_blk ].get_position_of_velocity();
            velocity = { streamPt( 2*istream, npos_vel[ 0 ] ), streamPt( 2*istream, npos_vel[ 1 ] ), 
                         streamPt( 2*istream, npos_vel[ 2 ] ) };
            ind_cell = facette.indices_cells.second;
            cur_point = intersect_data.first;
        }  // End for (istream)
        E_Int nfld = streamPt.getNfld();
        streamPt.reAllocMatSeq( 2*istream, nfld );
    }
}  // namespace


namespace K_POST
{
    //# #######################################################################
    //# #       Constructeur de la classe privée d'implémentation             #
    //# #######################################################################
    //_ _____________________________ Constructeur ____________________________
    ribbon_streamline::ribbon_streamline::Implementation::Implementation( 
        const point3d &init_pos, const std::vector<zone_data_view> &zones,
        E_Int max_vertices_for_streamline, E_Float width, bool is_bidirectional )
    {
        if (!is_bidirectional)
            build_ribbon_streamline(init_pos, zones, max_vertices_for_streamline, 
                                    this->fields, width);
        else
        {
            FldArrayF streamPt;
            build_ribbon_streamline(init_pos, zones, max_vertices_for_streamline, 
                                    streamPt, width);
            FldArrayF streamPtRev;
            build_ribbon_streamline(init_pos, zones, max_vertices_for_streamline, 
                                    streamPtRev, width, true);
            this->fields.malloc(streamPt.getSize()+streamPtRev.getSize()-2,
                                streamPt.getNfld());
            for ( E_Int f = 0; f < streamPt.getNfld(); ++f )
            {
                for ( E_Int i = 0; i < streamPt.getSize(); ++i )
                    this->fields(i + streamPtRev.getSize()-2, f+1) = streamPt(i,f+1);
                for ( E_Int i = 2; i < streamPtRev.getSize(); i += 2 )
                {
                    this->fields(streamPtRev.getSize()-1-i,f+1) = streamPtRev(i+1,f+1);
                    this->fields(streamPtRev.getSize()-2-i,f+1) = streamPtRev(i  ,f+1);
                }
            }
        }
    }
    //# #######################################################################
    //# #            Constructeur et méthodes de la classe publique           #
    //# #######################################################################
    //_ _____________________________ Constructeur ____________________________
    ribbon_streamline::ribbon_streamline( 
        const point3d& init_pos, const std::vector<zone_data_view>&  zones, 
        E_Int max_vertices_for_streamline, E_Float width, bool is_bidirectional) :
        ptr_impl(new Implementation(init_pos, zones, max_vertices_for_streamline, 
                                    width, is_bidirectional))
    {}
    //  _____________________________ Destructeur _____________________________
    ribbon_streamline::~ribbon_streamline()
    {
        delete ptr_impl;
    }
    //_ ________ Retourne le champ construit pour la ribbon stream ____________
    const FldArrayF& ribbon_streamline::field() const
    {
        return this->ptr_impl->fields;
    }
    //_ ________ Retourne le champ construit pour la ribbon stream ____________
    FldArrayF& ribbon_streamline::field()
    {
        return this->ptr_impl->fields;
    }
}  //! namespace K_POST
    
