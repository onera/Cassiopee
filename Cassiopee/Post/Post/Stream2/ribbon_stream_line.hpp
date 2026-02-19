#ifndef _POST_RIBBON_STREAM_LINE_HPP_
#define _POST_RIBBON_STREAM_LINE_HPP_
#include "stream_line.hpp"

namespace K_POST
{
    //# ################################################################
    //# #         Définition d'une classe ribbon streamline            #
    //# ################################################################

    // @brief      Describe a ribbon stream line                        
    class ribbon_streamline
    {
    public:
        //________________ Constructeurs et destructeur ________________

        // @brief      Constructs a new instance.                       
        //
        // @param[in]  init_pos          The initialize position        
        // @param[in]  zones             The zones                      
        // @param[in]  width             Width of the ribbon stream     
        // @param[in]  is_bidirectional  Indicates if bidirectional     
        // 
        // @details Construit une ribbon streamline à partir d'un  point
        //-         initial donné dans l'espace.                        
        //-         Si l'option is_bidirectional est vrai, on  construit
        //-         la ribbon streamline dans les deux sens, en aval  et
        //-         en amont du flux, sinon on ne construit que dans  le
        //-         sens amont.                                         
        ribbon_streamline( const point3d& init_pos, 
                           const std::vector<zone_data_view>&  zones, 
                           E_Int max_vertices_for_streamline, 
                           E_Float width=1., 
                           bool is_bidirectional=false);

        // @brief          Destructeur                                  
        // 
        // @details Détruit l'implémentation pointée par la classe      
        ~ribbon_streamline();

        //_________________ Accesseurs et modifieurs ___________________

        // @brief Retourne le champ correspondant à la streamline.      
        // 
        // @details  Le champ contient également les coordonnées        
        // 
        const FldArrayF& field() const;

        // @brief Retourne le champ correspondant à la streamline       
        // 
        //  @details Le champs contient également les coordonnées       
        // 
        FldArrayF& field();

        //~                  Partie privée de la classe                 
    private:
        struct Implementation;
        Implementation* ptr_impl;
    };
}

#endif
