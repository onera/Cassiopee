/**
 * Calcul du volume de solide en trois dimensions en utilisant le théorème de la divergence :
 #     → ⎛x⎞         ⎛x⎞ ⟶
 - ∭ div⎜y⎟.dΩ = ∯ ⎜y⎟.dΓ = 3.vol(Ω)
 -   Ω   ⎝z⎠       Γ ⎝z⎠
 * 
 */
#ifndef _POST_STREAM2_VOLUME_HPP_
#define _POST_STREAM2_VOLUME_HPP_
#include "point3d.hpp"
#include "face.hpp"

namespace K_POST {
/** 
 ! @brief Calcule le volume d'un tétraèdre défini par ses quatre sommets
 * @details Calcule le volume d'un tétraèdre défini par ses quatre sommets
 * 
 * On utilise le théorème de la divergence (voir en entête de fichier).
 * 
 * Soit le tétraèdre défini par les sommets p₁, p₂, p₃ et p₄ (définis dans l'ordre convenu dans CGNS)
 *
 *         p₄     
 *        ╱│╲     
 *       ╱ │ ╲p₃  
 *    p₁╱__│╱     
 *          p₂    
 *
 * Pour calculer le volume, on fait un changement de repère tel que p₁ est à l'origine
 * du nouveau repère.
 *
 * Sa fonction de transformation de l'espace paramétrique dans l'espace cartésien est :
 * 
 * ϕ : [0;1]³ ----> ℝ³                                           
 *     ψ,η,ζ  ----> ψ.E₁₂ + η.E₁₃ + ζ.E₁₄ avec 0 ≤ ψ + η + ζ ≤ 1 
 *  où Eᵢⱼ est le vecteur allant de pᵢ à pⱼ.                     
 *
 * Pour le triangle T₁₃₂ = (p₁, p₃, p₂), on a ζ = 0 et sa normale est (contenant la surface du triangle ) 
 * 
 *      n₁₃₂ = E₁₃ ∧ E₁₂  
 * 
 * Alors      
 *      ⎛x⎞ ⟶      1 1-η                                     
 * ∯   ⎜y⎟.dT₁₃₂ = ∫ ∫ (η.E₁₃ + ψ.E₁₂).(E₁₃ ∧ E₁₂) dψ dη = 0 
 * T₁₂₃ ⎝z⎠         0 0                                       
 * 
 * Pour le triangle T₁₂₄ = (p₁,p₂,p₄), on a η = 0 et sa normale   
 * 
 *      n₁₂₄ = E₁₂ ∧ E₁₄   
 * 
 * Alors   
 * 
 *      ⎛x⎞  ⟶     1 1-ζ                                      
 * ∯   ⎜y⎟.dT₁₂₄ = ∫ ∫ (ψ.E₁₂ + ζ.E₁₄).(E₁₂ ∧ E₁₄) dψ dζ = 0  
 * T₁₂₄ ⎝z⎠         0 0                                        
 * 
 * Pour le triangle T₂₃₄ = (p₂, p₃, p₄), on a ψ + η + ζ = 1 soit ζ = 1 - ψ - η et la normale :
 * 
 *      n₂₃₄ = E₂₃ ∧ E₂₄                                                                      
 * 
 * Alors    
 * 
 *      ⎛x⎞  ⟶     1 1-η                                              
 * ∯   ⎜y⎟.dT₂₃₄ = ∫ ∫ (ψ.E₁₂ + η.E₁₃ + (1-ψ-η).E₁₄).(E₂₃ ∧ E₂₄) dψ dη
 * T₂₃₄ ⎝z⎠         0 0                                                
 * 
 *                 1 1-η                                                             
 *               =∫ ∫ (ψ.E₄₂ + η.E₄₃ + E₁₄).(E₂₃ ∧ E₂₄) dψ dη = (1/2)E₁₄.(E₂₃ ∧ E₂₄)
 *                0 0                                                                
 * 
 * Si bien (puisque la dernière intégrale à calculer est nulle) qu'en définitif, le volume du tétraèdre vaut :
 * 
 * V = (1/6)E₁₄.(E₂₃ ∧ E₂₄) = (1/6)[E₁₄;E₂₃;E₂₄] où [;;] est le produit mixte
 *                
 * @param p1 Premier sommet du tétraèdre   
 * @param p2 Second sommet du tétraèdre    
 * @param p3 Troisième sommet du tétraèdre 
 * @param p4 Quatrième sommet du tétraèdre 
 * @return Retourne le volume du tétraèdre 
 */
E_Float compute_volume_tetrahedra( const point3d &p1, const point3d &p2, const point3d &p3, const point3d &p4 );

/**
 * @brief Calcule le volume d'une pyramide courbe (considéré ici comme un polytope et non un polyèdre)
 * @details Calcule le volume d'une pyramide dont la base est (p₁, p₄, p₃, p₂) et de sommet p₅ (convention CGNS) :
 * 
 * En suivant les recommandations de Bergot (thèse), on prend pour pyramide de référence la pyramide
 * dont la base est centrée en 0 (de section [-1;1]²) et de sommet (0,0,1)
 * 
 * La fonction permettant de passer de la pyramide de référence à la pyramide actuelle est une fonction de la forme :
 * 
 # F = p₁.ϕ₁ + p₂.ϕ₂ + p₃.ϕ₃ + p₄.ϕ₄ + p₅.ϕ₅                                   
 # avec ϕ₁ = (1/4)(1-ψ-η-ζ + (ψη)/(1-ζ) )                                      
 #      ϕ₂ = (1/4)(1+ψ-η-ζ - (ψη)/(1-ζ) )                                      
 #      ϕ₃ = (1/4)(1+ψ+η-ζ + (ψη)/(1-ζ) )                                      
 #      ϕ₄ = (1/4)(1-ψ+η-ζ - (ψη)/(1-ζ) )                                      
 #      ϕ₅ = ζ                                                                 
 # et (ψ,η,ζ) ∈ [-1;1]²⨯[0;1]                                                  
 * 
 * Pour la simplification des calculs, on effectue un changement de repère tel que : p₁+p₂+p₃+p₄ = 0. La fonction
 * F peut alors s'exprimer comme :
 * 
 # 4.F = ψ.(E₁₂-E₃₄) + η.(E₁₃+E₂₄) + ζ.(E₁₅+E₂₅+E₃₅+E₄₅) + ψη/(1-ζ).(E₂₃-E₁₄)  
 * 
 *  où Eᵢⱼ est le vecteur allant de pᵢ à pⱼ. 
 * 
 * Pour le quadrangle Q₁₄₃₂=(p₁, p₄, p₃, p₂), on impose ζ=0 et la normale est définie comme :
 * 
 #      n₁₄₃₂ = (1/16)[E₁₂-E₃₄ + η.(E₂₃-E₁₄)] ∧ [E₁₃+E₂₄ + ψ.(E₂₃-E₁₄)]        
 #            = (1/8)[E₁₃∧E₂₄ + ψ.E₁₂∧E₄₃ + η.E₂₃∧E₁₄]                        
 * 
 * Alors
 * 
 #       ⎛x⎞  ⟶            1  1                                                                                    
 # ∯    ⎜y⎟.dQ₁₄₃₂ =(1/32) ∫  ∫ [ψ.(E₁₂-E₃₄) + η.(E₁₃+E₂₄) + ψη.(E₂₃-E₁₄)].[E₁₃∧E₂₄ + ψ.E₁₂∧E₄₃ + η.E₂₃∧E₁₄] dψ dη 
 # Q₁₄₃₂ ⎝z⎠              -1 -1                                                                                      
 #                 = 0 (car il ne reste plus qu'un terme en ψη après simplification dans l'intégrale)               
 *                 
 * Pour le triangle T₁₂₅, on impose η=-1+ζ, et la normale vaut :
 * 
 #      n₁₂₅ = (1/2)E₁₂∧E₁₅                                                    
 *      
 * Alors
 * 
 #      ⎛x⎞  ⟶          1  1-ζ                                                                                   
 # ∯   ⎜y⎟.dT₁₂₅ =(1/8) ∫  ∫ [ψ.(E₁₃+E₄₂) + E₃₁+E₄₂ + ζ.(E₁₃+E₂₄+E₁₅+E₂₅+E₃₅+E₄₅) - ψ.(E₃₂+E₁₄)].[E₁₂∧E₁₅] dψ dζ 
 # T₁₂₅ ⎝z⎠             0 -1+ζ                                                                                    
 #               =(1/8)[E₁₂;E₁₃+E₁₄;E₁₅]                                                                         
 *               
 * Pour le triangle T₂₃₅, on impose ψ= 1-ζ et la normale vaut :
 * 
 #      n₂₃₅ = (1/2)E₂₃∧E₂₅                                                    
 * 
 * Alors
 * 
 #      ⎛x⎞  ⟶          1  1-ζ                                                                                   
 # ∯   ⎜y⎟.dT₂₃₅ =(1/8) ∫  ∫ [E₁₃+E₄₂ + η.(E₁₃+E₂₄) + ζ.(E₃₁+E₂₄+E₁₅+E₂₅+E₃₅+E₄₅) + η.(E₂₃-E₁₄)].[E₂₃∧E₂₅] dψ dζ 
 # T₂₃₅ ⎝z⎠             0 -1+ζ                                                                                    
 #               =(1/8){[E₁₂;E₁₃;E₁₅]+[E₂₃;E₂₄;E₂₅]}                                                             
 *               
 * Pour le triangle T₃₄₅, on impose η=+1-ζ et la normale vaut :
 * 
 #         n₃₄₅ = (1/2)E₃₄∧E₃₅                                                 
 * 
 * Alors
 * 
 #      ⎛x⎞  ⟶          1  1-ζ                                                                                  
 # ∯   ⎜y⎟.dT₃₄₅ =(1/8) ∫  ∫ [ψ.(E₁₂-E₃₄) + E₁₃+E₂₄ + ζ.(E₃₁+E₄₂+E₁₅+E₂₅+E₃₅+E₄₅) + ψ.(E₂₃-E₁₄)].[E₃₄∧E₃₅] dψ dζ
 # T₃₄₅ ⎝z⎠             0 -1+ζ                                                                                   
 #               = (1/8){[E₁₃;E₁₄;E₁₅]+[E₂₃;E₂₄;E₂₅]}                                                           
 * 
 * Pour le triangle T₄₁₅, on impose ψ=-1+ζ et la normale vaut :
 * 
 #         n₄₁₅ = (1/2)E₄₁∧E₄₅                                                 
 *         
 * Alors
 * 
 #      ⎛x⎞  ⟶           1  1-ζ                                                                                  
 # ∯   ⎜y⎟.dT₄₁₅ =(1/8) ∫  ∫ [E₂₁+E₃₄ + η.(E₁₃+E₂₄) + ζ.(E₁₂+E₄₃+E₁₅+E₂₅+E₃₅+E₄₅) - η.(E₃₂+E₁₄)].[E₄₁∧E₄₅] dη dζ 
 # T₄₁₅ ⎝z⎠              0 -1+ζ                                                                                   
 #               =(1/8)[E₁₂+E₁₃;E₁₄;E₁₅]                                                                         
 * 
 * Si bien que le volume du polytope pyramidoïde est :
 * 
 ? V(P) = (1/12){[E₁₂;E₁₃;E₁₅]+[E₁₂+E₁₃;E₁₄;E₁₅]+[E₂₃;E₂₄;E₂₅]}
 * 
 * où l'opérateur [u;v;w] est le produit mixte des vecteurs u, v et w.
 * 
 * @param p1 Premier sommet de la base de la pyramide
 * @param p2 Second sommet de la base de la pyramide
 * @param p3 Troisième sommet de la base de la pyramide
 * @param p4 Quatrième sommet de la base de la pyramide
 * @param p5 Cinquième sommet de la pyramide
 * @return Le volume de la pyramide
 */
E_Float compute_volume_pyramid( const point3d &p1, const point3d &p2, const point3d &p3,
                                const point3d &p4, const point3d &p5);
/**
 * @brief Calcule le volume d'un pentaedre
 * @details Calcule le volume d'un pentaedre en se basant sur le théorème de la divergence
 * 
 * p₄____p₆
 * │╲   ╱│   Soit le pentaèdre (p₁,p₂,p₃,p₄,p₅,p₆) dont le sommets sont ordonnés selon la
 * │ ╲ ╱ │   la norme CGNS. Par changement de repère, on se ramène à p₁ = (0,0,0)
 * │  p₅ │
 * │  │  │   On note Eᵢⱼ le vecteur ayant pour origine le sommet pᵢ et pour extrémité le sommet pⱼ.
 * p₁-│--p₃
 *  ╲ │ ╱    La fonction ϕ permettant de passer de l'espace paramétrique ( sur le pentaèdre de référence )
 *   ╲│╱    au pentaèdre dans l'espace cartésien est :
 *    p₂     ϕ : [0;1]³ ----> ℝ³
 *               ψ,η,ζ  ----> ψ.E₁₂ + η.E₁₃ + ζ.E₁₄ + ψζ.(E₄₅ + E₂₁) + ηζ.(E₄₆ + E₃₁) avec 0 ≤ ψ + η ≤ 1
 * Calculons le flux de volume passant par le quadrangle Q₁₂₅₄ = (p₁,p₂,p₅,p₄) :
 * 
 * Sur Q₁₂₅₄, on fixe η = 0 et la norme sur ce quadrangle est décrite par la fonction :
 * 
 *       n₁₂₅₄ = (E₁₂ + ζ.(E₄₅ + E₂₁)) ∧ (E₁₄ + ψ.(E₄₅ + E₂₁))
 *             = E₁₂∧E₁₄ + ψ.E₁₂∧E₄₅ + ζ.E₂₅∧E₁₄
 *             
 * Le flux vaut donc :
 * 
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₁₂₅₄ = ∫  ∫ [ψ.E₁₂ + ζ.E₁₄ + ψζ.(E₄₅ + E₂₁)].[E₁₂∧E₁₄ + ψ.E₁₂∧E₄₅ + ζ.E₂₅∧E₁₄] dψ dζ
 * Q₁₂₅₄ ⎝z⎠         0  0
 * 
 *                 = (1/4)[E₁₂;E₁₅;E₂₄]
 *                 
 * Calculons le flux du volume passant par le quadrangle Q₂₃₆₅ = (p₂,p₃,p₆,p₅) :
 * 
 * Sur Q₂₃₆₅, on fixe η = 1-ψ et la norme sur ce quadrangle est décrite par la fonction :
 * 
 *      n₂₃₆₅ = (E₃₆ + ψ.(E₆₅+E₂₃)) ∧ (E₃₂ + ζ.(E₆₅+E₂₃))
 *            = E₃₆∧E₃₂ + ψ.E₆₅∧E₃₂ + ζ.E₃₆∧E₂₅
 * 
 * Le flux vaut donc :
 * 
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₂₃₆₅ = ∫  ∫ [E₁₃ + ψ.E₃₂ + ζ.E₃₆ + ψζ.(E₆₅+E₂₃)].[E₃₆∧E₃₂ + ψ.E₆₅∧E₃₂ + ζ.E₃₆∧E₂₅] dψ dζ
 * Q₂₃₆₅ ⎝z⎠         0  0
 * 
 *                 = (1/2)[E₁₃;E₂₆;E₃₅] + (1/4)[E₂₆;E₆₅;E₃₅]
 *                 
 * Calculons le flux du volume passant par le quadrangle Q₃₁₄₆ = (p₃,p₁,p₄,p₆) :
 * 
 * Sur Q₃₁₄₆, on fixe ψ = 0 et la norme sur ce quadrangle est décrite par la fonction :
 * 
 *      n₃₁₄₆ = (E₁₄ + η.(E₄₆ + E₃₁)) ∧ (E₁₃ + ζ.(E₄₆ + E₃₁))
 *            = E₁₄∧E₁₃ + η.E₄₆∧E₁₃ + ζ.E₁₄∧E₃₆
 *            
 * Le flux vaut donc :
 * 
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₃₁₄₆ = ∫  ∫ [η.E₁₃ + ζ.E₁₄ + ηζ.(E₄₆ + E₃₁)].[E₁₄∧E₁₃ + η.E₄₆∧E₁₃ + ζ.E₁₄∧E₃₆] dη dζ
 * Q₃₁₄₆ ⎝z⎠         0  0
 * 
 *                 = (1/4)[E₁₃; E₃₄; E₁₆]
 *                 
 * Calculons le flux du volume passant par le triangle T₁₃₂ = (p₁,p₃,p₂) :
 * 
 * Sur T₁₃₂, on a ζ = 0 et la normale vaut (c'est une constante) :
 * 
 *      n₁₃₂ = E₁₃∧E₁₂
 *      
 * Le flux vaut donc :
 *      
 *      ⎛x⎞  ⟶     1  1-η
 * ∯   ⎜y⎟.dT₁₃₂ = ∫  ∫ [ψ.E₁₂ + η.E₁₃].[E₁₃∧E₁₂] dψ dη = 0
 * T₁₃₂ ⎝z⎠        0  0
 * 
 * Calculons le flux du volume passant par le triangle T₄₅₆ = (p₄,p₅,p₆) :
 * 
 * Sur T₄₅₆, on a ζ = 1 et la normale vaut (c'est une constante) :
 * 
 *      n₄₅₆ = E₄₅∧E₄₆
 *      
 *  Le flux vaut donc :
 *  
 *      ⎛x⎞  ⟶     1  1-η
 * ∯   ⎜y⎟.dT₄₅₆ = ∫  ∫ [E₁₄ + ψ.E₄₅ + η.E₄₆].[E₄₅∧E₄₆] dψ dη = (1/2)[E₁₄;E₄₅;E₄₆]
 * T₄₅₆ ⎝z⎠        0  0
 *  
 * Au total, on a donc le volume égal à :
 * 
 * V = (1/6){[E₁₄;E₄₅;E₄₆] + (1/2)[E₁₃; E₃₄; E₁₆] + [E₁₃;E₂₆;E₃₅] + (1/2)[E₂₆;E₆₅;E₃₅] + (1/2)[E₁₂;E₁₅;E₂₄]}
 * 
 * @param p1 Premier sommet du pentaèdre
 * @param p2 Second sommet du pentaèdre
 * @param p3 Troisième sommet du pentaèdre
 * @param p4 Quatrième sommet du pentaèdre
 * @param p5 Cinquième sommet du pentaèdre
 * @param p6 Sixième sommet du pentaèdre
 * @return Volume du polytope pentaédrique
 */
E_Float compute_volume_pentaedra( const point3d &p1, const point3d &p2, const point3d &p3,
                                  const point3d &p4, const point3d &p5, const point3d &p6);

/**
 * @brief Calcule le volume d'un polytope de type hexaèdre
 * @details Calcule le volume d'un polytope trilinéaire de type hexaèdre
 * 
 * L'hexaèdre H=(p₁,p₂,p₃,p₄,p₅,p₆,p₇,p₈) peut être défini par le fonction ϕ transformant l'hexaèdre
 * de référence [0;1]³ en l'hexaèdre H.
 * 
 *   p₈____________p₇       Pour simplifier les calculs, on effectue un changement de variable tel que
 *   ╱⁞           ╱│        le point p₁ coïncide avec l'origine du repère : p₁ = (0,0,0)
 *  ╱ ⁞          ╱ │
 * p₅___________p₆ │        ϕ : [0;1]³ ----> ℝ³
 * │  ⁞         │  │            ψ,η,ζ  ----> ψ.E₁₂ + η.E₁₄ + ζ.E₁₅ + ψη.(E₄₃+E₂₁) + ψζ.(E₅₆+E₂₁) +
 * │  p₄……………………│……p₃                        ηζ.(E₅₈+E₄₁) + ψηζ.(E₈₇+E₃₄+E₆₅+E₁₂)
 * │  /         │ ╱
 * │ /          │╱          avec Eᵢⱼ le vecteur ayant pour origine le sommet pᵢ et pour extrémité le sommet pⱼ.
 * p₁___________p₂
 * 
 * Calculons le flux du volume passant par le quadrangle Q₁₄₃₂ = (p₁,p₄,p₃,p₂) :
 * 
 * Sur Q₁₄₃₂, on fixe ζ=0 et la normale est en fonction de ψ et η :
 * 
 *      n₁₄₃₂ = (E₁₄ + ψ.(E₄₃+E₂₁))∧(E₁₂ + η.(E₄₃+E₂₁))
 *            = E₁₄∧E₁₂ + ψ.E₄₃∧E₁₂ + η.E₁₄∧E₂₃
 *      
 *  Le flux vaut donc :
 *  
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₁₄₃₂ = ∫  ∫ [ψ.E₁₂ + η.E₁₄ + ψη.(E₄₃+E₂₁)].[E₁₄∧E₁₂ + ψ.E₄₃∧E₁₂ + η.E₁₄∧E₂₃] dψ dη
 * Q₁₄₃₂ ⎝z⎠         0  0
 * 
 *                 = (1/4)[E₁₂;E₁₄;E₁₃]
 * 
 * Calculons le flux du volume passant par Q₁₂₆₅ = (p₁,p₂,p₆,p₅) :
 * 
 * Sur Q₁₂₆₅, on fixe η=0 et la normale est en fonction de ψ et ζ :
 * 
 *      n₁₂₆₅ = (E₁₂ + ζ.(E₅₆+E₂₁))∧(E₁₅ + ψ.(E₅₆+E₂₁))
 *            = E₁₂∧E₁₅ + ψ.E₁₂∧E₅₆ + ζ.E₂₆∧E₁₅
 * 
 *  Le flux vaut donc :
 *  
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₁₂₆₅ = ∫  ∫ [ψ.E₁₂ + ζ.E₁₅ + ψζ.(E₅₆+E₂₁)].[E₁₂∧E₁₅ + ψ.E₁₂∧E₅₆ + ζ.E₂₆∧E₁₅] dψ dζ
 * Q₁₂₆₅ ⎝z⎠         0  0
 * 
 *                 = (1/4)[E₁₅;E₁₂;E₁₆]
 * 
 * Calculons le flux du volume passant par le quadrangle Q₂₃₇₆ = (p₂,p₃,p₇,p₆) :
 * 
 * Sur Q₂₃₇₆, on fixe ψ = 1 et la normale est en fonction de η et ζ :
 * 
 *      n₂₃₇₆ = (E₂₃ + ζ.(E₆₇+E₃₂))∧(E₂₆ + η.(E₆₇+E₃₂))
 *            = E₂₃∧E₂₆ + η.E₂₃∧E₆₇ + ζ.E₃₇∧E₂₆
 *            
 *  Le flux vaut donc :
 *  
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₂₃₇₆ = ∫  ∫ [E₁₂ + η.E₂₃ + ζ.E₂₆ + ηζ.(E₆₇+E₃₂)].[E₂₃∧E₂₆ + η.E₂₃∧E₆₇ + ζ.E₃₇∧E₂₆] dη dζ
 * Q₂₃₇₆ ⎝z⎠         0  0
 * 
 *                 = (1/2)[E₁₂;E₆₃;E₂₇] + (1/4)[E₂₃;E₃₇;E₂₆]
 *
 * Calculons le flux du volume passant par le quadrangle Q₃₄₈₇ = (p₃,p₄,p₈,p₇) :
 * 
 * Sur Q₃₄₈₇, on fixe η = 1 et la normale s'exprime en fonction de ψ et ζ :
 * 
 *      n₃₄₈₇ = (E₄₈ + ψ.(E₈₇+E₃₄))∧(E₄₃ + ζ.(E₈₇+E₃₄))
 *            = E₄₈∧E₄₃ + ψ.E₈₇∧E₄₃ + ζ.E₄₈∧E₃₇
 *            
 *  Le flux vaut donc :
 *  
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₃₄₈₇ = ∫  ∫ [E₁₄ + ψ.E₄₃ + ζ.E₄₈ + ψζ.(E₈₇+E₃₄)].[E₄₈∧E₄₃ + ψ.E₈₇∧E₄₃ + ζ.E₄₈∧E₃₇] dψ dζ
 * Q₃₄₈₇ ⎝z⎠         0  0
 * 
 *                 = (1/2)[E₁₄;E₄₇;E₈₃] + (1/4)[E₄₃;E₄₈;E₃₇]  
 *            
 *  Calculons le flux du volume passant par le quadrangle Q₁₅₈₄ = (p₁,p₅,p₈,p₄) :
 *  
 *  Sur Q₁₅₈₄, on fixe ψ = 0 et la normale est en fonction de η et ζ :
 *  
 *      n₁₅₈₄ = (E₁₅ + η.(E₅₈+E₄₁))∧(E₁₄ + ζ.(E₅₈+E₄₁))
 *            = E₁₅∧E₁₄ + η.E₅₈∧E₁₄ + ζ.E₁₅∧E₄₈
 *            
 *  Le flux vaut donc :
 *  
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₁₅₈₄ = ∫  ∫ [η.E₁₄ + ζ.E₁₅ + ηζ.(E₅₈+E₄₁)].[E₁₅∧E₁₄ + η.E₅₈∧E₁₄ + ζ.E₁₅∧E₄₈] dη dζ
 * Q₁₅₈₄ ⎝z⎠         0  0
 * 
 *                 = (1/4)[E₁₅;E₅₈;E₁₄]
 *            
 *  Calculons le flux du volume passant par le quadrangle Q₅₆₇₈ = (p₅,p₆,p₇,p₈) :
 *  
 *  Sur Q₅₆₇₈, on fixe ζ=1 et la normale est en fonction de ψ et η :
 *  
 *      n₅₆₇₈ = (E₅₆ + η.(E₈₇+E₆₅))∧(E₅₈ + ψ.(E₈₇+E₆₅))
 *            = E₅₆∧E₅₈ + ψ.E₅₆∧E₈₇ + η.E₆₇∧E₅₈
 *            
 *  Le flux vaut donc :
 *  
 *       ⎛x⎞  ⟶      1  1
 * ∯    ⎜y⎟.dQ₅₆₇₈ = ∫  ∫ [E₁₅ + ψ.E₅₆ + η.E₅₈ + ψη.(E₈₇+E₆₅)].[E₅₆∧E₅₈ + ψ.E₅₆∧E₈₇ + η.E₆₇∧E₅₈] dψ dη
 * Q₅₆₇₈ ⎝z⎠         0  0
 *  
 *                 = (1/2)[E₁₅;E₈₆;E₅₇] + (1/4)[E₅₆;E₆₇;E₅₈] 
 * 
 * On en déduit le volume de l'hexaèdre polytope :
 * 
 * V = (1/6){[E₁₂;E₆₃;E₂₇] + [E₁₄;E₄₇;E₈₃] + [E₁₅;E₈₆;E₅₇]} + (1/12){[E₁₂;E₁₄;E₁₃] + [E₁₅;E₁₂;E₁₆] + [E₂₃;E₃₇;E₂₆] +
 *                                                                   [E₄₃;E₄₈;E₃₇] + [E₁₅;E₅₈;E₁₄] + [E₅₆;E₆₇;E₅₈]}
 *                                                                   
 * @param p1 Premier sommet de l'hexaèdre
 * @param p2 Second sommet de l'hexaèdre
 * @param p3 Troisième sommet de l'hexaèdre
 * @param p4 Quatrième sommet de l'hexaèdre
 * @param p5 Cinquième sommet de l'hexaèdre
 * @param p6 Sixième sommet de l'hexaèdre
 * @param p7 Septième sommet de l'hexaèdre
 * @param p8 Huitième sommet de l'hexaèdre
 * @return Le volume du polytope hexaèdrique.
 */
E_Float compute_volume_hexaedra ( const point3d &p1, const point3d &p2, const point3d &p3,
                                 const point3d &p4, const point3d &p5, const point3d &p6,
                                 const point3d &p7, const point3d &p8 );

/**
 * @brief Calcule le volume d'un N-Gon
 * @details Calcule le volume d'un polytope trilinéaire quelconque
  
   On suppose pour cette fonction que les faces du polytope ne sont constitués que de quadrangles et de triangles.
   Si une face est définie avec plus de quatre sommets, on suppose que l'utilisateur à d'abord tessaliser la face
   pour la décomposer en triangle. On suppose également que chaque face est orientée de telle sorte que la normale
   à la face soit orientée vers l'extérieur du polytope.
   
   L'idée pour calculer le volume est que la trace de la fonction "paramétrique" décrivant le polytope sur chaque face
   peut être décrit par une fonction paramétrique décrivant cette face (dans un repère fixé pour toutes les faces) et qu'à l'aide
   du théorème de la divergence, il est alors possible de calculer le volume.
   
   Ainsi, pour une face triangulaire Fᵢⱼₖ = (pᵢ,pⱼ,pₖ), la paramétrisation de la face est :
   
   ϕᵢⱼₖ : [0;1]² ----> ℝ³ avec 0 ≤ ψ+η ≤ 1
           ψ,η  ----> pᵢ + ψ.Eᵢⱼ + η.Eᵢₖ
           
   et la normale vaut :
   
        nᵢⱼₖ = Eᵢⱼ∧Eᵢₖ
        
   et le flux du volume au travers du triangle vaut : 
   
        ⎛x⎞  ⟶     1  1-η
   ∯   ⎜y⎟.dTᵢⱼₖ = ∫  ∫ [pᵢ + ψ.Eᵢⱼ + η.Eᵢₖ].[Eᵢⱼ∧Eᵢₖ] dψ dη = (1/2)[pᵢ;Eᵢⱼ;Eᵢₖ]
   Tᵢⱼₖ ⎝z⎠        0  0
   
   Pour une face quadrangulaire Fᵢⱼₖₗ = (pᵢ,pⱼ,pₖ,pₗ), la paramétrisation de la face est :
   
   ϕᵢⱼₖₗ : [0;1]² ----> ℝ³
            ψ,η  ----> pᵢ + ψ.Eᵢⱼ + η.Eᵢₗ + ψη.(Eⱼₖ + Eₗᵢ)
           
   et la normale est en fonctione de ψ et η :
   
        nᵢⱼₖₗ = (Eᵢⱼ + η.(Eⱼₖ + Eₗᵢ))∧(Eᵢₗ + ψ.(Eⱼₖ + Eₗᵢ))
             = Eᵢⱼ∧Eᵢₗ + ψ.Eᵢⱼ∧Eₗₖ + η.Eⱼₖ∧Eᵢₗ
             
   Le flux du volume au travers du quadrangle vaut :
   
        ⎛x⎞  ⟶     1  1
   ∯   ⎜y⎟.dQᵢⱼₖₗ = ∫  ∫ [pᵢ + ψ.Eᵢⱼ + η.Eᵢₗ + ψη.(Eⱼₖ + Eₗᵢ)].[Eᵢⱼ∧Eᵢₗ + ψ.Eᵢⱼ∧Eₗₖ + η.Eⱼₖ∧Eᵢₗ] dψ dζ
   Qᵢⱼₖₗ ⎝z⎠        0  0
   
                 = (1/2)[pᵢ;Eₗⱼ;Eᵢₖ] + (1/4)[Eᵢⱼ;Eⱼₖ;Eᵢₗ]
   
   Le volume sera égal au tier de la somme des flux du volume au travers des faces.
   
 * @param faces Une liste de faces (voir la structure @face) contenant pour 
                chacun les sommets
 * @return Le volume du polytope 
 */
E_Float compute_volume_ngon( const std::vector<face>& faces );
}


#endif
