/**
 * @brief Calcul du rotationnel de la vitesse pour différents types de cellules
 *
 * On calcul le rotationnel moyenné sur une cellule à l'aide du théorème du rotationnel :
 #     ⟶ →          →    →     →                                         
 # ∭ᵥ rot 𝑣 dV = -∬ₛ 𝑣 ∧ d𝑠 où d𝑠 est la normale dirigée vers l'extérieur.
 *
 * Puisque la vitesse est constante par cellule, on en déduit le rotationnel par :s
 *
 #  ⟶ →          n      →     →                                       
 # rot 𝑣  = -(1/V) ∑ ∬  (vₖ ∧ dSₖ)                                      
 #                k=1 Sₖ                                                
 * où V est le volume de la cellule considérée et Sₖ la surface de la kième interface de la cellule (qui contient n interfaces).
 *
 * Le rotationnel est donc considéré comme constant par cellule...
 */
#ifndef _K_POST_ROTATIONAL_HPP_
#define _K_POST_ROTATIONAL_HPP_
#include <vector>
#include "point3d.hpp"
#include "vector3d.hpp"
#include "face.hpp"
#include "Memory/vector_view.hpp"

namespace K_POST
{
    /**
     * @brief Calcule le rotationnel de la vitesse dans un tétraèdre
     * @details Calcule le rotationnel de la vitesse dans un tétraèdre
     * 
     * Soit T=(p₁, p₂, p₃, p₄) un tétraèdre défini dans l'espace (un élément de maillage par exemple)
     * avec des vitesses définies en chaque sommet du tétraèdre. On veut calculer le rotationnel de  
     * la vitesse dans le tétraèdre (supposé constant dans toute la cellule)
     * 
     * Le tétraèdre peut se définir par la fonction de passage du tétraèdre de référence au tétraèdre
     * physique comme :
     * 
     # p(ψ, η, ζ) = p₁ + ψ.E₁₂ + η.E₁₃ +  ζ.E₁₄; avec 0 ≤ ψ+η+ζ ≤ 1                                  
     #                   ⎛           ⎞ ⎛ψ⎞          ⎛ψ⎞                                                
     #            = p₁ + ⎜E₁₂;E₁₃;E₁₄⎟.⎜η⎟ = p₁ + F.⎜η⎟                                                
     #                   ⎝           ⎠ ⎝ζ⎠          ⎝ζ⎠                                                
     * 
     * En ce qui concerne le volume du tétraèdre, voir la documentation dans volume.hpp
     * 
     *  La vitesse étant définie uniquement sur les sommets du tétraèdre, on peut l'exprimer sous la forme
     *  d'une fonction linéaire en ψ,η etζ :
     * 
     # v(ψ, η, ζ) = v₁ + ψ.(v₂-v₁) + η.(v₃-v₁) +  ζ.(v₄-v₁); avec 0 ≤ ψ+η+ζ ≤ 1                      
     #                   ⎛              ⎞ ⎛ψ⎞          ⎛ψ⎞                                             
     #            = v₁ + ⎜Δv₁₂;Δv₁₃;Δv₁₄⎟.⎜η⎟ = v₁ + B.⎜η⎟                                             
     #                   ⎝              ⎠ ⎝ζ⎠          ⎝ζ⎠                                             
     *
     * Calculons les intégrales sur chaque face du tétraèdre :
     * 
     * Pour la première face, le triangle T₁₃₂ = (p₁, p₃, p₂), on a ζ = 0 et sa normale est (contenant la surface du triangle ) 
     * 
     #      n₁₃₂ = E₁₃ ∧ E₁₂                                                                 
     *
     * L'intégrale vaut donc :
     * 
     #    →  →            1 1-η                                                               
     # ∬ v∧n₁₃₂ dT₁₃₂ = ∫ ∫   [v₁+ψ.(v₂-v₁) + η.(v₃-v₁)]∧[E₁₃ ∧ E₁₂] dψ dη                  
     #  T₁₃₂            0 0                                                                   
     #                 = (1/2)(v₁∧E₁₃∧E₁₂) + (1/6)((v₂-v₁)^E₁₃∧E₁₂) + (1/6)((v₃-v₁)^E₁₃∧E₁₂)
     #                 = (1/6)((v₁+v₂+v₃)∧E₁₃∧E₁₂)                                           
     *
     * Pour la seconde face, le triangle T₁₂₄ = (p₁,p₂,p₄), on a η = 0 et sa normale   
     * 
     #      n₁₂₄ = E₁₂∧E₁₄                                                                    
     *
     * L'intégrale vaut donc :
     * 
     #    →  →            1 1-η                                                               
     # ∬ v∧n₁₂₄ dT₁₂₄ = ∫ ∫   [v₁+ψ.(v₂-v₁) + ζ.(v₄-v₁)]∧[E₁₂∧E₁₄] dψ dζ                    
     #  T₁₂₄            0 0                                                                   
     #                 = (1/2)(v₁∧E₁₂∧E₁₄) + (1/6)((v₂-v₁)^E₁₂∧E₁₄) + (1/6)((v₄-v₁)^E₁₂∧E₁₄)
     #                 = (1/6)((v₁+v₂+v₄)∧E₁₂∧E₁₄)                                           
     *
     * Pour la troisième face, le triangle T₂₃₄ = (p₂, p₃, p₄), on a ψ + η + ζ = 1 soit ζ = 1 - ψ - η et la normale :
     * 
     #      n₂₃₄ = E₂₃∧E₂₄                                                                    
     *
     * L'intégrale vaut donc :
     *
     #    →  →            1 1-η                                                               
     # ∬ v∧n₂₃₄ dT₂₃₄ = ∫ ∫   [v₁+ψ.(v₂-v₁)+ η.(v₃-v₁)+(1-ψ-η).(v₄-v₁)]∧[E₂₃∧E₂₄] dψ dη     
     #  T₂₃₄            0 0                                                                   
     #                 = (1/6)((v₂+v₃+v₄)∧E₂₃∧E₂₄)                                           
     *
     * Pour la dernière face, le triangle T₃₁₄ = (p₃,p₁,p₄), on a ψ = 0 et la normale :
     * 
     #         n₃₁₄ = E₃₁∧E₃₄                                                                 
     *
     * L'intégrale vaut donc :
     * 
     #    →  →            1 1-η                                                               
     # ∬ v∧n₃₁₄ dT₃₁₄ = ∫ ∫   [v₁+η.(v₃-v₁)+ζ.(v₄-v₁)]∧[E₃₁∧E₃₄] dη dζ                    
     #  T₃₁₄            0 0                                                                   
     #                 = (1/6)((v₁+v₃+v₄)∧E₃₁∧E₃₄)                                           
     *
     ? Le rotationnel vaut donc :
     #  →                                                                                     
     # rot(v) = (-1/(6V)).[(v₁+v₂+v₃)∧E₁₃∧E₁₂ + (v₁+v₂+v₄)∧E₁₂∧E₁₄ + (v₂+v₃+v₄)∧E₂₃∧E₂₄ +  
     #                     (v₁+v₃+v₄)∧E₃₁∧E₃₄]                                               
     #        = (-1/(6V)).[v₁∧E₃₂∧E₃₄ + v₂∧E₁₃∧E₃₄ + v₃∧E₂₁∧E₁₄ + v₄∧E₁₂∧E₁₃]             
     *                          
     * @param p1 Premier    sommet p₁ du tétraèdre
     * @param p2 Second     sommet p₂ du tétraèdre
     * @param p3 Troisième  sommet p₃ du tétraèdre
     * @param p4 Quatrième  sommet p₄ du tétraèdre
     * @param v1 Vitesse au sommet p₁ 
     * @param v2 Vitesse au sommet p₂
     * @param v3 Vitesse au sommet p₃
     * @param v4 Vitesse au sommet p₄
     * @return Le rotationnel constant au sein du tétraèdre
     */
    vector3d compute_rotational_on_tetrahedra( const point3d& p1, const point3d& p2,
                                               const point3d& p3, const point3d& p4,
                                               const vector3d& v1, const vector3d& v2,
                                               const vector3d& v3, const vector3d& v4 );
    /**
     * @brief Calcul du rotationnel de la vitesse sur une pyramide
     * @details Calcul du rotationel de la vitese sur une pyramide dont la vitesse est définie aux sommets
     * 
     * Soit une pyramide P = (p₁,p₂,p₃,p₄,p₅) défini par ses cinq sommets, où (p₁,p₂,p₃,p₄) définissent la base
     * de la pyramide. On définit la pyramide par la fonction transformant la pyramide de référence en la pyramide
     * considérée, comme dans :
     * 
     * @ref Higher-Order Finite Elements for Hybrid Meshes Using New Nodal Pyramidal Elements; Morgane Bergot, Gary Cohen 
     * @ref et Marc Duruflé, https://hal.archives-ouvertes.fr/hal-00454261, Février 2010
     * 
     # F(ψ,η,ζ) = (1/4)[ψ.(E₁₂-E₃₄) + η.(E₁₃+E₂₄) + ζ.(E₁₅+E₂₅+E₃₅+E₄₅) + ψη/(1-ζ).(E₂₃-E₁₄)]                
     *
     ! Remarque : on a choisit de faire un changement de repère en prenant pour origine le barycentre de la base
     !            de la pyramide.
     *
     * Pour le calcul du volume de la pyramide voir :
     * @ref volume.hpp
     * En ce qui concerne la vitesse définie en chaque sommet de la pyramide, on peut l'exprimer par une fonction
     * l'interpolant dans la pyramide :
     * 
     # v(ψ,η,ζ) = (1/4)[v₁+v₂+v₃+v₄ + ψ.(v₂+v₃-v₁-v₄) + η.(v₃+v₄-v₁-v₂) +      
     #                  ζ.(4.v₅-v₁-v₂-v₃-v₄) + ψη/(1-ζ).(v₃+v₁-v₂-v₄)]         
     * 
     * Calculons les intégrales sur chaque face de la pyramide :
     * 
     * Pour le quadrangle Q₁₄₃₂=(p₁, p₄, p₃, p₂), on impose ζ=0 et la normale est définie comme :
     * 
     #      n₁₄₃₂ = (1/16)[E₁₃+E₂₄ + ψ.(E₂₃+E₄₁)] ∧ [E₁₂+E₄₃ + η.(E₂₃+E₄₁)]   
     #            = (1/8)[E₂₄∧E₁₃ + ψ.E₄₃∧E₁₂ + η.E₁₄∧E₂₃]                   
     *
     * On a donc :
     * 
     #    →   →                 1  1                                           
     # ∬ v∧n₁₄₃₂ dQ₁₄₃₂ =(1/32)∫ ∫[v₁+v₂+v₃+v₄+ψ(v₂+v₃-v₁-v₄)+η(v₃+v₄-v₁-v₂) +
     #  Q₁₄₃₂                 -1 -1 ψη(v₃+v₁-v₂-v₄)]∧[E₂₄∧E₁₃ + ψ.E₄₃∧E₁₂   + 
     #                                                  η.E₁₄∧E₂₃] dψ dη       
     #                 = (1/8){v₁∧[E₂₄∧E₁₃+(1/3)E₃₄∧E₁₂+(1/3)E₄₁∧E₂₃)        +
     #                         v₂∧[E₂₄∧E₁₃+(1/3)E₄₃∧E₁₂+(1/3)E₄₁∧E₂₃]        +
     #                         v₃∧[E₂₄∧E₁₃+(1/3)E₄₃∧E₁₂+(1/3)E₁₄∧E₂₃]        +
     #                         v₄∧[E₂₄∧E₁₃+(1/3)E₃₄∧E₁₂+(1/3)E₁₄∧E₂₃]}  
     *
     * Pour le triangle T₁₂₅, on impose η=-1+ζ, et la normale vaut :
     * 
     #      n₁₂₅ = (1/2)E₁₂∧E₁₅                                                
     *
     * On a donc :
     * 
     #    →  →               1  1-ζ                                             
     # ∬ v∧n₁₂₅ dT₁₂₅=(1/4)∫  ∫ [v₁+v₂+ψ(v₂-v₁)+ζ(2v₅-v₁-v₂)]∧[E₁₂∧E₁₅] dψ dζ
     #  T₁₂₅               0 -1+ζ                                              
     #                 = (1/6)(v₁∧E₁₂∧E₁₅+v₂∧E₁₂∧E₁₅+v₅∧E₁₂∧E₁₅)             
     * 
     * Pour le triangle T₂₃₅, on impose ψ= 1-ζ et la normale vaut :
     * 
     #      n₂₃₅ = (1/2)E₂₃∧E₂₅                                                
     * 
     * On a donc :
     * 
     #    →   →              1  1-ζ                                            
     # ∬ v∧n₂₃₅ dT₂₃₅=(1/4)∫  ∫ [v₂+v₃+η(v₃-v₂)+ζ(2v₅-v₂-v₃)]∧[E₂₃∧E₂₅] dη dζ
     #  T₂₃₅               0 -1+ζ                                              
     #                =(1/6)(v₂∧E₂₃∧E₂₅+v₃∧E₂₃∧E₂₅+v₅∧E₂₃∧E₂₅)              
     * 
     * Pour le triangle T₃₄₅, on impose η=+1-ζ et la normale vaut :
     * 
     #     n₃₄₅ = (1/2)E₃₄∧E₃₅                                                 
     * 
     * On a donc :
     * 
     #    →  →               1  1-ζ                                            
     # ∬ v∧n₃₄₅ dT₃₄₅=(1/4)∫  ∫ [v₃+v₄+ψ(v₃-v₄)+ζ(2v₅-v₃-v₄)]∧[E₃₄∧E₃₅] dψ dζ
     #  T₃₄₅               0 -1+ζ                                              
     #                =(1/6)(v₃∧E₃₄∧E₃₅+v₄∧E₃₄∧E₃₅+v₅∧E₃₄∧E₃₅)
     * 
     * Pour le triangle T₄₁₅, on impose ψ=-1+ζ et la normale vaut :
     * 
     #     n₄₁₅ = (1/2)E₄₁∧E₄₅                                                 
     *         
     * On a donc :
     * 
     #    →  →                1  1-ζ                                            
     # ∬ v∧n₄₁₅.dT₄₁₅ =(1/4)∫  ∫ [v₁+v₄+η(v₄-v₁)+ζ(2v₅-v₁-v₄)]∧[E₄₁∧E₄₅] dη dζ
     #  T₄₁₅                0 -1+ζ                                              
     #                 =(1/6)(v₁∧E₄₁∧E₄₅+v₄∧E₄₁∧E₄₅+v₅∧E₄₁∧E₄₅)              
     * 
     ? On en déduit le rotationnel comme :
     * 
     #   →                                                                     
     # rot(v) = (-1/(24V)).[3(v₁+v₂+v₃+v₄)∧E₂₄∧E₁₃ + (v₁-v₂-v₃+v₄)∧E₃₄∧E₁₂  +
     #                       (v₁+v₂-v₃-v₄)∧E₄₁∧E₂₃ + 4(v₁+v₂+v₅)∧E₁₂∧E₁₅    +
     #                       4(v₂+v₃+v₅)∧E₂₃∧E₂₅ + 4(v₃+v₄+v₅)∧E₃₄∧E₃₅      +
     #                       4(v₁+v₄+v₅)∧E₄₁∧E₄₅]
     * 
     * où V est le volume de la pyramide
     * 
     * @param p1 Premier   sommet p₁ de la pyramide
     * @param p2 Second    sommet p₂ de la pyramide
     * @param p3 Troisième sommet p₃ de la pyramide
     * @param p4 Quatrième sommet p₄ de la pyramide
     * @param p5 Cinquième sommet p₅ de la pyramide
     * @param v1 Valeur de la vitesse au sommet p₁
     * @param v2 Valeur de la vitesse au sommet p₂
     * @param v3 Valeur de la vitesse au sommet p₃
     * @param v4 Valeur de la vitesse au sommet p₄
     * @param v5 Valeur de la vitesse au sommet p₅
     * @return La valeur du rotationnel (constant) dans la pyramide
     */
    vector3d compute_rotational_on_pyramid( const point3d& p1, const point3d& p2,
                                            const point3d& p3, const point3d& p4, const point3d& p5,
                                            const vector3d& v1, const vector3d& v2,
                                            const vector3d& v3, const vector3d& v4, const vector3d& v5 );
    /**
     * @brief Calcul du rotationnel de la vitesse dans une cellule pentaèdrique
     * @details Calcul le rotationnel de la vitesse dans une cellule pentaèdrique avec les vitesses définies aux sommets.
     * 
     * Soit un pentaèdre P=(p₁,p₂,p₃,p₄,p₅,p₆) dont les sommets sont ordonnés selon la norme CGNS.
     * Pour la simplification des calculs, on peut supposer que le sommet p₁ est à l'origine du repère,
     * quitte à faire une translation du repère cartésien pour le besoin.
     * 
     * Le pentaèdre peut être décrit à l'aide d'une fonction transformant le pentaèdre de référence dans le pentaèdre
     * P à l'aide de trois paramètres ψ,η,ζ :
     * 
     #         ϕ : [0;1]³ ----> ℝ³                                                                         
     #             ψ,η,ζ  ----> ψ.E₁₂ + η.E₁₃ + ζ.E₁₄ + ψζ.(E₄₅ + E₂₁) + ηζ.(E₄₆ + E₃₁) avec 0 ≤ ψ + η ≤ 1 
     *             
     * où Eᵢⱼ définit le vecteur ayant pour origine le sommet pᵢ et pour extrémité le sommet pⱼ.
     * 
     * La vitesse étant définie en chaque sommet du pentaèdre, on peut interpoler la vitesse au sein de pentaèdre
     * à l'aide de la fonction polynomiale :
     * 
     # v(ψ,η,ζ) = v₁ + ψ(v₂-v₁) + η(v₃-v₁) + ζ(v₄-v₁) + ψζ(v₁+v₅-v₂-v₄) + ηζ(v₁+v₆-v₃-v₄)                 
     * 
     * Nous allons calculer le rotationnel en utilisant le théorème du rotationnel.
     * Calculons les intégrales sur chaque face du pentaèdre :
     * 
     * Pour la face Q₁₂₅₄ = (p₁,p₂,p₅,p₄), on impose η = 0 et la normale s'exprime comme :
     * 
     #    n₁₂₅₄(ψ,ζ) = (E₁₂ + ζ.(E₄₅ + E₂₁))∧(E₁₄ + ψ.(E₄₅ + E₂₁))             
     #               = E₁₂∧E₁₄ + ψ.E₁₂∧E₄₅ + ζ.E₂₅∧E₁₄                        
     * 
     * L'intégrale vaut donc :
     * 
     #    →   →            1  1                                                
     # ∬ v∧n₁₂₅₄ dQ₁₂₅₄ = ∫ ∫[v₁+ψ(v₂-v₁)+ζ(v₄-v₁)+ ψζ(v₁+v₅-v₂-v₄)]∧[E₁₂∧E₁₄ +
     #                    0 0                    ψ.E₁₂∧E₄₅ + ζ.E₂₅∧E₁₄] dψ dζ 
     #                   = (1/12) [ v₁∧E₄₂∧(E₁₅+E₁₄) + v₂∧(E₁₂+E₄₂)∧E₁₅       +
     #                              v₄∧E₁₅∧(E₂₄+E₁₄) + v₅∧(E₁₂+E₁₅)∧E₂₄ ]     
     * 
     * Pour la face Q₂₃₆₅ = (p₂,p₃,p₆,p₅), on impose η = 1-ψ et la normale s'exprime comme :
     * 
     #    n₂₃₆₅(ψ,ζ) = (E₃₆ + ψ.(E₆₅ + E₂₃))∧(E₃₂ + ζ.(E₆₅ + E₂₃))             
     #               =  E₃₆∧E₃₂ + ψ.E₆₅∧E₃₂ + ζ.E₃₆∧E₂₅                       
     *    
     * L'intégrale vaut donc 
     * 
     #    →   →            1  1                                                
     # ∬ v∧n₂₃₆₅ dQ₂₃₆₅ = ∫ ∫[v₃+ψ(v₂-v₃)+ζ(v₆-v₃)+ψζ(v₅-v₆+v₃-v₂)]∧[E₃₆∧E₃₂ +
     #                    0 0                    ψ.E₆₅∧E₃₂ + ζ.E₃₆∧E₂₅] dψ dζ 
     #                   = (1/12) [ v₃∧(E₃₅+E₃₂)∧E₆₂ + v₂∧E₃₅∧(E₆₂+E₃₂)       +
     #                              v₅∧E₂₆∧(E₂₅+E₃₅) + v₆∧(E₃₆+E₂₆)∧E₃₅ ]     
     * 
     * Pour la face Q₃₁₄₆ = (p₃,p₁,p₄,p₆), on impose ψ=0 et la normale s'exprime comme :
     * 
     #    n₃₁₄₆(η,ζ) = (E₁₄ + η.(E₄₆ + E₃₁))∧(E₁₃ + ζ.(E₄₆ + E₃₁))             
     #               = E₁₄∧E₁₃ + η.E₄₆∧E₁₃ + ζ.E₁₄∧E₃₆                        
     * 
     * L'intégrale vaut donc :
     * 
     #    →   →            1  1                                                
     # ∬ v∧n₃₁₄₆ dQ₃₁₄₆ = ∫ ∫[v₁+η(v₃-v₁)+ζ(v₄-v₁)+ηζ(v₁+v₆-v₃-v₄)]∧[E₁₄∧E₁₃ +
     #                    0 0                    η.E₄₆∧E₁₃ + ζ.E₁₄∧E₃₆] dη dζ 
     #                   = (1/12) [ v₁∧E₃₄∧(E₁₆+E₁₃) + v₃∧(E₃₄+E₃₁)∧E₁₆       +
     #                              v₄∧(E₃₄+E₁₄)∧E₁₆ + v₆∧E₃₄∧(E₃₆+E₁₆) ]     
     * 
     * Pour la face T₁₃₂ = (p₁,p₃,p₂), on impose ζ=0 et la normale vaut :
     * 
     #    n₁₃₂ = E₁₃∧E₁₂                                                       
     * 
     * et l'intégrale vaut donc :
     * 
     #    →   →           1 1-η                                                
     # ∬ v∧n₁₃₂ dT₁₃₂ = ∫ ∫ [v₁ + ψ(v₂-v₁) + η(v₃-v₁) ]∧[E₁₃∧E₁₂] dψ dη      
     #                   0 0                                                   
     #                 = (1/6)(v₁+v₂+v₃)∧E₁₃∧E₁₂                              
     * 
     * Pour la face T₄₅₆ = (p₄,p₅,p₆), on impose ζ=1 et la normale vaut :
     * 
     #    n₄₅₆ = E₄₅∧E₄₆                                                       
     * 
     #    →   →           1 1-η                                                
     # ∬ v∧n₁₃₂ dT₁₃₂ = ∫ ∫ [v₄ + ψ(v₅-v₄) + η(v₆-v₄) ]∧[E₄₅∧E₄₆] dψ dη      
     #                   0 0                                                   
     #                 = (1/6)(v₄+v₅+v₆)∧E₄₅∧E₄₆                              
     * 
     * En sommant les cinq intégrales calculées ci--dessus et en multipliant par (-1/V),
     * on trouve bien le rotationnel recherché.
     * 
     * @param p1 Premier   sommet du pentaèdre
     * @param p2 Second    sommet du pentaèdre
     * @param p3 Troisième sommet du pentaèdre
     * @param p4 Quatrième sommet du pentaèdre
     * @param p5 Cinquième sommet du pentaèdre
     * @param p6 Sixième   sommet du pentaèdre
     * @param v1 Valeur de la vitesse au point p₁
     * @param v2 Valeur de la vitesse au point p₂
     * @param v3 Valeur de la vitesse au point p₃
     * @param v4 Valeur de la vitesse au point p₄
     * @param v5 Valeur de la vitesse au point p₅
     * @param v6 Valeur de la vitesse au point p₆
     * @return Le vecteur rotationnel obtenu
     */
    vector3d compute_rotational_on_pentahedra( const point3d& p1, const point3d& p2, const point3d& p3, 
                                               const point3d& p4, const point3d& p5, const point3d& p6,
                                               const vector3d& v1, const vector3d& v2, const vector3d& v3,
                                               const vector3d& v4, const vector3d& v5, const vector3d& v6);
    /**
     * @brief Calcul le rotationnel de la vitesse sur un hexaèdre
     * @details Calcul le rotationnel de la vitesse définie aux sommets d'un hexaèdre
     * 
     * Soit H = (p₁,p₂,p₃,p₄,p₅,p₆,p₇,p₈) un hexaèdre. On suppose, quitte à faire une translation du repère,
     * que p₁ se situe à l'origine du repère. On peut alors trouver une fonction permettant de passer de l'hexaèdre de
     * référence à H à l'aide d'une fonction trilinéaire :
     * 
     ?   p₈____________p₇       
     ?   ╱⁞           ╱│        
     ?  ╱ ⁞          ╱ │
     ? p₅___________p₆ │        ϕ : [0;1]³ ----> ℝ³
     ? │  ⁞         │  │            ψ,η,ζ  ----> ψ.E₁₂ + η.E₁₄ + ζ.E₁₅ + ψη.(E₄₃+E₂₁) + ψζ.(E₅₆+E₂₁) +
     ? │  p₄……………………│……p₃                        ηζ.(E₅₈+E₄₁) + ψηζ.(E₈₇+E₃₄+E₆₅+E₁₂)
     ? │  /         │ ╱
     ? │ /          │╱          avec Eᵢⱼ le vecteur ayant pour origine le sommet pᵢ et pour extrémité le sommet pⱼ.
     ? p₁___________p₂
     * 
     * Puisque le champs de vitesse est défini par rapport aux valeurs de la vitesse aux sommets de H, on peut l'exprimer
     * également à l'aide d'une fonction trilinéaire :
     * 
     # v(ψ,η,ζ) = v₁ + ψ.(v₂-v₁) + η.(v₄-v₁) + ζ.(v₅-v₁) + ψη.(v₃-v₄+v₁-v₂)    +
     #      ψζ.(v₆-v₅+v₁-v₂) + ηζ.(v₈-v₅+v₁-v₄) + ψηζ.(v₇-v₈+v₄-v₃+v₅-v₆+v₂-v₁)
     * 
     * Calculons les intégrales sur chaque face de l'hexaèdre pour ensuite utiliser le théorème du rotationnel.
     * On utilise pour cela la formule générale trouvée pour les n-gons dans le cadre de face quadrangulaire
     *   
     * @param p1 Premier sommet de l'hexaèdre                           
     * @param p2 second sommet de l'hexaèdre                            
     * @param p3 Troisième sommet de l'hexaèdre                         
     * @param p4 Quatrième sommet de l'hexaèdre                         
     * @param p5 Cinquième sommet de l'hexaèdre                         
     * @param p6 Sixième sommet de l'hexaèdre                           
     * @param p7 Septième sommet de l'hexaèdre                          
     * @param p8 Huitième sommet de l'hexaèdre                          
     * @param v1 Vecteur vitesse au premier sommet de l'hexaèdre        
     * @param v2 Vecteur vitesse au second sommet de l'hexaèdre         
     * @param v3 Vecteur vitesse au troisième sommet de l'hexaèdre      
     * @param v4 Vecteur vitesse au quatrième sommet de l'hexaèdre      
     * @param v5 Vecteur vitesse au cinquième sommet de l'hexaèdre      
     * @param v6 Vecteur vitesse au sixième sommet de l'hexaèdre        
     * @param v7 Vecteur vitesse au septième sommet de l'hexaèdre       
     * @param v8 Vecteur vitesse au huitième sommet de l'hexaèdre       
     * @return Retourne le rotationel moyen au sein de l'hexaèdre       
     */
    vector3d compute_rotational_on_hexaedra( const point3d& p1, const point3d& p2, const point3d& p3, const point3d& p4,
                                             const point3d& p5, const point3d& p6, const point3d& p7, const point3d& p8,
                                             const vector3d& v1, const vector3d& v2, const vector3d& v3, const vector3d & v4,
                                             const vector3d& v5, const vector3d& v6, const vector3d& v7, const vector3d& v8 );

    /**
     * @brief Calcul le rotationnel de la vitesse sur un polygone
     * @details Calcul le rotationnel de la vitesse sur un polygone quelconque à face triangulaire ou quadrangulaire
     * 
     * L'idée ici est de calculer chaque intégrale de surface en paramétrisant localement la face dans un repère
     * absolu.
     * 
     * Ainsi, pour une face triangulaire Fᵢⱼₖ = (pᵢ,pⱼ,pₖ), la paramétrisation de la face est :
     * 
     # ϕᵢⱼₖ : [0;1]² ----> ℝ³ avec 0 ≤ ψ+η ≤ 1                                 
     #         ψ,η  ----> pᵢ + ψ.Eᵢⱼ + η.Eᵢₖ                                   
     *         
     * et la normale (orientée vers l'extérieur du polygone) vaut :
     * 
     #      nᵢⱼₖ = Eᵢⱼ∧Eᵢₖ                                                     
     *
     * De même sur Fᵢⱼₖ, la vitesse peut être exprimer par une fonction linéaire :
     * 
     # v(ψ,η) = vᵢ + ψ.(vⱼ-vᵢ) + η.(vₖ-vᵢ)                                     
     *
     * La valeurs de l'intégrale est alors :
     * 
     #    →  →              1  1-η                                             
     # ∬ v∧nᵢⱼₖ dTᵢⱼₖ=(1/4)∫  ∫ [vᵢ + ψ.(vⱼ-vᵢ) + η.(vₖ-vᵢ)]∧[Eᵢⱼ∧Eᵢₖ] dψ dη   
     #  Tᵢⱼₖ               0  0                                                
     #               = (1/6)(vᵢ+vⱼvₖ)∧Eᵢⱼ∧Eᵢₖ                                  
     *
     * Pour une face quadrangulaire Fᵢⱼₖₗ = (pᵢ,pⱼ,pₖ,pₗ), la paramétrisation de la face est :
     * 
     # ϕᵢⱼₖₗ : [0;1]² ----> ℝ³                                                  
     #          ψ,η  ----> pᵢ + ψ.Eᵢⱼ + η.Eᵢₗ + ψη.(Eⱼₖ + Eₗᵢ)                   
     *         
     * et la normale est en fonction de ψ et η :
     * 
     #      nᵢⱼₖₗ = (Eᵢⱼ + η.(Eⱼₖ + Eₗᵢ))∧(Eᵢₗ + ψ.(Eⱼₖ + Eₗᵢ))                    
     #           = Eᵢⱼ∧Eᵢₗ + ψ.Eᵢⱼ∧Eₗₖ + η.Eⱼₖ∧Eᵢₗ                               
     *           
     * La valeurs de l'intégrale est alors :
     * 
     #    →   →           1 1                                                 
     # ∬ v∧nᵢⱼₖₗ dFᵢⱼₖₗ = ∫ ∫  [vᵢ+ψ(vⱼ-vᵢ)+η(vₗ-vᵢ)+ ψη(vᵢ+vₖ-vⱼ-vₗ)]∧[Eᵢⱼ∧Eᵢₗ   +
     #                   0 0                         ψ.Eᵢⱼ∧Eₗₖ + η.Eⱼₖ∧Eᵢₗ] dψ dζ
     #                = (1/12)[vᵢ∧Eₗⱼ∧(Eᵢₖ+Eᵢₗ) + vⱼ∧Eᵢₖ∧(Eⱼₗ+Eⱼᵢ)              +
     #                         vₖ∧Eₗⱼ∧(Eₗₖ+Eᵢₖ)  + vₗ∧Eᵢₖ∧(Eⱼₗ+Eᵢₗ) ]              
     *  
     *  Il suffit alors de sommer les intégrales  calculées  par face et
     *  de diviser par le volume du polygone
     @ ref Pour le volume, voir volume.hpp pour la méthode de calcul  du
     -     volume du polygone                                           
     * 
     * @param faces Le tableau des faces du polygone                    
     @ ref voir  face.hpp  pour  une  plus  ample   explication  sur  la
     -     structure et les méthodes de face                            
     * @return Le rotationnel (constant) calculé dans le polygône       
     */
    vector3d 
    compute_rotational_on_ngon( const std::vector<face>& faces, 
                                const std::array<K_MEMORY::vector_view<const E_Float>,3>& velocity );
}

#endif
