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
/**
 * @defgroup   TRIANGULATED_POLYHEDRON triangulated polyhedron
 *
 * @brief      This file implements triangulated polyhedron.
 *
 * Test robuste d'appartenance d'un point à un polyèdre
 * d'après un algorithme en python de Mark Dickinson :
 *              https://github.com/mdickinson/polyhedron
 *
 * Soit une surface fermée, orientée dans ℝ³ décrite par un maillage triangulaire.
 * Le code ci-dessous donne un algorithme robuste pour déterminer si un point
 * donné est dans, à la surface ou en dehors de la surface. L'algorithme devrait
 * donner des résultats corrects même dans les cas dégénérés, et s'applique également
 * à des polyèdres non connexes, des surfaces non simplement connexes, etc.
 * La surface n'a pas besoin d'être convexe, simple, connexe ou simplement connexe.
 * 
 * Plus précisément, la méthode proposée ici est basée sur le calcul de l'indice
 * (en analyse complexe, winding-number en anglais) d'une surface orientée fermée S
 * autour d'un point O qui n'est pas sur S. Grossièrement, l'indice d'une surface
 * orientée, fermée S autour d'un point O n'appartenant pas à S est le nombre de fois
 * que la surface enferme ce point; pour une surface simple, dont les normales sont orientées
 * à l'extérieur (par exemple, un polyèdre convexe), l'indice sera de un pour les points contenus
 * dans l'intérieur, et zéro pour les points extérieurs.
 * 
 * Pour une définition plus précise de l'indice, il faut faire appel à la topologie
 * algégrique : notre surface orientée est représentée par une collection de données
 * combinatoires définissant des sommets abstraits, des arêtes et des triangles, avec
 * une application des sommets vers ℝ³. Les données combinatoires décrivent
 * un complexe simplicial C, et en supposant que O n'est pas sur la surface, l'application
 * des sommets à ℝ³ nous donne une application continue de la réalisation géométrique de C
 * vers ℝ³ - {O}. Cela induit une application sur le second groupe d'homologie :
 * 
 *     ℍ²(C,Z) → ℍ²(ℝ³ - {O}, Z)
 *     
 * et en prenant l'orientation direct usuelle dans ℝ³, nous identifions
 * ℍ²(ℝ³ - {O}, Z) avec Z. L'image de [S] par cette application nous donne notre indice.
 * En particuliers, la bonne définition de notre indice ne dépend pas des propriétés
 * topologiques de l'enveloppe : ce n'est pas grave si la surface est auto-intersectante ou
 * possède des triangles dégénérés. La seule condition est que O n'appartienne pas à S.
 * 
 * Algorithme utilisé :
 * --------------------
 * 
 * L'algorithme est basé sur la méthode classique de lancer de rayon : on prend une
 * ligne verticale L passant par O et on compte le nombre d'intersections de cette ligne
 * avec les triangles de la surface, en prenant soin de garder trace des orientations au
 * fur et à mesure.. Ignorons les cas particuliers pour le moment et supposons que :
 * 
 *  ① O n'appartient pas à la surface, et
 *  ② Pour chaque triangle T (formant un sous ensemble fermé de ℝ³) intersectant notre
 *    ligne verticale L, L s'intersecte avec l'intérieur de T en exactement un point Q
 *  
 *  Il y a alors quatre possibilités pour chacun des triangles T :
 *  ⑴ T se trouve "au dessus" de O et est orienté "vers le haut" ("s'éloignant" de O).
 *  ⑵ T se trouve "au dessus" de O et est orienté "vers le bas"  ("se dirigeant" vers O).
 *  ⑶ T se trouve "au dessous" de O et est orienté "vers le bas" ("s'éloignant" de O).
 *  ⑷ T se trouve "en dessous" de O et est orienté "vers le haut"("se dirigeant" vers O).
 *  
 *  Notons respectivement N₁, N₂, N₃ et N₄ le nombre de triangles satisfaisant les conditions
 *  ⑴, ⑵, ⑶ et ⑷. Puisque notre surface est fermée, ces nombres ne sont pas indépendants;
 *  ils satisfont la relation :
 *  
 *      N₁ + N₄ =  N₂ + N₃
 *  
 *  Ainsi, le nombre de triangles dirigés "vers le haut" doit être égal au nombre de triangles
 *  dirigés "vers le bas". L'indice w est alors donné par la relation :
 *  
 *      w =  N₁ - N₂ = N₃ - N₄
 *  
 *  Dans le code donné ci--dessous, on calcule simplement 2w = (N₁ + N₃) - (N₂ + N₄), tel
 *  que chaque triangle orienté en "s'éloignant" de O contribue pour 1 à 2w tandis que chaque
 *  triangle orienté en "se dirigeant" vers O contribue pour -1.
 *  
 *  Pour rendre l'algorithme robuste :
 *  ................................
 *  
 *  Nous allons décrire comment améliorer l'algorithme basique décrit au dessus pour
 *  prendre en compte :
 *  
 *     - Le traitement correction des cas dégénérés (triangle vertical, cas où L s'intersecte
 *       avec une arête ou un sommet directement, etc.)
 *     - Détecter le cas où le point se trouve directement sur la surface.
 * 
 * Il se trouve que pour rendre l'algorithme robuste, tout ce qui est nécessaire est d'être
 * consistant et faire attention en classifiant les sommets, les arêtes et les triangles. Cela
 * se fait de la manière suivante :
 * 
 *     - Chaque sommet de la surface qui n'est pas égal à O est considéré comme *positif* si
 *       ses coordonnées sont lexicographiquement plus grandes que O, et *négatif* sinon.
 *       
 *     - Pour une arête PQ de la surface qui n'est pas colinéaire avec O, on décrit en premier
 *       dans la classification que P est négatif et Q positif, et on étend à des PQ arbitraires.
 *       
 *       Pour P négatif et Q positif, il y a deux cas :
 *       
 *       ⒈ P et Q ont des abcisses x distinctes. Dans ce cas, nous classifions l'arête PQ
 *         par ses intersections avec le plan passant par O et parallèle au plan YZ : l'arête
 *         est *positive* si le point d'intersection est positif, et *négatif* sinon.
 *         
 *       ⒉ P et Q ont la même abcisse x, et dans ce cas, ils doivent avoir une ordonnée y
 *         différente (si ils ont même abcisse et même ordonnée, alors PQ passe par O.) On classifie
 *         par l'intersection de PQ avec la ligne parallèle à l'axe des y passant par O.
 *         
 *         Pour P positif et Q négatif, nous classifions de la même manière mais en inversant le signe.
 *         Pour P et Q ayant même signes, la classification n'est pas utilisée.
 *       
 *         Du point de vue calculatoire, dans le cas ⒈ au dessus, l'ordonnée du point d'intersection est :
 *       
 *            P𝑦 + (Q𝑦 - P𝑦) * (O𝑥-P𝑥)/(Q𝑥-P𝑥)
 *          
 *         et c'est plus grand que Oᵧ si et seulement si :
 *       
 *            (P𝑦-O𝑦) * (Q𝑥 - O𝑥) - (P𝑥 - O𝑥) * (Q𝑦 - O𝑦)
 *          
 *         est positif. Ainsi le signe de l'arête et le signe de l'expression ci-dessus.
 *       
 *         De la même manière, si cette quantité est nulle alors nous avons besoins de regarder la
 *         coordonnée z  de l'intersection et le signe de l'arête est donnée par :
 *       
 *            (P𝑧-O𝑧) * (Q𝑥-O𝑥) - (P𝑥-O𝑥) * (Q𝑧-O𝑧)
 *          
 *         Dans le cas deux,les deux quantités sont toutes les deux nulles, et le signe de l'arête
 *         est le signe de :
 *       
 *            (P𝑧-O𝑧) * (Q𝑦 - O𝑦) - (P𝑦-O𝑦) * (Q𝑧-O𝑧)
 *          
 *         On peut voir cela d'une autre manière : si P, Q et O ne sont pas colinéaires,
 *         alors la matrice :
 *       
 *             ⎛ P𝑥 Q𝑥 O𝑥 ⎞
 *             ⎜ P𝑦 Q𝑦 O𝑦 ⎜
 *             ⎜ P𝑧 Q𝑧 O𝑧 |
 *             ⎝  1  1  1 ⎠
 *           
 *         est de rang 3. Il s'en suit qu'au moins un des trois sous-déterminants 3x3 :
 * 
 *             ⎜ P𝑥 Q𝑥 O𝑥 ⎜ ⎜ P𝑥 Q𝑥 O𝑥 ⎜  ⎜ P𝑦 Q𝑦 O𝑦 ⎜
 *             ⎜ P𝑦 Q𝑦 O𝑦 ⎜ ⎜ P𝑧 Q𝑧 O𝑧 |  ⎜ P𝑧 Q𝑧 O𝑧 | 
 *             ⎜  1  1  1 ⎜ ⎜  1  1  1 ⎜  ⎜  1  1  1 ⎜   
 *           
 *         est non nul. Nous définissons le signe de PQ comme l'opposé du signe du premier
 *         sous-déterminant non nul de cette liste.
 *       
 *     - Chaque triangle PQR de la surface qui n'est pas coplanaire avec O est considéré *positif*
 *       si sa normale pointe en s'éloignant de O, et *négatif" si sa normale pointe vers O.
 *       
 *       De manière calculatoire, le signe du triangle PQR est le signe du déterminant 4x4 :
 *          
 *              | P𝑥 Q𝑥 R𝑥 O𝑥 |
 *              | P𝑦 Q𝑦 R𝑦 O𝑦 ⎜
 *              | P𝑧 Q𝑧 R𝑧 O𝑧 |
 *              |  1  1  1  1 ⎜
 *              
 *       ce qui est équivalent à calculer le déterminant 3x3 :
 *       
 *              | P𝑥-O𝑥 Q𝑥-O𝑥 R𝑥-O𝑥 |
 *              | P𝑦-O𝑦 Q𝑦-O𝑦 R𝑦-O𝑦 ⎜
 *              | P𝑧-O𝑧 Q𝑧-O𝑧 R𝑧-O𝑧  ⎜
 *              
 * Ainsi pour calculer la contribution d'un triangle donné à l'indice total, il faut :
 * 
 * ⓵ Classifier les sommets du triangle. En même temps, on peut vérifier qu'aucun des sommets
 *    est confondu avec O. Si tous les sommets ont le même signe, alors la contribution pour
 *    l'indice du triangle est nulle.
 *    
 * ⓶ Si les sommets n'ont pas tous le même signe, deux des trois arêtes sont connectées
 *    à des sommets de signes opposés. Classifier les deux arêtes (et vérifier en même temps qu'elles
 *    ne passent pas par O). Si les arêtes ont une classification opposée, alors la
 *    contribution du triangle à l'indice est nulle.
 * 
 * ⓷ Si les deux des arêtes ont le même signe : il faut classifier le triangle lui-même. Si le
 *    triangle est positifi, il contribue à ½ de l'indice total; si négatif, il contribue
 *    à -½. En pratique, nous comptons les contributions comme 1 et -1, et divisons par deux
 *    à la fin.
 *    
 * Remarquez qu'une arête entre deux sommets de même signe ne peut pas passer par O, alors il
 * n'y a pas besoin de vérifier la troisième arête à l'étape ⓶. De même, un triangle donc
 * le cycle des arêtes est trivial ne peut pas contenir O à l'intérieur.
 * 
 * POur bien comprendre l'algorithme au dessus, il est conseillé de revenir dans le monde de
 * l'homologie de nouveau. L'homologie de ℝ³ - {O} peut être identifié avec celui
 * de la sphère S² par déformation contractante, et nous pouvons décomposer la sphère comme un CW-complexe
 * consistant de six cellules, comme suit :
 * 
 * Des 0-cellules B et F, où B = (-1, 0, 0) et F = (1, 0, 0)
 * Des 1-cellules L et R, où
 *     L = {(cos(t),sin(t),0) ; -π ≤ t ≤ 0}
 *     R = {(cos(t),sin(t),0) ;  0 ≤ t ≤ π}
 * Des 2-cellules U et D, où U est la moitié supérieure de la sphère (z ≥ 0) et D est la demi-sphère
 * inférieure (z ≤ 0), les deux orientées vers l'extérieur.
 * 
 * L'homologie du CW-complexe est maintenant représentable en termes d'homologie cellulaire :
 * 
 *                  d              d
 *     Z[U] + Z[D] ⟶ Z[L] + Z[R] ⟶ Z[B] + Z[F]
 * 
 * avec les applications frontière :
 *     d[U] = [L] + [R]; d[D] = -[L] - [R]
 *     d[R] = [B] - [F]; d[L] = [F] - [B]
 *     
 * Maintenant l'application d'origine C ⟶ ℝ³ - {O} de la réalisation géométrique
 * du complexe simplicial est homotopique à l'application C ⟶ S² qui retourne :
 * 
 *  - Chaque sommet positif de F et chaque sommet négatif de B
 *  - Chaque arête avec frontière [F] - [B] à L si l'arête est négative, et -R si l'arête
 *    est positive
 *  - Chaque arête avec frontière [B] - [F] à R si l'arête est positive, et -L si l'arête
 *    est négative
 *  - Toutes les autres arêtes à 0
 *  - Chaque triangle avec la frontière [L]+[R] à soit U soit -D selon
 *    que le triangle est positif ou négatif
 *  - Chaque triangle avec la frontière -[L] -[R] à soit D soit -U
 *    selon que le triangle est positif ou négatif
 *  - Tous les autres triangles à 0
 *  
 *  En évaluant chaque triangle de la surface avec l'application et en sommant les
 *  résultats dans la seconde homologie, nous finissons avec :
 *     
 *       (indice)*([U] + [D])
 *  
 *
 * @author     Xavier Juvigny
 * @date       2021
 */
#include <stdexcept>
#include "point3d.hpp"
#include "vector3d.hpp"
#include "triangulated_polyhedron.hpp"

namespace
{
    template<typename T> int sign(const T& val)
    {
        return (val>0) - (val<0);
    }

    /**
     * @brief      Le ou logique de python
     *
     * @param[in]  i     premier entier
     * @param[in]  j     second entier
     *                                               i por  j = ?
     * @return    retourne le ou logique de python : 0 por  0 = 0
     *                                               0 por  1 = 1
     *                                               0 por -1 = -1
     *                                              -1 por -1 = -1
     *                                              -1 por  0 = -1
     *                                              -1 por  1 = -1
     *                                               1 por -1 = 1 (!)
     *                                               1 por  0 = 1
     *                                               1 por  1 = 1
     *            En gros, si i != 0, python prend i sinon il prend j dans son opérateur or
     */
    inline int por(int i, int j)
    {
        return (i!=0 ? i : j);
    }

    /**
     * @brief      Signe du sommet p par rapport à o, comme défini dans le commentaire plus haut
     *
     * @param[in]  p     Le point dont on recherche le signe dans l'homologie S²
     * @param[in]  o     Le point d'origine
     *
     * @return     Le signe du sommet dans l'homologie S²
     */
    int vertex_sign(const K_POST::point3d& p, const K_POST::point3d& o)
    {
        // Attention, il faut utiliser ici l'opération logique or de python (qui me semble un peu étrange...)
        int result = por(por(sign(p.x-o.x),sign(p.y-o.y)), sign(p.z-o.z));
        if (!result) throw std::invalid_argument("Le sommet coïncide avec l'origine !");
        return result;
    }


    /**
     * @brief      Signe de l'arête pq par rapport à o, comme défini au dessus
     *
     * @param[in]  p     Le point origine de l'arête
     * @param[in]  q     Le point extrémité de l'arête
     * @param[in]  o     Le point origine
     *
     * @return     Le signe de l'arête par rapport à o dans l'homologie S²
     */
    int edge_sign( const K_POST::point3d& p, const K_POST::point3d& q, const K_POST::point3d& o)
    {
        K_POST::vector3d op(o,p);
        K_POST::vector3d oq(o,q);
        int result =  por(por(sign(op.y * oq.x - op.x * oq.y),sign(op.z * oq.x - op.x * oq.z)),sign(op.z * oq.y - op.y * oq.z ));
        if (!result) throw std::invalid_argument("Les sommets sont colinéaires avec l'origine");
        return result;
    }

    /**
     * @brief      Retourne le signe du triangle PQR par rapport à o, comme décrit ci-dessus
     *
     * @param[in]  p     Le point P
     * @param[in]  q     Le point Q
     * @param[in]  r     Le point R
     * @param[in]  o     L'origine
     *
     * @return     REtourne le signe du triangle par rapport à l'origine dans l'homologie S²
     */
    int triangle_sign( const K_POST::point3d& p, const K_POST::point3d& q, const K_POST::point3d& r,
                       const K_POST::point3d& o)
    {
        K_POST::vector3d o_p(o,p);
        K_POST::vector3d o_q(o,q);
        K_POST::vector3d o_r(o,r);

        int result = sign( (o_p.x * o_q.y - o_p.y * o_q.x ) * o_r.z +
                           (o_q.x * o_r.y - o_q.y * o_r.x ) * o_p.z +
                           (o_r.x * o_p.y - o_r.y * o_p.x ) * o_q.z );
        if (!result)
            throw std::invalid_argument("Sommets coplanaires avec l'origine");
        return result;
    }

    /**
     * @brief      Retourne la contribution du rectangle v₁v₂v₃ à l'indice
     *
     * @param[in]  v1       Le point v₁
     * @param[in]  v2       Le point v₂
     * @param[in]  v3       Le point v₃
     * @param[in]  origine  Le point origine
     *
     * @return     Retourne la contribution à l'indice de v₁v₂v₃
     */
    int triangle_chain( const K_POST::point3d& v1, const K_POST::point3d& v2, const K_POST::point3d& v3, 
                        const K_POST::point3d& origine )
    {
        int v1_sign = vertex_sign(v1, origine);
        int v2_sign = vertex_sign(v2, origine);
        int v3_sign = vertex_sign(v3, origine);

        int face_boundary = 0;
        if (v1_sign != v2_sign)
            face_boundary += edge_sign(v1, v2, origine);
        if (v2_sign != v3_sign)
            face_boundary += edge_sign(v2, v3, origine);
        if (v3_sign != v1_sign)
            face_boundary += edge_sign(v3, v1, origine);
        if (!face_boundary)
            return 0;

        return triangle_sign(v1, v2, v3, origine);
    }
}
// ================================================================================================
K_POST::triangulated_polyhedron::triangulated_polyhedron(
                                 const std::vector<triangle_type>& trig_faces, 
                                 const std::array<vector_view<const double>,3>& crds) :
    triangular_faces(trig_faces),
    coords(crds)
{}
// ---------------------------------------------------------------------------------------------------    
double
K_POST::triangulated_polyhedron::volume() const
{
    double acc = 0;
    for ( const auto& ind_trig : this->triangular_faces)
    {
        const auto& pos = this->triangle_position(ind_trig);
        // Deux fois l'aire de la projection  sur le plan XY
        double det = ((pos[1][1] - pos[2][1]) * (pos[0][0] - pos[2][0]) -
                      (pos[1][0] - pos[2][0]) * (pos[0][1] - pos[2][1]) );
        // Trois fois la hauteur moyenne :
        double height = pos[0][2] + pos[1][2] + pos[2][2];
        acc += det * height;
    }
    return acc / 6.0;
}

int
K_POST::triangulated_polyhedron::winding_number(const point3d &pt) const
{
    int acc = 0;
    for ( const auto& ind_trig : this->triangular_faces)
    {
        const auto& pos = this->triangle_position(ind_trig);
        acc += triangle_chain(pos[0], pos[1], pos[2], pt);
    }
    return acc/2;
}

// A décommenter si on veut tester les routines en stand-alone (ou sinon rajouter -D__MAIN_TEST__ dans les options de compilation)
// #define __MAIN__TEST__
#if defined(__MAIN__TEST__)
#include <iostream>

using K_POST::triangulated_polyhedron;
using triangle_type = triangulated_polyhedron::triangle_type;
using polyhedron_data_type = std::pair<std::vector<triangle_type>, std::array<std::vector<double>,3>>;


enum {
    IS_OUTSIDE = 0,
    IS_INSIDE  = 1,
    IS_ON_BOUNDARY= 2,
    IS_DOUBLED = 3
};

inline triangulated_polyhedron
make_triangulated_polyhedron( const polyhedron_data_type&  data)
{
    return triangulated_polyhedron( data.first, 
                                    {
                                        K_MEMORY::vector_view<const double>(data.second[0].begin(), data.second[0].end()),
                                        K_MEMORY::vector_view<const double>(data.second[1].begin(), data.second[1].end()),
                                        K_MEMORY::vector_view<const double>(data.second[2].begin(), data.second[2].end())
                                    } );
}

int main()
{
    std::cout << "Testing triangulated_polyhedron..." << std::endl;
    // ##################### LES TESTS #########################
    // Tetraèdre régulier
    // =========================================================
    polyhedron_data_type tetrahedron_data(
     // Connectivité
     {    {0,1,3},
          {0,2,1},
          {0,3,2},
          {1,2,3}
     },
     // Coordonnées des sommets
     { std::vector<double>{0., 0., 1., 1.}, 
                          {0., 1., 0., 1.}, 
                          {0., 1., 1., 0.} 
     }
    );
    triangulated_polyhedron tetrahedron = make_triangulated_polyhedron(tetrahedron_data);
    auto tetrahedron_classify = [](const K_POST::point3d& pt) {
        if ( (pt.x+pt.y+pt.z < 2) && (pt.x+pt.y > pt.z) && (pt.x+pt.z > pt.y) && (pt.y + pt.z > pt.x) )
            return IS_INSIDE;
        if ( (pt.x+pt.y+pt.z <= 2) && (pt.x+pt.y >= pt.z) && (pt.x+pt.z >= pt.y) && (pt.y + pt.z >= pt.x) )
            return IS_ON_BOUNDARY;
        return IS_OUTSIDE;
    };
    std::vector<K_POST::point3d> points;
    points.reserve(49);
    for ( int ix = -1; ix < 6; ++ix )
        for  (int iy = -1; iy < 6; ++iy )
            for ( int iz = -1; iz < 6; ++iz )
            {
                points.emplace_back(0.25*ix,0.25*iy,0.25*iz);                
            }
    for (const auto& point : points )
    {
        auto classe = tetrahedron_classify(point);
        try{
            int  winding_number = tetrahedron.winding_number(point);
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Tetrahedron failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Tetrahedron failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Tetrahedron failed : point detected on boundary but would not be there !" << std::endl;
        }
    }

    // Octaèdre régulier avec les sommets sur les axes
    // =========================================================
    polyhedron_data_type octahedron_data(
        // Connectivité
        {    {0,2,1},
             {0,4,2},
             {0,1,5},
             {0,5,4},
             {3,1,2},
             {3,5,1},
             {3,2,4},
             {3,4,5} 
        },
        // Coordonnées des sommets
        { std::vector<double>{-1.,  0.,  0.,  1.,  0.,  0.},
                             { 0., -1.,  0.,  0.,  1.,  0.},
                             { 0.,  0., -1.,  0.,  0.,  1.}
        }
    );
    triangulated_polyhedron octahedron = make_triangulated_polyhedron(octahedron_data);
    auto octahedron_classify = [](const K_POST::point3d& pt) {
            double s = std::abs(pt.x) + std::abs(pt.y) + std::abs(pt.z);
            if (s < 1.)
                return IS_INSIDE;
            else if (s == 1.)
                return IS_ON_BOUNDARY;
            return IS_OUTSIDE;
    };

    // Quarter-integer boundaries from -1.25 to 1.25 inclusive.
    points.clear();
    points.reserve(1331);
    for ( int iv = -5; iv < 6; ++iv)
        for ( int jv = -5; jv < 6; ++jv)
            for ( int kv = -5; kv < 6; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = octahedron_classify(point);
        try{
            int winding_number = octahedron.winding_number(point);
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Octahedron failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Octahedron failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Octahedron failed : point detected on boundary but would not be there !" << std::endl;
        }
    }

    // Cube with vertices at (±1, ±1, ±1).
    // =========================================================
    polyhedron_data_type cube_data{
        // Connectivité
        {
            {1,3,2},
            {1,0,4},
            {1,5,7},
            {2,0,1},
            {2,6,4},
            {2,3,7},
            {4,5,1},
            {4,0,2},
            {4,6,7},
            {7,3,1},
            {7,6,2},
            {7,5,4}
        },
        // Coordonnées des sommets
        {
            std::vector<double>{-1., -1., -1., -1., +1., +1., +1., +1.},
                               {-1., -1., +1., +1., -1., -1., +1., +1.},
                               {-1., +1., -1., +1., -1., +1., -1., +1.}
        }
    };
    auto cube = make_triangulated_polyhedron(cube_data);
    assert(std::abs(cube.volume()-8.)<1.E-14);
    auto cube_classify = [](const K_POST::point3d& pt) {
        if ( (-1 < pt.x) && (pt.x < 1) && (-1 < pt.y) && (pt.y < 1) && (-1 < pt.z) && (pt.z < 1) )
            return IS_INSIDE;
        if ((-1 <= pt.x) && (pt.x <= 1) && (-1 <= pt.y) && (pt.y <= 1) && (-1 <= pt.z) && (pt.z <= 1))
            return IS_ON_BOUNDARY;
        return IS_OUTSIDE;
    };

    points.clear();
    points.reserve(1331);
    for ( int iv = -5; iv < 6; ++iv)
        for ( int jv = -5; jv < 6; ++jv)
            for ( int kv = -5; kv < 6; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = cube_classify(point);
        try{
            int winding_number = cube.winding_number(point);
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Cube failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Cube failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Cube failed : point detected on boundary but would not be there !" << std::endl;
        }
    }


    // Paire de cube qui partagent un même sommet à l'origine.
    // =========================================================
    polyhedron_data_type pair_of_cubes_data{
        // Connectivité
        { { 1, 3, 2}, { 1, 0, 4}, { 1, 5, 7},
          { 2, 0, 1}, { 2, 6, 4}, { 2, 3, 7},
          { 4, 5, 1}, { 4, 0, 2}, { 4, 6, 7},
          { 7, 3, 1}, { 7, 6, 2}, { 7, 5, 4},
          { 8,10, 9}, { 8, 7,11}, { 8,12,14},
          { 9, 7, 8}, { 9,13,11}, { 9,10,14},
          {11,12, 8}, {11, 7, 9}, {11,13,14},
          {14,10, 8}, {14,13, 9}, {14,12,11}
        }
        ,
        // Coordonnées des sommets
        { std::vector<double>{-1.,-1.,-1.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1.},
                             {-1.,-1., 0., 0.,-1.,-1., 0., 0., 0., 1., 1., 0., 0., 1., 1.},
                             {-1., 0.,-1., 0.,-1., 0.,-1., 0., 1., 0., 1., 0., 1., 0., 1.}
        }
    };
    auto pair_of_cubes = make_triangulated_polyhedron(pair_of_cubes_data);
    assert(std::abs(pair_of_cubes.volume()-2.0) < 1.E-14);
    auto pair_of_cubes_classify = [] (const K_POST::point3d& pt) {
        if ((-1<pt.x) && (pt.x<0) && (-1<pt.y) && (pt.y<0) && (-1<pt.z) && (pt.z<0))
            return IS_INSIDE;
        if ((0<pt.x) && (pt.x<1) && (0<pt.y) && (pt.y<1) && (0<pt.z) && (pt.z<1))
            return IS_INSIDE;
        if ((-1<=pt.x) && (pt.x<=0) && (-1<=pt.y) && (pt.y<=0) && (-1<=pt.z) && (pt.z<=0))
            return IS_ON_BOUNDARY;
        if ((0<=pt.x) && (pt.x<=1) && (0<=pt.y) && (pt.y<=1) && (0<=pt.z) && (pt.z<=1))
            return IS_ON_BOUNDARY;
        return IS_OUTSIDE;
    };

    points.clear();
    points.reserve(1331);
    for ( int iv = -5; iv < 6; ++iv)
        for ( int jv = -5; jv < 6; ++jv)
            for ( int kv = -5; kv < 6; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = pair_of_cubes_classify(point);
        try{
            int winding_number = pair_of_cubes.winding_number(point);
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Cube failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Cube failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Cube failed : point detected on boundary but would not be there !" << std::endl;
        }
    }

    // Deux pavés droits empilés, un directement au dessus de l'autre
    // ==============================================================
    polyhedron_data_type aligned_stacked_cuboids_data{
        // Connectivité
        {   // 1er cuboid :
            {1,3,2}, {1,0,4}, {1,5,7}, {2,0,1},
            {2,6,4}, {2,3,7}, {4,5,1}, {4,0,2},
            {4,6,7}, {7,3,1}, {7,6,2}, {7,5,4},
            // 2nd cuboid :
            { 9,11,10}, { 9, 8,12}, { 9,13,15}, {10, 8, 9},
            {10,14,12}, {10,11,15}, {12,13, 9}, {12, 8,10},
            {12,14,15}, {15,11, 9}, {15,14,10}, {15,13,12}
        }
        ,
        // Coordonnées des sommets
        { std::vector<double>{ 0., 0., 0., 0., 3., 3., 3., 3., 0., 0., 0., 0., 3., 3., 3., 3.},
                             { 0., 0., 3., 3., 0., 0., 3., 3., 0., 0., 3., 3., 0., 0., 3., 3.},
                             { 0., 1., 0., 1., 0., 1., 0., 1., 2., 3., 2., 3., 2., 3., 2., 3.} }
    };
    auto aligned_stacked_cuboids = make_triangulated_polyhedron(aligned_stacked_cuboids_data);
    assert(std::abs(aligned_stacked_cuboids.volume()-18.0) < 1.E-14);
 
    auto aligned_stacked_cuboids_classify = [](const K_POST::point3d& point) {
        if ((0. < point.x) && (point.x < 3) && (0 < point.y) && (point.y < 3) && (0 < point.z) && (point.z < 1))
            return IS_INSIDE;
        if ((0. < point.x) && (point.x < 3) && (0 < point.y) && (point.y < 3) && (2 < point.z) && (point.z < 3))
            return IS_INSIDE;
        if ((0. <= point.x) && (point.x <= 3) && (0 <= point.y) && (point.y <= 3) && (0 <= point.z) && (point.z <= 1))
            return IS_ON_BOUNDARY;
        if ((0. <= point.x) && (point.x <= 3) && (0 <= point.y) && (point.y <= 3) && (2 <= point.z) && (point.z <= 3))
            return IS_ON_BOUNDARY;
        return IS_OUTSIDE;       
    };
    points.clear();
    points.reserve(3375);
    for ( int iv = -1; iv < 14; ++iv )
        for ( int jv = -1; jv < 14; ++jv )
            for (int kv = -1; kv < 14; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = aligned_stacked_cuboids_classify(point);
        try{
            int winding_number = aligned_stacked_cuboids.winding_number(point);
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Aligned stacked cuboids failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Aligned stacked cuboids failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Aligned stacked cuboids failed : point detected on boundary but would not be there !" << std::endl;
        }
    }

    // Pareil que le test précédent mais avec les pavés non alignés
    // ==============================================================
    polyhedron_data_type misaligned_stacked_cuboids_data{
        // Connectivité
        {   // 1er cuboid :
            {1,3,2}, {1,0,4}, {1,5,7}, {2,0,1},
            {2,6,4}, {2,3,7}, {4,5,1}, {4,0,2},
            {4,6,7}, {7,3,1}, {7,6,2}, {7,5,4},
            // 2nd cuboid :
            { 9,11,10}, { 9, 8,12}, { 9,13,15}, {10, 8, 9},
            {10,14,12}, {10,11,15}, {12,13, 9}, {12, 8,10},
            {12,14,15}, {15,11, 9}, {15,14,10}, {15,13,12}
        }
        ,
        // Coordonnées des sommets
        { std::vector<double>{ 0., 0., 0., 0., 2., 2., 2., 2., 1., 1., 1., 1., 3., 3., 3., 3.},
                             { 0., 0., 2., 2., 0., 0., 2., 2., 1., 1., 3., 3., 1., 1., 3., 3.},
                             { 0., 1., 0., 1., 0., 1., 0., 1., 2., 3., 2., 3., 2., 3., 2., 3.} }
    };
    auto misaligned_stacked_cuboids = make_triangulated_polyhedron(misaligned_stacked_cuboids_data);
    assert(std::abs(misaligned_stacked_cuboids.volume()-8.0) < 1.E-14);
 
    auto misaligned_stacked_cuboids_classify = [](const K_POST::point3d& point) {
        if ((0. < point.x) && (point.x < 2) && (0 < point.y) && (point.y < 2) && (0 < point.z) && (point.z < 1))
            return IS_INSIDE;
        if ((1. < point.x) && (point.x < 3) && (1 < point.y) && (point.y < 3) && (2 < point.z) && (point.z < 3))
            return IS_INSIDE;
        if ((0. <= point.x) && (point.x <= 2) && (0 <= point.y) && (point.y <= 2) && (0 <= point.z) && (point.z <= 1))
            return IS_ON_BOUNDARY;
        if ((1. <= point.x) && (point.x <= 3) && (1 <= point.y) && (point.y <= 3) && (2 <= point.z) && (point.z <= 3))
            return IS_ON_BOUNDARY;
        return IS_OUTSIDE;       
    };
    points.clear();
    points.reserve(3375);
    for ( int iv = -1; iv < 14; ++iv )
        for ( int jv = -1; jv < 14; ++jv )
            for (int kv = -1; kv < 14; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = misaligned_stacked_cuboids_classify(point);
        try{
            int winding_number = misaligned_stacked_cuboids.winding_number(point);
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Misaligned stacked cuboids failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Misaligned stacked cuboids failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Misaligned stacked cuboids failed : point detected on boundary but would not be there !" << std::endl;
        }
    }
    // Cube creux : surface qui consiste de deux cubes, un orienté vers l'extérieur et un orienté vers l'intérieur
    // ===========================================================================================================
    polyhedron_data_type hollow_cube_data{
        // Connectivité
        {   // Premier cube :
            {1,3,2}, {1,0,4}, {1,5,7}, {2,0,1},
            {2,6,4}, {2,3,7}, {4,5,1}, {4,0,2},
            {4,6,7}, {7,3,1}, {7,6,2}, {7,5,4},
            // Second cube avec orientation inversée :
            {10,11, 9}, {12, 8, 9}, {15,13, 9}, { 9, 8,10},
            {12,14,10}, {15,11,10}, { 9,13,12}, {10, 8,12},
            {15,14,12}, { 9,11,15}, {10,14,15}, {12,13,15}
        }
        ,
        // Coordonnées des sommets
        { std::vector<double>{ 0., 0., 0., 0., 3., 3., 3., 3., 1., 1., 1., 1., 2., 2., 2., 2.},
                             { 0., 0., 3., 3., 0., 0., 3., 3., 1., 1., 2., 2., 1., 1., 2., 2.},
                             { 0., 3., 0., 3., 0., 3., 0., 3., 1., 2., 1., 2., 1., 2., 1., 2.} }
    };
    auto hollow_cube = make_triangulated_polyhedron(hollow_cube_data);
    assert(std::abs(hollow_cube.volume()-26.0) < 1.E-14);

    auto hollow_cube_classify = [](const K_POST::point3d& point)
    {
        if ((1 < point.x) && (point.x < 2) && (1 < point.y) && (point.y < 2) && (1 < point.z) && (point.z < 2))
            return IS_OUTSIDE;
        if ((1 <= point.x) && (point.x <= 2) && (1 <= point.y) && (point.y <= 2) && (1 <= point.z) && (point.z <= 2))
            return IS_ON_BOUNDARY;
        if ((0 < point.x) && (point.x < 3) && (0 < point.y) && (point.y < 3) && (0 < point.z) && (point.z < 3))
            return IS_INSIDE;
        if ((0 <= point.x) && (point.x <= 3) && (0 <= point.y) && (point.y <= 3) && (0 <= point.z) && (point.z <= 3))
            return IS_ON_BOUNDARY;
        return IS_OUTSIDE;
    };
    points.clear();
    points.reserve(3375);
    for ( int iv = -1; iv < 14; ++iv )
        for ( int jv = -1; jv < 14; ++jv )
            for (int kv = -1; kv < 14; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = hollow_cube_classify(point);
        try{
            int winding_number = hollow_cube.winding_number(point);
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Hollow cube failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Hollow cube failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Hollow cube failed : point detected on boundary but would not be there !" << std::endl;
        }
    }

    // Cubes emboîtés : la surface consiste de deux cubes, les deux dirigés vers l'extérieur
    //    Les points dans le cube interne devront avoir un indice homologique valant deux.
    // =====================================================================================
    polyhedron_data_type nested_cube_data{
        // Connectivité
        {   // 1er cube :
            {1,3,2}, {1,0,4}, {1,5,7}, {2,0,1},
            {2,6,4}, {2,3,7}, {4,5,1}, {4,0,2},
            {4,6,7}, {7,3,1}, {7,6,2}, {7,5,4},
            // 2ème cube :
            { 9,11,10}, { 9, 8,12}, { 9,13,15}, {10, 8, 9},
            {10,14,12}, {10,11,15}, {12,13, 9}, {12, 8,10},
            {12,14,15}, {15,11, 9}, {15,14,10}, {15,13,12}
        }
        ,
        // Coordonnées des sommets
        { std::vector<double>{ 0., 0., 0., 0., 3., 3., 3., 3., 1., 1., 1., 1., 2., 2., 2., 2.},
                             { 0., 0., 3., 3., 0., 0., 3., 3., 1., 1., 2., 2., 1., 1., 2., 2.},
                             { 0., 3., 0., 3., 0., 3., 0., 3., 1., 2., 1., 2., 1., 2., 1., 2.} }
    };
    auto nested_cube = make_triangulated_polyhedron(nested_cube_data);
    assert(std::abs(nested_cube.volume()-28.0) < 1.E-14);

    auto nested_cube_classify = [] (const K_POST::point3d& point)
    {
        if ((1 < point.x) && (point.x < 2) && (1 < point.y) && (point.y < 2) && (1 < point.z) && (point.z < 2))
            return IS_DOUBLED;
        if ((1 <= point.x) && (point.x <= 2) && (1 <= point.y) && (point.y <= 2) && (1 <= point.z) && (point.z <= 2))
            return IS_ON_BOUNDARY;
        if ((0 < point.x) && (point.x < 3) && (0 < point.y) && (point.y < 3) && (0 < point.z) && (point.z < 3))
            return IS_INSIDE;
        if ((0 <= point.x) && (point.x <= 3) && (0 <= point.y) && (point.y <= 3) && (0 <= point.z) && (point.z <= 3))
            return IS_ON_BOUNDARY;
        return IS_OUTSIDE;        
    };
    points.clear();
    points.reserve(3375);
    for ( int iv = -1; iv < 14; ++iv )
        for ( int jv = -1; jv < 14; ++jv )
            for (int kv = -1; kv < 14; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = nested_cube_classify(point);
        try{
            int winding_number = nested_cube.winding_number(point);
            if ((classe == IS_DOUBLED) && (winding_number != 2))
                std::cerr << "Nested cube failed : point must be twice time inside the polyhedron but homologie indice is "
                          << winding_number << " instead two !" << std::endl;
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Nested cube failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Nested cube failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Nested cube failed : point detected on boundary but would not be there !" << std::endl;
        }
    }
    // Tore
    // ===============================================================
    polyhedron_data_type torus_data{
        // Connectivité
        {   // Facette haute 
            { 8,14,12}, {14, 8,10}, {10,15,14}, {15,10,11},
            {11,13,15}, {13,11, 9}, { 9,12,13}, {12, 9, 8},
            // Facette basse
            { 4, 6, 0}, { 2, 0, 6}, { 6, 7, 2}, { 3, 2, 7},
            { 7, 5, 3}, { 1, 3, 5}, { 5, 4, 1}, { 0, 1, 4},
            // Facettes extérieures
            { 0, 2,10}, {10, 8, 0}, { 2, 3,11}, {11,10, 2},
            { 3, 1, 9}, { 9,11, 3}, { 1, 0, 8}, { 8, 9, 1},
            // Facettes internes
            { 4,12,14}, {14, 6, 4}, { 6,14,15}, {15, 7, 6},
            { 7,15,13}, {13, 5, 7}, { 5,13,12}, {12, 4, 5}
        }
        ,
        // Coordonnées sommets ｢ Dehors/bas  ｣ ｢ Dedans/bas  ｣｢ Dehors/haut  ｣｢ Dedans/haut  ｣(Toutes des facettes carrées)
        { std::vector<double>{ 0., 0., 3., 3., 1., 1., 2., 2., 0., 0., 3., 3., 1., 1., 2., 2.},
                             { 0., 3., 0., 3., 1., 2., 1., 2., 0., 3., 0., 3., 1., 2., 1., 2.},
                             { 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1.} }
    };
    auto torus = make_triangulated_polyhedron(torus_data);
    assert(std::abs(torus.volume()-8.0) < 1.E-14);

    auto torus_classify = [] (const K_POST::point3d& point)
    {
        if ((0 < point.z) && (point.z < 1) && ( ((0 < point.x) && (point.x < 1)) || ((2 < point.x) && (point.x < 3)) ) &&
            (0 < point.y) && (point.y < 3))
            return IS_INSIDE;
        if ((0 < point.z) && (point.z < 1) && ( ((0 < point.y) && (point.y < 1)) || ((2 < point.y) && (point.y < 3)) ) &&
            (0 < point.x) && (point.x < 3))
            return IS_INSIDE;
        if ((0 <= point.z) && (point.z <= 1) && ( ((0 <= point.x) && (point.x <= 1)) || ((2 <= point.x) && (point.x <= 3)) ) &&
            (0 <= point.y) && (point.y <= 3))
            return IS_ON_BOUNDARY;
        if ((0 <= point.z) && (point.z <= 1) && ( ((0 <= point.y) && (point.y <= 1)) || ((2 <= point.y) && (point.y <= 3)) ) &&
            (0 <= point.x) && (point.x <= 3))
            return IS_ON_BOUNDARY;
        return IS_OUTSIDE;
    };
    points.clear();
    points.reserve(1575);
    for ( int iv = -1; iv < 14; ++iv )
        for ( int jv = -1; jv < 14; ++jv )
            for (int kv = -1; kv < 6; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = torus_classify(point);
        try{
            int winding_number = torus.winding_number(point);
            if ((classe == IS_INSIDE) && (winding_number != 1))
                std::cerr << "Torus failed : point must be inside the polyhedron but homologie indice is "
                          << winding_number << " instead one !" << std::endl;
            else if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Torus failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Torus failed : point detected on boundary but would not be there !" << std::endl;
        }
    }
    // Une surface nulle
    // =========================================================================
    polyhedron_data_type empty_data = { {}, 
                                        { std::vector<double>{},
                                                             {},
                                                             {} }
                                      };
    auto empty = make_triangulated_polyhedron(empty_data);
    assert(std::abs(empty.volume()) < 1.E-14);
    points.clear();
    points.reserve(3375);
    for ( int iv = -1; iv < 14; ++iv )
        for ( int jv = -1; jv < 14; ++jv )
            for (int kv = -1; kv < 14; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        int winding_number = empty.winding_number(point);
        if (winding_number != 0)
        {
            std::cerr << "Empty test fail : found a non null winding number : " << winding_number << std::endl;
        }
    }
    // Triangle avec une double face (volume nul, donc). C'est la surface non triviale la plus simple
    // ===============================================================================================
    polyhedron_data_type triangle_data = { {
                                            {0, 1, 2},
                                            {2, 1, 0},
                                           }, 
                                           { std::vector<double>{1., 0., 0.},
                                                                {0., 1., 0.},
                                                                {0., 0., 1.} }
                                         };
    auto triangle = make_triangulated_polyhedron(triangle_data);
    assert(std::abs(triangle.volume()) < 1.E-14);
    auto triangle_classify = [] (const K_POST::point3d& point)
    {
        if ( (0 > point.x) || (point.x > 1) || (0 > point.y) || (point.y > 1) ||
             (0 > point.z) || (point.z > 1) || (point.x + point.y + point.z != 1.) )
            return IS_OUTSIDE;
        return IS_ON_BOUNDARY;
    };

    points.clear();
    points.reserve(7*7*7);
    for ( int iv = -1; iv < 6; ++iv )
        for ( int jv = -1; jv < 6; ++jv )
            for (int kv = -1; kv < 6; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = triangle_classify(point);
        try{
            int winding_number = triangle.winding_number(point);
            if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Triangle failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Triangle failed : point detected on boundary but would not be there !" << std::endl;
        }
    }
    // Un hexadecaèdre qui est topologiquement équivalent à deux pyramides octogonales
    // reliées par la base, mais avec une application dans ℝ³ qui enveloppe la surface
    // deux fois autour de l'origine si bien que son image dans ℝ³ ressemble à la surface
    // d'un octaèdre. Tous les points dans la surface ont un indice de 2.
    polyhedron_data_type hexadecaedre_data = { 
        // Connectivité
        {
            {0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5},
            {0, 5, 6}, {0, 6, 7}, {0, 7, 8}, {0, 8, 1},
            {9, 2, 1}, {9, 3, 2}, {9, 4, 3}, {9, 5, 4},
            {9, 6, 5}, {9, 7, 6}, {9, 8, 7}, {9, 1, 8}        
        },
        // Coordonnées des sommets 
        { std::vector<double>{0., 1., 0.,-1., 0., 1., 0.,-1., 0., 0.},
                             {0., 0., 1., 0.,-1., 0., 1., 0.,-1., 0.},
                             {1., 0., 0., 0., 0., 0., 0., 0., 0.,-1.} }
    };
    auto hexadecaedre = make_triangulated_polyhedron(hexadecaedre_data);
    assert(std::abs(hexadecaedre.volume() - 8./3.) < 1.E-14);
    auto hexadecaedre_classify = [] (const K_POST::point3d& point)
    {
        double s = std::abs(point.x) + std::abs(point.y) + std::abs(point.z);
        if (s==1) return IS_ON_BOUNDARY;
        if (s >1) return IS_OUTSIDE;
        return IS_DOUBLED;
    };
    points.clear();
    points.reserve(11*11*11);
    for ( int iv = -5; iv < 6; ++iv )
        for ( int jv = -5; jv < 6; ++jv )
            for (int kv = -5; kv < 6; ++kv )
                points.emplace_back(0.25*iv, 0.25*jv, 0.25*kv);
    for ( const auto& point : points )
    {
        auto classe = hexadecaedre_classify(point);
        try{
            int winding_number = hexadecaedre.winding_number(point);
            if ((classe == IS_OUTSIDE) && (winding_number != 0))
                std::cerr << "Hexadecaedre failed : point must be outside the polyhedron but homologie indice is "
                          << winding_number << " instead zero !" << std::endl;
            if ((classe == IS_DOUBLED) && (winding_number != 2))
                std::cerr << "Hexadecaedre failed : point must be doubled inside the polyhedron but homologie indice is "
                          << winding_number << " instead two !" << std::endl;
        } catch(std::invalid_argument& err)
        {
            if (classe != IS_ON_BOUNDARY)
                std::cerr << "Hexadecaedre failed : point detected on boundary but would not be there !" << std::endl;
        }
    }
}
#endif
