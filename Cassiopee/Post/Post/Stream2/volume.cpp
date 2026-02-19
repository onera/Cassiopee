
#include "volume.hpp"

inline E_Float triple_product( const K_POST::vector3d& v1, const K_POST::vector3d& v2, const K_POST::vector3d& v3)
{
    return  (v1|(v2^v3));   
}

namespace K_POST
{

    E_Float compute_volume_tetrahedra( const point3d &p1, const point3d &p2, const point3d &p3, const point3d &p4 )
    {
        constexpr const E_Float onesixth = 1./6.;

        vector3d e14(p1,p4);
        vector3d e23(p2,p3); vector3d e24(p2,p4);

        return onesixth*triple_product(e14, e23, e24);
    }

    E_Float compute_volume_pyramid( const point3d &p1, const point3d &p2, const point3d &p3,
                                  const point3d &p4, const point3d &p5)
    {
        constexpr const E_Float onetwelth = 1./12.;
        // V(P) = (1/12){[E₁₂;E₁₃;E₁₅]+[E₁₂+E₁₃;E₁₄;E₁₅]+[E₂₃;E₂₄;E₂₅]}
        vector3d e12(p1,p2); vector3d e13(p1,p3); 
        vector3d e15(p1,p5); vector3d e14(p1,p4);
        vector3d e23(p2,p3); vector3d e24(p2,p4); vector3d e25(p2,p5);

        return onetwelth*(triple_product(e12, e13, e15)+triple_product(e12+e13,e14,e15)+triple_product(e23,e24,e25));
    } 

    E_Float compute_volume_pentaedra( const point3d &p1, const point3d &p2, const point3d &p3,
                                     const point3d &p4, const point3d &p5, const point3d& p6)
    {
        constexpr const E_Float onesixth = 1./6.;

        // V = (1/6){[E₁₄;E₄₅;E₄₆] + (1/2)[E₁₃; E₃₄; E₁₆] + [E₁₃;E₂₆;E₃₅] + (1/2)[E₂₆;E₆₅;E₃₅] + (1/2)[E₁₂;E₁₅;E₂₄]}
        vector3d e12(p1,p2); vector3d e14(p1,p4); vector3d e13(p1,p3); 
        vector3d e15(p1,p5); vector3d e16(p1,p6);
        vector3d e24(p2,p4); vector3d e26(p2,p6);
        vector3d e34(p3,p4); vector3d e35(p3,p5);
        vector3d e45(p4,p5); vector3d e46(p4,p6);
        vector3d e65(p6,p5);

        return onesixth*(triple_product(e14, e45, e46)+triple_product(e13, e26, e35) +
                    0.5*(triple_product(e13, e34, e16)+triple_product(e26, e65, e35) + triple_product(e12, e15, e24)) );
    }

    E_Float compute_volume_hexaedra ( const point3d &p1, const point3d &p2, const point3d &p3,
                                     const point3d &p4, const point3d &p5, const point3d &p6,
                                     const point3d &p7, const point3d &p8 )
    {
        constexpr const E_Float onesixth = 1./6.;
        // V = (1/6){[E₁₂;E₆₃;E₂₇] + [E₁₄;E₄₇;E₈₃] + [E₁₅;E₈₆;E₅₇]} + (1/12){[E₁₂;E₁₄;E₁₃] + [E₁₅;E₁₂;E₁₆] + [E₂₃;E₃₇;E₂₆] +
        //                                                                   [E₄₃;E₄₈;E₃₇] + [E₁₅;E₅₈;E₁₄] + [E₅₆;E₆₇;E₅₈]}
        vector3d e12(p1,p2); vector3d e13(p1,p3); vector3d e14(p1,p4);
        vector3d e15(p1,p5); vector3d e16(p1,p6);
        vector3d e23(p2,p3); vector3d e26(p2,p6); vector3d e27(p2,p7);
        vector3d e37(p3,p7);
        vector3d e43(p4,p3); vector3d e47(p4,p7); vector3d e48(p4,p8);
        vector3d e56(p5,p6); vector3d e57(p5,p7); vector3d e58(p5,p8);
        vector3d e63(p6,p3); vector3d e67(p6,p7);
        vector3d e83(p8,p3); vector3d e86(p8,p6);        

        return onesixth*(triple_product(e12, e63, e27)+triple_product(e14, e47, e83)+triple_product(e15,e86,e57) +
                         0.5*(triple_product(e12, e14, e13)+triple_product(e15,e12,e16)+triple_product(e23,e37,e26)+
                              triple_product(e43, e48, e37)+triple_product(e15,e58,e14)+triple_product(e56,e67,e58)) );
    }

    E_Float compute_volume_ngon( const std::vector<face>& faces )
    {
        // (1/2)[pᵢ;Eᵢⱼ;Eᵢₖ] pour le triangle
        // (1/2)[pᵢ;Eₗⱼ;Eᵢₖ] + (1/4)[Eᵢⱼ;Eⱼₖ;Eᵢₗ] pour le quadrangle
        E_Float vol = 0.;
        for (const auto& f : faces)
        {
            const point3d& pi = f.get_vertex(0);
            const point3d& pj = f.get_vertex(1);
            const point3d& pk = f.get_vertex(2);
            vector3d e0i(pi.x,pi.y,pi.z);
            vector3d eij(pi,pj);
            vector3d eik(pi,pk);
            if (f.number_of_vertices() == 3)// C'est un triangle
            {
                //const point3d& pk = f.get_vertex(2);
                vol += 0.5*triple_product(e0i,eij,eik);
            }
            else if (f.number_of_vertices() == 4)// C'est un quadrangle
            {
                const point3d& pl = f.get_vertex(3);
                vector3d elj(pl,pj);
                vector3d ejk(pj,pk);
                vector3d eil(pi,pl);
                vol += 0.5*(triple_product(e0i, elj, eik)+0.5*triple_product(eij,ejk,eil));
            }
        }
        return vol/3.;
    }
}
