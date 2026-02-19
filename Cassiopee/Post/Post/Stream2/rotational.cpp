#include "rotational.hpp"
#include "volume.hpp"

namespace K_POST
{

    vector3d compute_part_of_rot_on_triangle( const point3d & pi, const point3d & pj, const point3d & pk,
                                              const vector3d& vi, const vector3d& vj, const vector3d& vk)
    {
        // # (1/6)(vᵢ+vⱼ+vₖ)∧(Eᵢⱼ∧Eᵢₖ)                                  
        vector3d E_ij (pi,pj);
        vector3d E_ik (pi,pk);
        return (1./6.)*( (vi+vj+vk) ^ (E_ij ^ E_ik) );
    }

    vector3d compute_part_of_rot_on_quadrangle( const point3d & pi, const point3d & pj, const point3d & pk, const point3d & pl,
                                                const vector3d& vi, const vector3d& vj, const vector3d& vk, const vector3d& vl)
    {
        // # (1/12)[vᵢ∧Eₗⱼ∧(Eᵢₖ+Eᵢₗ) + vⱼ∧Eᵢₖ∧(Eⱼₗ+Eⱼᵢ) +                
        // #        vₖ∧Eₗⱼ∧(Eₗₖ+Eᵢₖ)  + vₗ∧Eᵢₖ∧(Eⱼₗ+Eᵢₗ)  ]                
        vector3d E_jl(pj,pl), E_ik(pi,pk), E_il(pi,pl), E_ij(pi,pj), E_kl(pk,pl);
        vector3d r = (1./12.)*( (vi^((E_ik+E_il)^E_jl)) + (vj^(E_ik^(E_jl-E_ij))) +
                          (vk^((E_ik-E_kl)^E_jl)) + (vl^(E_ik^(E_il+E_jl))));;
        return r;
    }

    vector3d compute_rotational_on_tetrahedra( const point3d& p1, const point3d& p2,
                                               const point3d& p3, const point3d& p4,
                                               const vector3d& v1, const vector3d& v2,
                                               const vector3d& v3, const vector3d& v4 )
    {
        E_Float vol = compute_volume_tetrahedra(p1, p2, p3, p4);
        return (-1./vol)*(compute_part_of_rot_on_triangle(p1, p3, p2, v1, v3, v2) +
                          compute_part_of_rot_on_triangle(p1, p2, p4, v1, v2, v4) +
                          compute_part_of_rot_on_triangle(p2, p3, p4, v2, v3, v4) +
                          compute_part_of_rot_on_triangle(p3, p1, p4, v3, v1, v4));
    }

    vector3d compute_rotational_on_pyramid( const point3d& p1, const point3d& p2,
                                            const point3d& p3, const point3d& p4, const point3d& p5,
                                            const vector3d& v1, const vector3d& v2,
                                            const vector3d& v3, const vector3d& v4, const vector3d& v5 )
    {
        // Volume de la pyramide :
        E_Float vol = compute_volume_pyramid(p1, p2, p3, p4, p5);

        vector3d E_12(p1,p2), E_13(p1,p3), E_14(p1,p4), E_15(p1,p5);
        vector3d E_23(p2,p3), E_24(p2,p4), E_34(p3,p4);
    // intégrale rot au travers de la base de la pyramide :        
    // #                 = (1/8){v₁∧[E₂₄∧E₁₃+(1/3)E₃₄∧E₁₂+(1/3)E₄₁∧E₂₃)       +
    // #                         v₂∧[E₂₄∧E₁₃+(1/3)E₄₃∧E₁₂+(1/3)E₄₁∧E₂₃]       +
    // #                         v₃∧[E₂₄∧E₁₃+(1/3)E₄₃∧E₁₂+(1/3)E₁₄∧E₂₃]       +
    // #                         v₄∧[E₂₄∧E₁₃+(1/3)E₃₄∧E₁₂+(1/3)E₁₄∧E₂₃]}       
        vector3d E_2413(3*E_24^E_13), E_1234(E_12^E_34), E_1423(E_14^E_23);
        vector3d rot = (1./24.)*((v1^(E_2413-E_1234-E_1423)) +
                                 (v2^(E_2413+E_1234-E_1423)) +
                                 (v3^(E_2413+E_1234+E_1423)) +
                                 (v4^(E_2413-E_1234+E_1423)) );
    // Intégrale sur les quatres autres faces (triangulaires)
    //  T₁₂₅ -- T₂₃₅ -- T₃₄₅ -- T₄₁₅
        rot += compute_part_of_rot_on_triangle(p1,p2,p5,v1,v2,v5) +
               compute_part_of_rot_on_triangle(p2,p3,p5,v2,v3,v5) +
               compute_part_of_rot_on_triangle(p3,p4,p5,v3,v4,v5) +
               compute_part_of_rot_on_triangle(p4,p1,p5,v4,v1,v5);

        return (-1./vol)*rot;
    }

    vector3d compute_rotational_on_pentahedra( const point3d& p1, const point3d& p2, const point3d& p3, 
                                               const point3d& p4, const point3d& p5, const point3d& p6,
                                               const vector3d& v1, const vector3d& v2, const vector3d& v3,
                                               const vector3d& v4, const vector3d& v5, const vector3d& v6)
    {
        // Calcul volume :
        E_Float vol = compute_volume_pentaedra(p1, p2, p3, p4, p5, p6);
        // Q₁₂₅₄ - Q₂₃₆₅ - Q₃₁₄₆ - T₁₃₂ - T₄₅₆
        return (-1./vol)*(compute_part_of_rot_on_quadrangle(p1,p2,p5,p4,v1,v2,v5,v4) +
                          compute_part_of_rot_on_quadrangle(p2,p3,p6,p5,v2,v3,v6,v5) +
                          compute_part_of_rot_on_quadrangle(p3,p1,p4,p6,v3,v1,v4,v6) +
                          compute_part_of_rot_on_triangle(p1,p3,p2,v1,v3,v2)         +
                          compute_part_of_rot_on_triangle(p4,p5,p6,v4,v5,v6) );
    }

    vector3d compute_rotational_on_hexaedra( const point3d& p1, const point3d& p2, const point3d& p3, const point3d& p4,
                                             const point3d& p5, const point3d& p6, const point3d& p7, const point3d& p8,
                                             const vector3d& v1, const vector3d& v2, const vector3d& v3, const vector3d & v4,
                                             const vector3d& v5, const vector3d& v6, const vector3d& v7, const vector3d& v8 )
    {
        // Calcul du volume
        E_Float vol = compute_volume_hexaedra(p1, p2, p3, p4, p5, p6, p7, p8);
        // Q₁₄₃₂ -- Q₁₂₆₅ -- Q₂₃₇₆ -- Q₃₄₈₇ -- Q₁₅₈₄ -- Q₅₆₇₈
        return (-1./vol)*(compute_part_of_rot_on_quadrangle(p1,p4,p3,p2,v1,v4,v3,v2) +
                          compute_part_of_rot_on_quadrangle(p1,p2,p6,p5,v1,v2,v6,v5) +
                          compute_part_of_rot_on_quadrangle(p2,p3,p7,p6,v2,v3,v7,v6) +
                          compute_part_of_rot_on_quadrangle(p3,p4,p8,p7,v3,v4,v8,v7) +
                          compute_part_of_rot_on_quadrangle(p1,p5,p8,p4,v1,v5,v8,v4) +
                          compute_part_of_rot_on_quadrangle(p5,p6,p7,p8,v5,v6,v7,v8) );
    }

    vector3d compute_rotational_on_ngon( const std::vector<face>& faces,  
                                         const std::array<K_MEMORY::vector_view<const E_Float>,3>& vel )
    {
        vector3d rot{0.,0.,0.};

        for (const auto& f : faces)
        {
            E_Int i = f.indices_vertices[0];
            const point3d&  pi = f.get_vertex(0);
            const vector3d& vi{vel[0][i],vel[1][i],vel[2][i]};
            E_Int j = f.indices_vertices[1];
            const point3d&  pj = f.get_vertex(1);
            const vector3d& vj{vel[0][j],vel[1][j],vel[2][j]};
            E_Int k = f.indices_vertices[2];
            const point3d&  pk = f.get_vertex(2);
            const vector3d& vk{vel[0][k],vel[1][k],vel[2][k]};
            if (f.number_of_vertices() == 3)// C'est un triangle
            {
                rot += compute_part_of_rot_on_triangle(pi, pj, pk, vi, vj, vk);
            }
            else if (f.number_of_vertices() == 4)// C'est un quadrangle
            {
                E_Int l = f.indices_vertices[3];
                const point3d&  pl = f.get_vertex(3);
                const vector3d& vl{vel[0][k],vel[1][k],vel[2][l]};
                rot += compute_part_of_rot_on_quadrangle(pi, pj, pk, pl, vi, vj, vk, vl);
            }
        }
        E_Float vol = compute_volume_ngon(faces);
        return (-1./vol)*rot;
    }

}
