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
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include "linear_algebra.hpp" 

double K_POST::det(const matrix_2x2_type& A)
{
    return A[0][0]*A[1][1] - A[0][1]*A[1][0];
}
// ----------------------------------------------------------------------------------------------------
auto K_POST::operator * ( const matrix_2x2_type& A, const vector2d& u) -> vector2d
{
    return {
        A[0][0]*u[0] + A[0][1]*u[1],
        A[1][0]*u[0] + A[1][1]*u[1]
           };
}
// ----------------------------------------------------------------------------------------------------
auto K_POST::inverse(const matrix_2x2_type& A) -> matrix_2x2_type
{
    // On utilise la formule des cofacteurs pour une matrice 2x2
    double detA = det(A);
    if (std::abs(detA)<1.E-14) throw std::invalid_argument("The matrix A is not inversible !");
    double inv_det = 1./detA;
    return {
        std::array<double,2>{ inv_det*A[1][1], -inv_det*A[0][1]},
                            {-inv_det*A[1][0],  inv_det*A[0][0]}
           };
}
// ----------------------------------------------------------------------------------------------------
auto
K_POST::inverse_linear_system(const matrix_2x2_type& A, const vector2d& u) -> vector2d
{
    // On utilise la formule des cofacteurs pour une matrice 2x2
    double detA = det(A);
    if (std::abs(detA)<1.E-14) throw std::invalid_argument("The matrix A is not inversible !");
    matrix_2x2_type coMat = {
        std::array<double,2>{ A[1][1], -A[0][1]},
                            {-A[1][0],  A[0][0]}
                            };
    double inv_det = 1./detA;
    return {
            inv_det*(coMat[0][0]*u[0]+coMat[0][1]*u[1]),
            inv_det*(coMat[1][0]*u[0]+coMat[1][1]*u[1])
           };
}
// ----------------------------------------------------------------------------------------------------
std::ostream& 
K_POST::operator << ( std::ostream& out, const vector2d& u)
{
    out << "( " << u[0] << ", " << u[1] << " )";
    return out;
}
// ----------------------------------------------------------------------------------------------------
std::string K_POST::to_string(const matrix_2x2_type& A)
{
    std::ostringstream sout;
    constexpr const int w = 9;
    constexpr const int p = 4;
    sout << "(" << std::fixed 
         << std::setprecision(p) << std::setw(w) << A[0][0] 
         << std::setprecision(p) << std::setw(w) << A[0][1] 
         << " )" << std::endl;
    sout << "(" << std::fixed
         << std::setprecision(p) << std::setw(w) << A[1][0]
         << std::setprecision(p) << std::setw(w) << A[1][1] << " )" << std::endl;
    return sout.str();
}
// ====================================================================================================
double K_POST::det(const matrix_3x3_type& A)
{
    return A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] -
           A[0][2]*A[1][1]*A[2][0] - A[1][2]*A[2][1]*A[0][0] - A[2][2]*A[0][1]*A[1][0];
}
// ----------------------------------------------------------------------------------------------------
auto K_POST::operator * ( const matrix_3x3_type& A, const vector3d& u) -> vector3d
{
    return {
        A[0][0]*u[0] + A[0][1]*u[1] + A[0][2]*u[2],
        A[1][0]*u[0] + A[1][1]*u[1] + A[1][2]*u[2],
        A[2][0]*u[0] + A[2][1]*u[1] + A[2][2]*u[2]
           };
}
// ----------------------------------------------------------------------------------------------------
auto K_POST::inverse(const matrix_3x3_type& A) -> matrix_3x3_type
{
    // On utilise la formule des cofacteurs pour une matrice 3x3
    double detA = det(A);
    if (std::abs(detA)<1.E-14) throw std::invalid_argument("The matrix A is not inversible !");
    double inv_det = 1./detA;
    return {
        std::array<double,3>{inv_det*(A[1][1]*A[2][2]-A[1][2]*A[2][1]), inv_det*(A[2][1]*A[0][2]-A[0][1]*A[2][2]), inv_det*(A[0][1]*A[1][2]-A[1][1]*A[0][2])},
                            {inv_det*(A[2][0]*A[1][2]-A[1][0]*A[2][2]), inv_det*(A[0][0]*A[2][2]-A[2][0]*A[0][2]), inv_det*(A[1][0]*A[0][2]-A[0][0]*A[1][2])},
                            {inv_det*(A[1][0]*A[2][1]-A[2][0]*A[1][1]), inv_det*(A[2][0]*A[0][1]-A[0][0]*A[2][1]), inv_det*(A[0][0]*A[1][1]-A[1][0]*A[0][1])}
           };
}
// ----------------------------------------------------------------------------------------------------
auto
K_POST::inverse_linear_system(const matrix_3x3_type& A, const vector3d& u) -> vector3d
{
    // On utilise la formule des cofacteurs pour une matrice 3x3
    double detA = det(A);
    if (std::abs(detA)<1.E-14) throw std::invalid_argument("The matrix A is not inversible !");
    matrix_3x3_type coMat = {
        std::array<double,3>{A[1][1]*A[2][2]-A[1][2]*A[2][1], A[2][1]*A[0][2]-A[0][1]*A[2][2], A[0][1]*A[1][2]-A[1][1]*A[0][2]},
                            {A[2][0]*A[1][2]-A[1][0]*A[2][2], A[0][0]*A[2][2]-A[2][0]*A[0][2], A[1][0]*A[0][2]-A[0][0]*A[1][2]},
                            {A[1][0]*A[2][1]-A[2][0]*A[1][1], A[2][0]*A[0][1]-A[0][0]*A[2][1], A[0][0]*A[1][1]-A[1][0]*A[0][1]}
                            };
    double inv_det = 1./detA;
    return {
            inv_det*(coMat[0][0]*u[0]+coMat[0][1]*u[1]+coMat[0][2]*u[2]),
            inv_det*(coMat[1][0]*u[0]+coMat[1][1]*u[1]+coMat[1][2]*u[2]),
            inv_det*(coMat[2][0]*u[0]+coMat[2][1]*u[1]+coMat[2][2]*u[2])
           };
}
// ----------------------------------------------------------------------------------------------------
std::string K_POST::to_string(const matrix_3x3_type& A)
{
    std::ostringstream sout;
    constexpr const int w = 9;
    constexpr const int p = 4;
    sout << "(" << std::fixed 
         << std::setprecision(p) << std::setw(w) << A[0][0] 
         << std::setprecision(p) << std::setw(w) << A[0][1] 
         << std::setprecision(p) << std::setw(w) << A[0][2] 
         << " )" << std::endl;
    sout << "|" << std::fixed
         << std::setprecision(p) << std::setw(w) << A[1][0] 
         << std::setprecision(p) << std::setw(w) << A[1][1] 
         << std::setprecision(p) << std::setw(w) << A[1][2] 
         << " |" << std::endl;
    sout << "(" << std::fixed
         << std::setprecision(p) << std::setw(w) << A[2][0]
         << std::setprecision(p) << std::setw(w) << A[2][1] 
         << std::setprecision(p) << std::setw(w) << A[2][2] << " )" << std::endl;
    return sout.str();
}
// ====================================================================================================
auto K_POST::operator * ( const matrix_4x4_type& A, const vector4d& u) -> vector4d
{
    return {
        A[0][0]*u[0]+A[0][1]*u[1]+A[0][2]*u[2]+A[0][3]*u[3],
        A[1][0]*u[0]+A[1][1]*u[1]+A[1][2]*u[2]+A[1][3]*u[3],
        A[2][0]*u[0]+A[2][1]*u[1]+A[2][2]*u[2]+A[2][3]*u[3],
        A[3][0]*u[0]+A[3][1]*u[1]+A[3][2]*u[2]+A[3][3]*u[3]
           };
}
// ----------------------------------------------------------------------------------------------------
auto K_POST::factorize(const matrix_4x4_type& A) -> factorized_matrix_4x4_type
{
    std::array<int,4> pivot{0,1,2,3};// Pivot au depart = identite
    matrix_4x4_type LU(A);

    for (int i = 0; i < 4; i++) 
    {
        // Recherche du plus grand coefficient de la premiere colonne de la sous-matrice
        int imax = i;
        for (int j = i+1; j < 4; ++j )
            if (std::abs(LU[pivot[j]][i]) > std::abs(LU[pivot[imax]][i])) imax = j;
        if (imax != i)
        {
           std::swap(pivot[i], pivot[imax]);
           std::swap(LU[i],LU[imax]);
        }
        double Apiv = LU[i][i];
        if (std::abs(Apiv) < 1.E-14) throw std::invalid_argument("La matrice A est non inversible !");
        double inv_piv = 1./Apiv;
        for (int j = i + 1; j < 4; j++) {
            LU[j][i] *= inv_piv;

            for (int k = i + 1; k < 4; k++)
                LU[j][k] -= LU[j][i] * LU[i][k];
        }
    }
    return {LU,pivot};
}
// ----------------------------------------------------------------------------------------------------
auto K_POST::inverse_linear_system(const factorized_matrix_4x4_type& PLU, const vector4d& u) -> vector4d
{
    const matrix_4x4_type&   LU = PLU.first;
    const std::array<int,4>& pivot = PLU.second;
    vector4d v{u[pivot[0]],u[pivot[1]],u[pivot[2]],u[pivot[3]]};
    // Descente
    for (int i = 0; i < 4; i++) 
    {
        for (int k = 0; k < i; k++)
            v[i] -= LU[i][k] * v[k];
    }
    // Remontee :
    for (int i = 3; i >= 0; i--) 
    {
        for (int k = i + 1; k < 4; k++)
            v[i] -= LU[i][k] * v[k];

        v[i] /= LU[i][i];
    }
    return v;
}
// ----------------------------------------------------------------------------------------------------
std::ostream& 
K_POST::operator << ( std::ostream& out, const vector4d& u)
{
    out << "( " << u[0] << ", " << u[1] << ", " << u[2] << ", " << u[3] << " )";
    return out;
}
// ----------------------------------------------------------------------------------------------------
std::string K_POST::to_string(const matrix_4x4_type& A)
{
    std::ostringstream sout;
    constexpr const int w = 9;
    constexpr const int p = 4;
    sout << "(" << std::fixed 
         << std::setprecision(p) << std::setw(w) << A[0][0] 
         << std::setprecision(p) << std::setw(w) << A[0][1] 
         << std::setprecision(p) << std::setw(w) << A[0][2] 
         << std::setprecision(p) << std::setw(w) << A[0][3] 
         << " )" << std::endl;
    sout << "|" << std::fixed
         << std::setprecision(p) << std::setw(w) << A[1][0] 
         << std::setprecision(p) << std::setw(w) << A[1][1] 
         << std::setprecision(p) << std::setw(w) << A[1][2] 
         << std::setprecision(p) << std::setw(w) << A[1][3] 
         << " |" << std::endl;
    sout << "|" << std::fixed
         << std::setprecision(p) << std::setw(w) << A[2][0] 
         << std::setprecision(p) << std::setw(w) << A[2][1] 
         << std::setprecision(p) << std::setw(w) << A[2][2] 
         << std::setprecision(p) << std::setw(w) << A[2][3] 
         << " |" << std::endl;
    sout << "(" << std::fixed
         << std::setprecision(p) << std::setw(w) << A[3][0]
         << std::setprecision(p) << std::setw(w) << A[3][1] 
         << std::setprecision(p) << std::setw(w) << A[3][2]
         << std::setprecision(p) << std::setw(w) << A[3][3] << " )" << std::endl;
    return sout.str();

}
// ----------------------------------------------------------------------------------------------------
std::string K_POST::to_string(const factorized_matrix_4x4_type& PLU)
{
    const matrix_4x4_type&   LU = PLU.first;
    const std::array<int,4>& pivot = PLU.second;    
    std::ostringstream sout;
    constexpr const int w = 9;
    constexpr const int p = 4;
    sout << "(" << std::fixed 
         << std::setprecision(p) << std::setw(w) << (pivot[0] == 0 ? 1 : 0) 
         << std::setprecision(p) << std::setw(w) << (pivot[0] == 1 ? 1 : 0)
         << std::setprecision(p) << std::setw(w) << (pivot[0] == 2 ? 1 : 0)
         << std::setprecision(p) << std::setw(w) << (pivot[0] == 3 ? 1 : 0)
         << " ) ("
         << std::setprecision(p) << std::setw(w) << 1. 
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << 0. 
         << " ) ("
         << std::setprecision(p) << std::setw(w) << LU[0][0] 
         << std::setprecision(p) << std::setw(w) << LU[0][1] 
         << std::setprecision(p) << std::setw(w) << LU[0][2] 
         << std::setprecision(p) << std::setw(w) << LU[0][3] 
         << " )" << std::endl;
    sout << "|" << std::fixed 
         << std::setprecision(p) << std::setw(w) << (pivot[1] == 0 ? 1 : 0) 
         << std::setprecision(p) << std::setw(w) << (pivot[1] == 1 ? 1 : 0)
         << std::setprecision(p) << std::setw(w) << (pivot[1] == 2 ? 1 : 0)
         << std::setprecision(p) << std::setw(w) << (pivot[1] == 3 ? 1 : 0)
         << " | |"
         << std::setprecision(p) << std::setw(w) << LU[1][0]
         << std::setprecision(p) << std::setw(w) << 1.
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << 0. 
         << " | |"
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << LU[1][1] 
         << std::setprecision(p) << std::setw(w) << LU[1][2] 
         << std::setprecision(p) << std::setw(w) << LU[1][3] 
         << " |" << std::endl;
    sout << "|" << std::fixed 
         << std::setprecision(p) << std::setw(w) << (pivot[2] == 0 ? 1 : 0) 
         << std::setprecision(p) << std::setw(w) << (pivot[2] == 1 ? 1 : 0)
         << std::setprecision(p) << std::setw(w) << (pivot[2] == 2 ? 1 : 0)
         << std::setprecision(p) << std::setw(w) << (pivot[2] == 3 ? 1 : 0)
         << " |.|"
         << std::setprecision(p) << std::setw(w) << LU[2][0]
         << std::setprecision(p) << std::setw(w) << LU[2][1]
         << std::setprecision(p) << std::setw(w) << 1. 
         << std::setprecision(p) << std::setw(w) << 0. 
         << " |.|"
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << LU[2][2] 
         << std::setprecision(p) << std::setw(w) << LU[2][3] 
         << " |" << std::endl;
    sout << "(" << std::fixed 
         << std::setprecision(p) << std::setw(w) << (pivot[3] == 0 ? 1 : 0) 
         << std::setprecision(p) << std::setw(w) << (pivot[3] == 1 ? 1 : 0)
         << std::setprecision(p) << std::setw(w) << (pivot[3] == 2 ? 1 : 0)
         << std::setprecision(p) << std::setw(w) << (pivot[3] == 3 ? 1 : 0)
         << " ) ("
         << std::setprecision(p) << std::setw(w) << LU[3][0]
         << std::setprecision(p) << std::setw(w) << LU[3][1]
         << std::setprecision(p) << std::setw(w) << LU[3][2]
         << std::setprecision(p) << std::setw(w) << 1. 
         << " ) ("
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << 0. 
         << std::setprecision(p) << std::setw(w) << LU[3][3] 
         << " )" << std::endl;

    return sout.str();
}

//#define __MAIN_TEST__
#if defined(__MAIN_TEST__)
using namespace K_POST;
int main()
{
    matrix_2x2_type A2{ std::array<double,2>{ 3,1},
                                            {-1,3}
                      };
    std::cout << "A2 : " << std::endl << to_string(A2) << std::endl;
    vector2d u2{1.,2.};
    std::cout << "u2 = " << u2 << std::endl;
    double det_A2 = det(A2);
    assert(std::abs(det_A2-10.)<1.E-14);
    vector2d v2 = A2*u2;
    std::cout << "v2 = A2.u2 = " << v2 << std::endl;
    assert(std::abs(v2[0]-5.)<1.E-14);
    assert(std::abs(v2[1]-5.)<1.E-14);

    vector2d w2_1 = inverse_linear_system(A2, v2);
    std::cout << "w2_1 = A2^-1.v2 = " << w2_1 << std::endl;
    assert(std::abs(w2_1[0]-u2[0])<1.E-14);
    assert(std::abs(w2_1[0]-u2[0])<1.E-14);

    auto invA2 = inverse(A2);
    std::cout << "A2⁻¹ : " << std::endl << to_string(invA2) << std::endl;
    vector2d w2_2 = invA2*v2;
    std::cout << "w2_2 = A2^-1.v2 = " << w2_2 << std::endl;
    assert(std::abs(w2_2[0]-u2[0])<1.E-14);
    assert(std::abs(w2_2[0]-u2[0])<1.E-14);
    // ##############################################################
    matrix_3x3_type A{ std::array<double,3>{ 3, 1, 1},
                                           {-1, 3, 1},
                                           {-1,-1, 3}
                     };
    std::cout << "A : " << std::endl << to_string(A) << std::endl;
    vector3d u{1.,2.,3.};
    std::cout << "u = " << std::string(u) << std::endl;
    double det_A = det(A);
    assert(std::abs(det_A-36.)<1.E-14);
    vector3d v = A*u;
    std::cout << "v = A.u = " << std::string(v) << std::endl;
    assert(std::abs(v.x-8.)<1.E-14);
    assert(std::abs(v.y-8.)<1.E-14);
    assert(std::abs(v.z-6.)<1.E-14);

    vector3d w = inverse_linear_system(A, v);
    std::cout << "w = A^-1.v = " << std::string(w) << std::endl;
    assert(std::abs(w.x-u.x)<1.E-14);
    assert(std::abs(w.y-u.y)<1.E-14);
    assert(std::abs(w.z-u.z)<1.E-14);

    auto invA = inverse(A);
    std::cout << "A^-1 : " << std::endl << to_string(invA) << std::endl;
    vector3d w2 = invA*v;
    std::cout << "w2 = A^-1.v = " << std::string(w2) << std::endl;
    assert(std::abs(w2.x-u.x)<1.E-14);
    assert(std::abs(w2.y-u.y)<1.E-14);
    assert(std::abs(w2.z-u.z)<1.E-14);
    // ##############################################################
    matrix_4x4_type B{ std::array<double,4>{0.5, 1, 1, 1},
                                           {-1, 4, 1, 1},
                                           {-1,-1, 4, 1},
                                           {-1,-1,-1, 2}
                     };
    std::cout << "B : " << std::endl << to_string(B) << std::endl;
    vector4d u4{1.,2.,3.,4.};
    std::cout << "u4 = " << u4 << std::endl;
    vector4d v4 = B*u4;
    std::cout << "v4 = B.u4 = " << v4 << std::endl;

    auto PLU = factorize(B);
    std::cout << "PLU(B) : " << std::endl << to_string(PLU) << std::endl;
    auto w4 = inverse_linear_system(PLU, v4);
    std::cout << "w4 = B^-1.v4 = " << w4 << std::endl;
    assert(std::abs(w4[0]-u4[0]) < 1.E-14);
    assert(std::abs(w4[1]-u4[1]) < 1.E-14);
    assert(std::abs(w4[2]-u4[2]) < 1.E-14);
    assert(std::abs(w4[3]-u4[3]) < 1.E-14);
    return EXIT_SUCCESS;
}
#endif
