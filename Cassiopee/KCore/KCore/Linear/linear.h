/*    
    Copyright 2013-2026 ONERA.

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

#ifndef _KCORE_LINEAR_H
#define _KCORE_LINEAR_H
#include "Def/DefTypes.h"
#include "Fld/FldArray.h"
#include "Def/DefFunction.h"
#include <vector>
#include "Math/math.h"

// Les matrices n x m (n lignes, m colonnes) sont numerotes : i + j*n 

namespace K_LINEAR
{
  /* Solve A x = B square systems
     IN: A: n x n matrix
     IN: B: n x m vector: m seconds membres
     IN: method: "gauss" ou "bicg"
     OUT: X: n x m : m solutions.
  */
  E_Int solve(E_Int n, E_Int m, 
              E_Float* A, E_Float* B, E_Float* X, char* method);
  
  /* Produit matrice vecteur
     IN: A: n x m matrix
     IN: B: m vector
     OUT: C: n vector
  */
  void prodv(E_Int n, E_Int m, E_Float* A, E_Float* B, E_Float* C);
  
  /* Produit matrice matrice
     IN: A: n x m matrix
     IN: B: m x p matrix
     OUT: C: n x p matrix
  */
  void prodm(E_Int n, E_Int m, E_Int p, E_Float* A, E_Float* B, E_Float* C);

  /* Inverse a 2x2 matrix.
     IN: A: 2x2 matrix
     OUT: B: 2x2 matrix
     Retourne 0 si la matrice n'est pas inversible */
  E_Int inv2(E_Float* A, E_Float* B);
  
  /* Inverse a 3x3 matrix.
     IN: A: 3x3 matrix
     OUT: B: 3x3 matrix
     Retourne 0 si la matrice n'est pas inversible */
  E_Int inv3(E_Float* A, E_Float* B);

  /* Valeurs et vecteurs propre de matrices symetriques */
  /* IN: a00, ...: matrix coefficients (symmetric)
     OUT: lambda0, lambda1: eigen values
     OUT: v0, v1: eigen vectors */
  E_Int eigen2(E_Float a00, E_Float a11, E_Float a10, 
               E_Float& lambda0, E_Float& lambda1, E_Float* v0, E_Float* v1);
  /* IN: a00, ...: matrix coefficients (symmetric)
     OUT: lambda0, lambda1, lambda2: eigen values */
  E_Int root3(E_Float a00, E_Float a01, E_Float a02,
              E_Float a11, E_Float a12, E_Float a22,
              E_Float& lambda0, E_Float& lambda1, E_Float& lambda2);
  /* IN: a00, ...: matrix coefficients (symmetric)
     OUT: lambda0, lambda1, lambda2: eigen values
     OUT: v0, v1, v3: eigen vectors */
  E_Int eigen3(E_Float a00, E_Float a01, E_Float a02,
               E_Float a11, E_Float a12, E_Float a22,
               E_Float& lambda0, E_Float& lambda1, E_Float& lambda2, 
               E_Float* v0, E_Float* v1, E_Float* v2);
  /* Fonctions internes pour eigen3bis */
  void tred2(E_Float V[3][3], E_Float d[3], E_Float e[3]);
  void tql2(E_Float V[3][3], E_Float d[3], E_Float e[3]);
  /* IN: a00, ...: matrix coefficients symmetric
     OUT: lambda0, lambda1, lambda2: eigen values
     OUT: v0, v1, v3: eigen vectors */
  E_Int eigen3bis(E_Float A[3][3], E_Float V[3][3], E_Float d[3]);
  
  /* Internes */
  void getComplement31(E_Float u0, E_Float u1, E_Float u2,
                       E_Float v0, E_Float v1, E_Float v2,
                       E_Float& w0, E_Float& w1, E_Float& w2);
  void getComplement32(E_Float u0, E_Float u1, E_Float u2,
                       E_Float& v0, E_Float& v1, E_Float& v2,
                       E_Float& w0, E_Float& w1, E_Float& w2);
  E_Int rank3(E_Float a00, E_Float a01, E_Float a02,
              E_Float a11, E_Float a12, E_Float a22,
              E_Float& b00, E_Float& b01, E_Float& b02,
              E_Float& b10, E_Float& b11, E_Float& b12,
              E_Float& b20, E_Float& b21, E_Float& b22);


  
  /* Inversion of a symmetric positive definite matrix using Cholesky decomposition
     IN: A  : n x n    matrix
     IN: rows: #rows   vector, rows to compute A^-1
     OUT: A_1 : n*#rows  rows of the inverse 
     OUT: ok 1: inversion completed
             0: A is singular  */
  E_Int inv(E_Int n, E_Float* A, E_Float* A_1, const std::vector<E_Int>& rows);
  
  /* Cholesky decomposition : A = RT R
     IN : A: n x n symmetric positive definite matrix 
     OUT: R: n x n upper triangular matrix
     OUT: ok: 1  A is positive definite
              0  A is not positive definite  */
  E_Int cholesky(E_Int n, E_Float* A, E_Float* R);

  /* Back substitution algorithm : Solve Ax=b
     IN : A: n x n upper triangular matrix
     IN : b: n     vector
     OUT: b: n     solution
     OUT: ok 1: back substitution completed
             0: A is singular  */
  E_Int backSubstitution(E_Int n, E_Float* A, E_Float* b);


  /* Solve linear system A.x = b within tolerance tol
     IN: A: n x n matrix
     IN: b: n vector
     IN: n: size of square matrix A
     OUT: x: n solution vector
  */
  E_Int BiCGStab(const E_Float *, const E_Float *, const E_Int, E_Float *, const E_Float tol = 1e-6);
  
  /* Eigenvalues and eigenvectors of 3x3 symmetric matrix
     IN: A: 3x3 symmetric matrix
     OUT: L: eigenvalues
     OUT: v1, v2, v3: eigenvectors
  */
  void sym3mat_eigen(const E_Float [6], E_Float L[3], E_Float v1[3],
    E_Float v2[3], E_Float v3[3], const E_Float tol = 1e-12);
}
#endif
