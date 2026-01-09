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
#include <vector>
#include "linear.h"
#include "iostream"

using namespace K_FUNC;
using namespace std;


//=============================================================================
/* Inversion of a symmetric positive definite matrix using Cholesky 
   decomposition
   IN:  A: n x n matrix
   IN: rows: #rows vector, rows to compute A^-1
   OUT: A_1: n*#rows  rows of the inverse 
   OUT: ok 1: inversion completed
           0: A is singular
*/
//=============================================================================
E_Int K_LINEAR::inv(E_Int n, E_Float* A, E_Float* A_1, const vector<E_Int>& rows)
{

  E_Int ok = 1;

  vector<E_Float> R(n*n);
  vector<E_Float> R_1(n*n);

  for (E_Int ii=0; ii<n*n; ii++)
  {R[ii]=0; R_1[ii]=0;} 

  ok = cholesky(n, A, &R[0]);
  
  if (ok==0) return ok;

  vector<E_Float> ej(n);
  // R^-1
  for (E_Int j=0; j<n; j++)
  {
    // jth column of In
    //E_Float ej[n];
    for (E_Int ii=0; ii<n; ii++) ej[ii]=0; 
    ej[j] = 1;
    
    ok=backSubstitution(n, &R[0], &ej[0]);
    if (ok==0) return ok;
    
    //ej is now the jth column of R^-1 
    for (E_Int i=0; i<n; i++) {R_1[i+j*n]=ej[i];}   
  }
    
  // A^-1 = R^-1 RT^-1 
  // Only row in vector rows of A^-1 are actually computed
  for (unsigned int i=0; i < rows.size(); i++)
    for (E_Int j=0; j<n; j++) 
    {
      A_1[i+j*rows.size()]=0;
      // R_1 is also upper triangular
      for (E_Int k=j; k<n; k++)
      {
        A_1[i+j*rows.size()] += R_1[rows[i]+k*n]*R_1[j+k*n];
      }
    }

  return ok;
}




//=============================================================================
/* Cholesky decomposition : A = RT R
   IN : A: n x n symmetric positive definite matrix 
   OUT: R: n x n upper triangular matrix
   OUT ok : 1 if A is positive definite
            0 if A is not positive definite
*/
//=============================================================================

E_Int K_LINEAR::cholesky(E_Int n, E_Float* A, E_Float* R)
{
  E_Int ok = 1;
  
  if (R == NULL)
    ok=0;
  

  for (E_Int i = 0; i < n; i++)
    for (E_Int j = 0; j < (i+1); j++) 
    {
      E_Float s = 0.;
      for (E_Int k = 0; k < j; k++)
        s += R[i * n + k] * R[j * n + k];
      if ( i==j )
      {
        if (A[i * n + i] - s < 0.) 
        {
          ok=0;
          //printf("Cholesky: Matrix is not positive definite 1");
        }
        R[i * n + j] = sqrt(A[i * n + i] - s);
      }
      else
      {
        if (K_FUNC::fEqual(R[j * n + j], 0., 1.e-8)==true)
        {
          ok=0;
          //printf("Cholesky: Matrix is not positive definite 2");
        }
        R[i * n + j] = (1.0 / R[j * n + j] * (A[i * n + j] - s));
      }
    }
 
  return ok;
}


//=============================================================================
/* Back substitution algorithm : Solve Ax=b
   IN : A: n x n upper triangular matrix
   IN : b: n     vector
   OUT: b: n     solution
   OUT: ok 1 : back substitution completed
           0: A is singular
*/
//=============================================================================

E_Int K_LINEAR::backSubstitution(E_Int n, E_Float* A, E_Float* b)
{  
  
  E_Int ok=1;

  for (E_Int i=n-1; i>=0; i--)
  {
    for (E_Int j = n-1; j>i; j--)
    {
      b[i] = b[i]-A[i+j*n]*b[j];
    }
    if(A[i+i*n]!=0) {b[i] = b[i]/A[i+i*n];}
    else 
    { 
      ok=0;
      //printf("BackSubstitution : Singular Matrix !\n");
    }
  }

  return ok;
}
