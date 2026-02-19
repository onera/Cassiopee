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
# include <vector>
# include "Interp/Interp.h"
# include "Linear/linear.h"
# include "Metric/metric.h"
#include "String/kstring.h"
using namespace std;
using namespace K_FLD;
using namespace K_ARRAY;

extern "C"
{
  void k6compsurfofstructcell_(E_Int& ni, E_Int& nj, E_Int& nk,
                                E_Int& indcell, E_Float* x,
                                E_Float* y, E_Float* z,
                                E_Float& surface);
}
//=============================================================================
/* On calcule le vecteur cfLoc tel que
  u_tilde(X) = cfLoc(X) Ui

  avec :
  cfLoc(X) = p(X).T A^-1 B

  A = P.TWP
  B = P.TW
  Matrice P | (P)ij = pj(xi)
  Matrice W | (W)ij = W(X, xi) delta_ij

  Pour diminuer le conditionnement on calcule les fonctions de la base
  au point 2*(xi-X)/radius

  IN: order: ordre de la formule
  IN: dimPb: dimension du probleme
  IN: pt: point ou l'on calcule les coefficients d'interpolation
  IN: xtDnr, ytDnr, ztDnr: coordonnees x, y, z du maillage donneur
  IN: dnrIndices: indices des points du stencil donneur
  IN: radius: longueurs des demi-axes de l'ellipse
  IN: axis: direction des demi-axes de l'ellipse
  IN: axisConst: direction(s) constante(s)
  OUT: cfloc: coefficients d'interpolation
  Retourne 1 : succes,
          -1 : pas de base de polynomes creee
          -2 : matrice singuliere
          -3 : matrice non SDP
          -4 : matrice B non creee */
//=============================================================================
E_Int K_INTERP::getInterpCoefMLS(E_Int order, E_Int dimPb, E_Int sizeBasis,
                                 E_Float* pt,
                                 E_Float* xtDnr, E_Float* ytDnr, E_Float* ztDnr,
                                 vector<E_Int>& dnrIndices,
                                 E_Float* radius, E_Float* axis, E_Int* axisConst,
                                 E_Float* cfLoc)
{
  // Taille de l'echantillon de donnees
  E_Int nDnrPts = dnrIndices.size();

  // declaration de
  vector<E_Float> B(sizeBasis*nDnrPts);
  vector<E_Float> A(sizeBasis*sizeBasis);
  vector<E_Int> rows; // i | pi(0) !=0
  rows.reserve(sizeBasis); // <---- Optimisation Xavier Juvigny
  E_Float zeros[3];
  E_Int ok;

  // Matrix B = P.TW | (P.TW)ij = pi(xj)*W(X,xj)
  ok = matrixB(nDnrPts, sizeBasis, dnrIndices, xtDnr, ytDnr, ztDnr, order, dimPb, pt, radius, axis, axisConst, &B[0]);
  if (ok == -1) return -4;

  // Matrix A = P.TW.P | (P.TW.P)ij = sum(k, 1, nDnrPts) pi(xk)*W(X,xk) pj(xk)
  //                      P.TW.P    = sum(k, 1, nDnrPts) p(xk) p(xk).T W(X,xk)
  ok = matrixA(nDnrPts, sizeBasis, dnrIndices, xtDnr, ytDnr, ztDnr, order, dimPb, pt, radius, axis, axisConst, &A[0]);
  if (ok == -1) return -3;

  // Find rows such as basis(0) != 0
  zeros[0] = 0.; zeros[1] = 0.; zeros[2] = 0.;

  // evaluate basis
  vector<E_Float> basis;
  ok = polyBasis(order, dimPb, sizeBasis, axisConst, zeros, basis);
  if ( ok == -1) return -1;

  for (E_Int i = 0; i < sizeBasis; i++)
    if (basis[i] != 0) rows.push_back(i);

  vector<E_Float> A_1(sizeBasis*rows.size());  //Only rows such as basis(0) is not null of A^-1 will be computed

  // Inversion of A: compute rows of A^-1 such as basis(0)[row]!=0
  if (order >= 2)
  {
    ok = K_LINEAR::inv(sizeBasis, &A[0], &A_1[0], rows);
    if (ok == 0)
    {
      // printf("Error: setInterpDataLS: singular matrix.\n");
      // for (E_Int i=0; i<sizeBasis ; i++)
      // {
      //   for (E_Int j=0; j<sizeBasis; j++)
      //   {
      //     printf(SF_F_ " ", A[i+j*sizeBasis]);
      //   }
      //   printf("\n");
      // }
      return -2;
    }
  }
  else A_1[0] = 1./A[0];

  // cfLoc = p(0).T A^-1 B = sum(j in rows) pj * (A^-1 B) [j][:]
  E_Int nrows = rows.size();
  for (E_Int j = 0; j < nDnrPts; j++)
  {
    cfLoc[j] = 0.;
    for (size_t i = 0; i < rows.size(); i++)
      for (E_Int k = 0; k < sizeBasis; k++)
      {
        cfLoc[j] += basis[rows[i]]*A_1[i+k*nrows]*B[k+j*sizeBasis];
      }
  }
  return 1;
}
//=============================================================================
/*  Matrix B = P.TW | (P.TW)ij = pi(xj)*W(X,xj)                              */
//=============================================================================
namespace {
    struct Vec3D{ E_Float crds[3]; };
}
E_Int K_INTERP::matrixB(
  E_Int nDnrPts, E_Int sizeBasis, vector<E_Int>& dnrIndices,
  E_Float *xtDnr, E_Float *ytDnr, E_Float *ztDnr,
  E_Int order, E_Int dimPb,
  E_Float *pt, E_Float *radius, E_Float *axis, E_Int *axisConst, E_Float *B)
{
  vector<Vec3D> dnrIn(nDnrPts);
  //vector<E_Float[3]> dnrIn(nDnrPts);//[3];
  E_Int ind, ok;

  for (E_Int ii = 0; ii < nDnrPts; ii++)
  {
    ind = dnrIndices[ii];
    dnrIn[ii].crds[0] = xtDnr[ind];
    dnrIn[ii].crds[1] = ytDnr[ind];
    dnrIn[ii].crds[2] = ztDnr[ind];
  }

  E_Float ptLoc[3];
  vector<E_Float> basisLoc;
  E_Float dx, dy, dz;

  E_Float r0 = 2./radius[0];
  E_Float r1 = 2./radius[1];
  E_Float r2 = 2./radius[2];
  for (E_Int i = 0; i < sizeBasis ; i++)
  {
    for (E_Int j = 0; j < nDnrPts; j++)
    {
      // ptLoc is dnrIn in the reference frame of the ellipse
      // Reduce the condition number
      dx = dnrIn[j].crds[0]-pt[0];
      dy = dnrIn[j].crds[1]-pt[1];
      dz = dnrIn[j].crds[2]-pt[2];

      ptLoc[0] = dx*axis[0] + dy*axis[1] + dz*axis[2];
      ptLoc[1] = dx*axis[3] + dy*axis[4] + dz*axis[5];
      ptLoc[2] = dx*axis[6] + dy*axis[7] + dz*axis[8];

      ptLoc[0] = ptLoc[0]*r0;
      ptLoc[1] = ptLoc[1]*r1;
      ptLoc[2] = ptLoc[2]*r2;
      ok = polyBasis(order, dimPb, sizeBasis, axisConst, ptLoc, basisLoc);
      if ( ok == -1) return -1;
      B[i+j*sizeBasis] = basisLoc[i] *weightFunction(pt,dnrIn[j].crds,radius,axis);
    }
  }
  return 1;
}
//=============================================================================
/*
  Matrix A = P.TWP | (P.TWP)ij = sum(k, 1, N) pi(xk)*W(X,xk)*pj(xk)
                      P.TWP    = sum(k, 1, N) p(xk) p(xk).T W(X,xk)
*/
//=============================================================================
E_Int K_INTERP::matrixA(E_Int nDnrPts, E_Int sizeBasis, vector<E_Int>& indicesIn,
                        E_Float *xtDnr, E_Float *ytDnr, E_Float *ztDnr,
                        E_Int order, E_Int dimPb,
                        E_Float *pt, E_Float *radius, E_Float *axis, E_Int *axisConst, E_Float *A)
{
  vector<Vec3D> dnrIn(nDnrPts);//[3];
  E_Int ind;
  E_Int ok = 1;
  vector<E_Float> basis;

  for (E_Int ii = 0; ii < nDnrPts; ii++)
  {
    ind = indicesIn[ii];
    dnrIn[ii].crds[0] = xtDnr[ind];
    dnrIn[ii].crds[1] = ytDnr[ind];
    dnrIn[ii].crds[2] = ztDnr[ind];
  }

  E_Float ptLoc[3];
  E_Float sum;
  E_Float r0 = 2./radius[0];
  E_Float r1 = 2./radius[1];
  E_Float r2 = 2./radius[2];
  E_Float dx, dy, dz;

  for (E_Int i = 0; i < sizeBasis ; i++)
  {
    for (E_Int j = 0; j < sizeBasis; j++)
    {
      sum = 0.;
      for (E_Int k = 0; k < nDnrPts; k++)
      {
        // ptLoc is dnrIn in the reference frame of the ellipse
        // Reduce the condition number
        dx = dnrIn[k].crds[0]-pt[0];
        dy = dnrIn[k].crds[1]-pt[1];
        dz = dnrIn[k].crds[2]-pt[2];
        ptLoc[0] = dx*axis[0] + dy*axis[1] + dz*axis[2];
        ptLoc[1] = dx*axis[3] + dy*axis[4] + dz*axis[5];
        ptLoc[2] = dx*axis[6] + dy*axis[7] + dz*axis[8];

        ptLoc[0] = ptLoc[0]*r0;
        ptLoc[1] = ptLoc[1]*r1;
        ptLoc[2] = ptLoc[2]*r2;
        ok = polyBasis(order, dimPb, sizeBasis, axisConst, ptLoc, basis);
        if ( ok != 1 ) return -1;
        sum += basis[i]*basis[j]*weightFunction(pt, dnrIn[k].crds, radius, axis);
      }
      A[i+j*sizeBasis] = sum;
    }
  }

  // Test si un element diagonal est nul : matrice non SDP
  // E_Int pb = 0;
  for (E_Int i = 0; i < sizeBasis; i++)
  {
    if (K_FUNC::E_abs(A[i+i*sizeBasis]) < 1.e-13) return -2;
  }
  return 1;
}
//=============================================================================
/* Construction de la base polynomiale et evaluation au point X
   IN: order: ordre maximal de la base
   IN: dimPb: dimension du maillage donneur
   IN: X: point ou on evalue la base
   basis est videe en debut de fonction
   retourne 1 si ok, -1 si echec */
//=============================================================================
E_Int K_INTERP::polyBasis(E_Int order, E_Int dimPb, E_Int sizeBasis,
                          E_Int *axisConst, E_Float* X,
                          vector<E_Float>& basis)
{
  basis.clear(); basis.push_back(1.);
  if (order >= 2)
  {
    for (E_Int dir = 0; dir < 3; dir++)
    {
      if (axisConst[dir] == 0) basis.push_back(X[dir]);
    }
  }

  // Formule d'interpolation d'ordre trois: base quadratique
  if (order >= 3)
  {
    for (E_Int dir = 0; dir < 3; dir++)
    {
      if (axisConst[dir] == 0) basis.push_back(X[dir]*X[dir]);

      // Termes croises
      for (E_Int dir2 = 0; dir2 < dir; dir2++)
      {
        if (axisConst[dir] == 0 && axisConst[dir2] == 0)
          basis.push_back(X[dir]*X[dir2]);
      }
    }
  }

  // Formule d'interpolation d'ordre quatre : base cubique
  if (order == 4)
  {
    for (E_Int dir = 0; dir < 3; dir++)
    {
      if (axisConst[dir] == 0) basis.push_back(X[dir]*X[dir]*X[dir]);

      // Termes croises
      for (E_Int dir2 = 0; dir2 < dir; dir2++)
      {
        if (axisConst[dir] == 0 && axisConst[dir2] == 0)
        {
          basis.push_back(X[dir]*X[dir2]*X[dir2]);
          basis.push_back(X[dir]*X[dir]*X[dir2]);
        }
      }
    }

    // Terme xyz en 3D
    if (dimPb == 3 && axisConst[0]==0 && axisConst[1]==0 && axisConst[2]==0)
      basis.push_back(X[0]*X[1]*X[2]);
  }

  E_Int basisSize = basis.size();
  if (basisSize != sizeBasis)
  {
    // printf("Error: K_INTERP:getInterpCoefsMLS: basis generation failed " SF_D3_ "\n",
    //        axisConst[0], axisConst[1], axisConst[2]);
    return -1;
  }
  return 1;
}
//=============================================================================
/* Fonction de poids pour la methode MLS

   w = 1 - 6*d**2 + 6*d**3   si     d<0.5
       2*(1-d)**3            si 0.5<d<1.
       0                     sinon
*/
//=============================================================================
E_Float K_INTERP::weightFunction(E_Float *pt, E_Float *xi,
                                    E_Float *radius, E_Float *axis)
{
  E_Float w = 0;
  E_Float t;
  // Calcul de la distance entre pt et xi divisee par r
  E_Float d = 0.;
  d = distEllipse(xi[0], xi[1], xi[2], pt, radius, axis);

  if (d <= 1./2.) { t = d*d; w = 1.-6.*t+6.*t*d; }
  else if (d > 0.5 && d < 1.) { t = (1.-d); w = 2.*t*t*t; }
  else if (d >= 1.) w = 0.;
  return w;
}
//=============================================================================
/* Equation de l'ellipse de centre center, d'axes de longueur
   radius le long des axes axis
   IN: x, y, z: point a considerer
   IN: center: centre de l'ellipse
   IN: radius: longueurs des demi-axes de l'ellipse
   IN: axis: axes de l'ellipse
   OUT: d: sqrt((X/a)**2+(Y/b)**2+(Z/c)**2) dans le bon repere
           si d==1 on est sur l'ellipse, d<1 on est dedans, d>1 on est dehors
*/
//=============================================================================
E_Float K_INTERP::distEllipse(E_Float x, E_Float y, E_Float z,
                              E_Float *center, E_Float *radius, E_Float *axis)

{
  E_Float d;
  E_Float X, Y, Z;
  E_Float cx, cy, cz;
  cx = x-center[0]; cy = y-center[1]; cz = z-center[2];
  // X Y Z are the coordinates in the reference frame of the ellipse
  X = cx*axis[0] + cy*axis[1] + cz*axis[2];
  Y = cx*axis[3] + cy*axis[4] + cz*axis[5];
  Z = cx*axis[6] + cy*axis[7] + cz*axis[8];
  X = X/radius[0];
  Y = Y/radius[1];
  Z = Z/radius[2];

  // Canonical equation of the ellipse
  d = X*X + Y*Y + Z*Z;
  return sqrt(d);
}
//=============================================================================
/* Construction du stencil carre autour du point donneur le plus proche
   du receveur (pour les maillages structures)
   IN: dimPb: dimension du maillage donneur
   IN: indi: indice du point donneur le plus proche du receveur
   IN: depth: demi-largeur du stencil
   IN: ni, nj, nk: dimension du maillage donneur
   OUT: indicesIn: vecteur des indices du stencil
*/
//=============================================================================
void K_INTERP::structStencil(E_Int dimPb, E_Int indi, E_Int depth,
                                E_Int ni, E_Int nj, E_Int nk,
                                vector<E_Int>& indicesIn)
{
  indicesIn.clear();

  // indi: indice global du sommet le plus proche de pt dans
  // la cellule contenant pt
  // indi = i_ + j_*ni + k_*ni*nj
  E_Int j_, k_; // conversion en indices locaux de indi
  E_Int nij = ni*nj;
  E_Int nijk = nij*nk;
  k_ = indi/nij;
  j_ = (indi % nij)/ni;
  //i_ = (indi % nij)%ni;

  indicesIn.push_back(indi); // On place de le sommet le plus proche en premier

  E_Int hj, hk;   // extensions suivant j et k
  hk = hj = depth;
  if (dimPb == 2) hk = 0; // On reste dans le meme plan en k
  else if (dimPb == 1) hk = hj = 0; // On reste dans le meme plan en k et sur la meme ligne j

  for (E_Int k = -hk; k <= hk; k++)
  {
    if (indi+nij*k >= 0 && indi+nij*k < nijk)
    {
      for (E_Int j = -hj; j <= hj; j++)
      {
        if (indi+ni*j >= k_*nij && indi+ni*j < (k_+1)*nij)   //(indi+ni*j)/(ni*nj)==k_   On reste dans le bon plan k_
          for (E_Int i = -depth; i <= depth; i++)
            if (indi+i >= j_*ni+k_*nij && indi+i < (j_+1)*ni+k_*nij && indi+i+ni*j+nij*k != indi)      //(indi+i)/ni ==j_+k_*nj   On reste sur la bonne ligne j_
            {
              if (indi+i+ni*j+nij*k < 0)
              {
                printf("indices negatif i=" SF_D_ " j=" SF_D_ " k=" SF_D_ "\n", i, j, k);
              }
              else indicesIn.push_back(indi+i+ni*j+nij*k);
            }
      }
    }
  }
  return;
}

//=============================================================================
/* Construction du stencil carre autour du point donneur le plus proche du receveur (maillage NGON)
   IN: indi: indice du point donneur le plus proche du receveur
   IN: depth: nombre de voisinages souhaites
   IN: cVN: connectivite sommet-sommet
   OUT: indIn: vecteur des indices du stencil
*/
//=============================================================================
void K_INTERP::NGONStencil(
  E_Int dimPb, E_Int indi,
  E_Int depth, vector< vector<E_Int> >& cVN,
  E_Float* xtDnr, E_Float* ytDnr, E_Float* ztDnr, vector<E_Int>& indicesIn)
{
  indicesIn.clear();
  E_Int ind, ind2, size;
  indicesIn.push_back(indi);  //point
  E_Float zref = ztDnr[indi];

  // premier voisinage
  size = cVN[indi].size();
  if (dimPb == 2) // ne garde que les pts de meme zref
  {
    for (E_Int i = 0; i < size; i++)
    {
      ind = cVN[indi][i]-1;
      if (K_FUNC::fEqualZero(ztDnr[ind] -zref) == true) indicesIn.push_back(ind);
    }
  }
  else // dim=3
  {
    for (E_Int i = 0; i < size; i++)
      indicesIn.push_back(cVN[indi][i]-1);
  }

  vector<E_Int> neighbour;
  E_Int sizeN1=indicesIn.size(); //taille du 1er voisinage
  E_Int sizeN2=0; //taille du 2eme voisinage
  E_Int sizeN3=0; //taille du 3eme voisinage
  E_Int sizeN4=0; //taille du 4eme voisinage

  if (depth >= 2)  //deuxieme voisinage
  {
    if (dimPb == 2)
    {
      for (E_Int n1 = 1; n1 < sizeN1; n1++)
      {
        ind = indicesIn[n1];
        for (size_t i = 0; i < cVN[ind].size(); i++)
        {
          ind2 = cVN[ind][i]-1;
          if (K_FUNC::fEqualZero(ztDnr[ind2]-zref) == true) neighbour.push_back(ind2);
        }
      }
    }
    else // dim=3
    {
      for (E_Int n1 = 1; n1 < sizeN1; n1++)
      {
        ind = indicesIn[n1];
        for (size_t i = 0; i < cVN[ind].size(); i++)
          neighbour.push_back(cVN[ind][i]-1);
      }
    }
    // On enleve les doublons et on determine la taille du deuxieme voisinage
    for (size_t ii = 0; ii < neighbour.size(); ii++)
    {
      if (find(indicesIn.begin(), indicesIn.end(), neighbour[ii]) == indicesIn.end())
      {
        indicesIn.push_back(neighbour[ii]);
        sizeN2++;
      }
    }
  }

  neighbour.clear();
  if (depth >= 3) //troisieme voisinage
  {
    if (dimPb == 2)
    {
      for (E_Int n2 = 0; n2 < sizeN2; n2++)
      {
        ind = indicesIn[n2+sizeN1];
        for (size_t i = 0; i < cVN[ind].size(); i++)
        {
          ind2 = cVN[ind][i]-1;
          if (K_FUNC::fEqualZero(ztDnr[ind2]-zref) == true) neighbour.push_back(ind2);
        }
      }
    }
    else // dim=3
    {
      for (E_Int n2 = 0; n2 < sizeN2; n2++)
      {
        ind = indicesIn[n2+sizeN1];
        for (size_t i = 0; i < cVN[ind].size(); i++)
          neighbour.push_back(cVN[ind][i]-1);
      }
    }

    // On enleve les doublons
    for (size_t ii = 0; ii < neighbour.size(); ii++)
    {
      if (find(indicesIn.begin(), indicesIn.end(), neighbour[ii]) == indicesIn.end())
      {
        indicesIn.push_back(neighbour[ii]);
        sizeN3++;
      }
    }
  }

  neighbour.clear();
  if (depth >= 4) //quatrieme voisinage
  {
    if (dimPb == 2)
    {
      for (E_Int n3 = 0; n3 < sizeN3; n3++)
      {
        ind = indicesIn[n3+sizeN2+sizeN1];
        for (size_t i = 0; i < cVN[ind].size(); i++)
        {
          ind2 = cVN[ind][i]-1;
          if (K_FUNC::fEqualZero(ztDnr[ind2]-zref) == true) neighbour.push_back(ind2);
        }
      }
    }
    else // dim=3
    {
      for (E_Int n3 = 0; n3 < sizeN3; n3++)
      {
        ind = indicesIn[n3+sizeN2+sizeN1];
        for (size_t i = 0; i < cVN[ind].size(); i++)
          neighbour.push_back(cVN[ind][i]-1);
      }
    }
    // On enleve les doublons
    for (size_t ii = 0; ii < neighbour.size(); ii++)
    {
      if (find(indicesIn.begin(), indicesIn.end(), neighbour[ii]) == indicesIn.end())
      {
        indicesIn.push_back(neighbour[ii]);
        sizeN4++;
      }
    }
  }

  return;
}

//=============================================================================
/*
   Determine les demi-axes du stencil elliptique grace a l'OBB du stencil
   carre
*/
//=============================================================================
void K_INTERP::findRadius(
  E_Int dimPb, vector<E_Int>& indicesIn, E_Int depth, E_Int order,
  E_Float *pt, E_Float *xtDnr, E_Float *ytDnr, E_Float *ztDnr,
  E_Float *axis, E_Float *radius, E_Int *axisConst)

{
  E_Int nDnrPts = indicesIn.size();
  E_Float bbox[6];

  vector<E_Float> xDnrIn(nDnrPts);
  vector<E_Float> yDnrIn(nDnrPts);
  vector<E_Float> zDnrIn(nDnrPts);
  E_Int ind;

  // Tableau des points du stencil
  for (E_Int ii = 0; ii < nDnrPts; ii++)
  {
    ind = indicesIn[ii];
    xDnrIn[ii] = xtDnr[ind];
    yDnrIn[ii] = ytDnr[ind];
    zDnrIn[ii] = ztDnr[ind];
  }

  // Oriented Bounding Box
  OBbox(dimPb, nDnrPts, &xDnrIn[0], &yDnrIn[0], &zDnrIn[0], axis, bbox);

  E_Float Pe1, Pe2, Pe3; // Projection du centre (pt) sur les axes de l'ellipse
  Pe1 = pt[0]*axis[0] + pt[1]*axis[1] + pt[2]*axis[2];
  Pe2 = pt[0]*axis[3] + pt[1]*axis[4] + pt[2]*axis[5];
  Pe3 = pt[0]*axis[6] + pt[1]*axis[7] + pt[2]*axis[8];

  // On met factor > 1.01 pour etre sur de prendre en compte les points extremes des axes de l 'ellipse
  E_Float factor=1.;
  if (order==2)
  {
    if (depth == 1)      factor=1.05;
    else if (depth == 2) factor=1.8;
    else factor=1.015;
  }
  else if (order == 3)
  {
    if (depth == 1)      factor=1.05;
    else if (depth == 2) factor=1.15;
    else factor=1.015;
  }
  else if (order == 4)
  {
    // CB: depth = 1?
    if (depth==2)      factor=1.015;
    else if (depth==3) factor=1.015;
    else if (depth>=4) factor=1.015;
  }

  radius[0] = factor*K_FUNC::E_max(bbox[3]-Pe1, Pe1-bbox[0]);
  radius[1] = factor*K_FUNC::E_max(bbox[4]-Pe2, Pe2-bbox[1]);
  radius[2] = factor*K_FUNC::E_max(bbox[5]-Pe3, Pe3-bbox[2]);
  //printf("radius ellipse %g %g %g\n", radius[0], radius[1], radius[2]);
  //printf("axis " SF_D3_ "\n", axisConst[0], axisConst[1], axisConst[2]);

  // Si le rayon est nul c'est que la direction de l'axe est constante
  // On met alors le rayon a 1 pour ne pas diviser par 0 ensuite
  for (E_Int i = 0; i < 3; i++)
  {
    /*
    if (K_FUNC::fEqual(radius[i], 0., 1.e-13) == true)
    {
      radius[i] = 1.;
      //axisConst[i] = 1;
    }
    //else axisConst[i] = 0;
    */
    // Tentative pour rendre le code plus portable
    if (K_FUNC::fEqual(radius[i], 0., 1.e-8) == true)
    {
      radius[i] = 1.e-8;
    }
  }

}

//=============================================================================
/*
  On ne garde que les points contenus dans l'ellipse de centre pt et d'axes axis
*/
//=============================================================================
void K_INTERP::indicesInEllipse(
  E_Float *pt, E_Float *radius, E_Float *axis,
  E_Float* xtDnr, E_Float* ytDnr, E_Float* ztDnr,
  vector<E_Int> &indicesIn)
{
  const E_Int nDnr = indicesIn.size();
  E_Int pos = 0;
  E_Int ind; E_Float d;
  for (E_Int ii = 0; ii < nDnr; ii++)
  {
    // Calcul de la distance entre pt et le donneur
    ind = indicesIn[pos];
    d = distEllipse2(
      xtDnr[ind], ytDnr[ind], ztDnr[ind],
      pt, radius, axis);

    if (d > 1.) indicesIn.erase(indicesIn.begin()+pos);
    else pos++;
  }
}

//=============================================================================
/* Enleve les points contenus dans le stencil ayant un cellN egal a 0 ou 2
   IN: cellNtDnr: cellN des points donneurs
   IN/OUT: indicesIn: indices des points du stencil
*/
//=============================================================================
void K_INTERP::eraseCellN0or2(E_Float *cellNtDnr, vector<E_Int> &indicesIn)
{
  const E_Int nDnr = indicesIn.size();

  E_Int pos = 0;
  for (E_Int ii = 0; ii < nDnr; ii++)
  {
    if (cellNtDnr[indicesIn[pos]] == 0 || cellNtDnr[indicesIn[pos]] == 2)
      indicesIn.erase(indicesIn.begin()+pos);
    else pos++;
  }
}

//=============================================================================
/* Enleve les points contenus dans le stencil ayant un cellN egal a 0
   IN: cellNtDnr: cellN des points donneurs
   IN/OUT: indicesIn: indices des points du stencil
*/
//=============================================================================
void K_INTERP::eraseCellN0(E_Float *cellNtDnr, vector<E_Int> &indicesIn)
{
  const E_Int nDnr = indicesIn.size();
  E_Int pos = 0;
  for (E_Int ii = 0; ii < nDnr; ii++)
  {
    if (cellNtDnr[indicesIn[pos]] == 0)
      indicesIn.erase(indicesIn.begin()+pos);
    else pos++;
  }
}

//=============================================================================
/* Enleve les points contenus dans le stencil a cause du cellN
   IN: cellNtDnr: cellN des points donneurs
   IN/OUT: indicesIn: indices des points du stencil
*/
//=============================================================================
void K_INTERP::eraseCellN(E_Int nature,
                             E_Float *cellNtDnr, vector<E_Int> &indicesIn)
{
  if (nature == 0) eraseCellN0(cellNtDnr, indicesIn);
  else eraseCellN0or2(cellNtDnr, indicesIn);
}

//=============================================================================
/* Equation de l'ellipse de centre center, d'axes de longueur
   radius le long des axes axis au carre
   IN: x, y, z: point a considerer
   IN: center: centre de l'ellipse
   IN: radius: longueurs des demi-axes de l'ellipse
   IN: axis: axes de l'ellipse
   OUT: d: (X/a)**2+(Y/b)**2+(Z/c)**2 dans le bon repere
   si d==1 on est sur l'ellipse, d<1 on est dedans, d>1 on est dehors
*/
//=============================================================================
E_Float K_INTERP::distEllipse2(
  E_Float x, E_Float y, E_Float z,
  E_Float *center, E_Float *radius, E_Float *axis)

{
  E_Float d;
  E_Float X, Y, Z;
  E_Float cx, cy, cz;
  cx = x-center[0]; cy = y-center[1]; cz = z-center[2];
  // X Y Z are the coordinates in the reference frame of the ellipse
  X = cx*axis[0] + cy*axis[1] + cz*axis[2];
  Y = cx*axis[3] + cy*axis[4] + cz*axis[5];
  Z = cx*axis[6] + cy*axis[7] + cz*axis[8];
  X = X/radius[0];
  Y = Y/radius[1];
  Z = Z/radius[2];

  // Canonical equation of the ellipse (au carre)
  d = X*X + Y*Y + Z*Z;
  return d;
}

//=============================================================================
/* Principal component analysis
   IN: n: size of the dataset
   IN: x, y, z: coordinates of the dataset
   OU: axis[8]: principal components: eigenvectors of the covariance matrix
*/
//=============================================================================
void K_INTERP::PCA(E_Int dimPb, E_Int n,
                      E_Float *x, E_Float *y, E_Float *z,
                      E_Float *axis)
{
  E_Float c00, c01, c11;       //covariance matrix
  E_Float c02, c12, c22;
  E_Float dx, dy, dz;
  E_Float mu[3];               //means
  E_Float lambda[3];           //eigenvalues
  E_Float e0[3], e1[3], e2[3]; //eigenvectors
  E_Float nfi = 1./n;
  E_Float nfi1 = 1./(n-1.);

  // Computation of the vector of means
  mu[0] = 0.; mu[1] = 0.; mu[2] = 0.;
  for(E_Int k = 0; k < n; k++)
  {
    mu[0] += x[k];
    mu[1] += y[k];
    mu[2] += z[k];
  }
  mu[0] *= nfi; mu[1] *= nfi; mu[2] *= nfi;

  // Computation of the covariance matrix
  c00 = 0.; c01 = 0.; c11 = 0.;
  c02 = 0.; c12 = 0.; c22 = 0.;
  for (E_Int k = 0; k < n; k++)
  {
    dx = x[k]-mu[0];
    dy = y[k]-mu[1];
    dz = z[k]-mu[2];
    c00 += dx*dx;
    c01 += dx*dy;
    c11 += dy*dy;
    c02 += dx*dz;
    c12 += dy*dz;
    c22 += dz*dz;
  }
  c00 *= nfi1; c01 *= nfi1; c11 *= nfi1;
  c02 *= nfi1; c12 *= nfi1; c22 *= nfi1;

  // Eigen values and vectors of covariance matrix
  if (dimPb == 2)
  {
    K_LINEAR::eigen2(c00, c11, c01,
		     lambda[0], lambda[1],
		     e0, e1);
    lambda[2] = 1.;
    e0[2] = 0.; e1[2] = 0.;
    e2[0] = 0.; e2[1] = 0.; e2[2] = 1.;
  }
  else
  {
    K_LINEAR::eigen3(c00, c01, c02,
		     c11, c12, c22,
		     lambda[0], lambda[1], lambda[2],
		     e0, e1, e2);
  }

  axis[0] = e0[0], axis[1] = e0[1], axis[2] = e0[2];
  axis[3] = e1[0], axis[4] = e1[1], axis[5] = e1[2];
  axis[6] = e2[0], axis[7] = e2[1], axis[8] = e2[2];
}

//=============================================================================
/* Oriented Bounding Box of a dataset
   IN: dimPb: dimension of the problem
   IN: n: size of the dataset
   IN: x, y, z: coordinates of the dataset
   OUT: axis: axis of the Oriented Bounding Box
   OUT: bbox: dimension of the OBB along each axis
*/
//=============================================================================
void K_INTERP::OBbox(
  E_Int dimPb, E_Int n, E_Float *x, E_Float *y, E_Float *z,
  E_Float *axis, E_Float *bbox)
{
  E_Float Pe1, Pe2, Pe3;   //Projection on the axis e1 e2 e3

  // On ne prend pas tous les points pour la PCA si le nuage est carre
  E_Int n_= n;
  E_Float val;
  if (dimPb == 2)
  {
    val = sqrt(n);
  }
  else if (dimPb == 3)
  {
    val = pow(n, 1./3.);
  }
  else val = E_Float(n);

  #ifdef E_ADOLC
  if (val == floor(val))   // si le nuage de point est carre on enleve une ligne (cas 2D)
    n_ = n-E_Int(val.value());
#else
  if (val == floor(val))   // si le nuage de point est carre on enleve une ligne (cas 2D)
    n_ = n-E_Int(val);
#endif

  // PCA of the stencil to find the axis
  PCA(dimPb, n_, x, y, z, axis);   // A MODIFIER ?
  //printf("axis ellipse 1: %g %g %g\n", axis[0], axis[1], axis[2]);
  //printf("axis ellipse 2: %g %g %g\n", axis[3], axis[4], axis[5]);
  //printf("axis ellipse 3: %g %g %g\n", axis[6], axis[7], axis[8]);

  bbox[0] = K_CONST::E_MAX_FLOAT; bbox[1] = K_CONST::E_MAX_FLOAT; bbox[2] = K_CONST::E_MAX_FLOAT;
  bbox[3] = -1.*K_CONST::E_MAX_FLOAT; bbox[4] = -1.*K_CONST::E_MAX_FLOAT; bbox[5] = -1.*K_CONST::E_MAX_FLOAT;

  for (E_Int k = 0; k < n; k++)
  {
    // Projection on the axis
    Pe1 = (x[k]*axis[0] + y[k]*axis[1] + z[k]*axis[2]);
    Pe2 = (x[k]*axis[3] + y[k]*axis[4] + z[k]*axis[5]);
    Pe3 = (x[k]*axis[6] + y[k]*axis[7] + z[k]*axis[8]);

    //Find extrema along each axis
    if (Pe1 < bbox[0])       bbox[0] = Pe1;      //minimun along axis 1
    else if (Pe1 > bbox[3])  bbox[3] = Pe1;      //maximun along axis 1

    if (Pe2 < bbox[1])       bbox[1] = Pe2;
    else if (Pe2 > bbox[4])  bbox[4] = Pe2;

    if (Pe3 < bbox[2])      bbox[2] = Pe3;
    else if (Pe3 > bbox[5]) bbox[5] = Pe3;
  }

}

//==============================================================================
E_Int K_INTERP::buildStencil(
  E_Int dimPb, E_Int order, E_Int nature, E_Int sizeBasis,
  E_Int center, E_Int resl, E_Int ni, E_Int nj, E_Int nk,
  vector< vector<E_Int> >& cVN,
  E_Float *xtDnr, E_Float *ytDnr, E_Float *ztDnr,
  E_Float *cellNtDnr, E_Float *pt, E_Float *radius,
  E_Float *axis, E_Int *axisConst, vector<E_Int>& indicesIn)
{
  indicesIn.clear();
  E_Int depth=0; //demi largeur du stencil

  vector<E_Float> A(sizeBasis*sizeBasis);
  vector<E_Float> A_1(sizeBasis);
  vector<E_Int> rows(1); rows[0] = 1;
  E_Int invOk = 1;  // etat de l'inversion de A

  E_Int tol = 2;
  if (order == 1 || order == 2) tol = 1;

  // Si on a pas assez de points, on augmente la demi-largeur du stencil
  while ((E_Int)indicesIn.size() <= sizeBasis+tol && depth < 10)
  {
    depth++;
    indicesIn.clear();
    if (resl == 1) // structure
      structStencil(dimPb, center, depth, ni, nj, nk, indicesIn);
    else if (resl == 2) // non structure
      NGONStencil(dimPb, center, depth, cVN, xtDnr, ytDnr, ztDnr, indicesIn);

    if ((E_Int)indicesIn.size() > sizeBasis+tol) // OK
    {
      // OBB et determination des demi axes de l'ellipse
      findRadius(dimPb, indicesIn, depth, order, pt,
                 xtDnr, ytDnr, ztDnr, axis, radius, axisConst);
      // Stencil elliptique
      indicesInEllipse(pt, radius, axis, xtDnr, ytDnr, ztDnr, indicesIn);
      // Attention: prendre en compte la valeur du cellN des points dans la sphere...
      eraseCellN(nature, cellNtDnr, indicesIn);
    }
  }

  if (depth >= 10) return -1;

  // On verifie que la matrice A du MLS est inversible
  // Si elle n'est pas inversible on augmente la taille du stencil
  K_INTERP::matrixA(indicesIn.size(), sizeBasis, indicesIn,
                    xtDnr, ytDnr, ztDnr,
                    order, dimPb,
                    pt, radius, axis, axisConst, &A[0]);
  invOk = K_LINEAR::inv(sizeBasis, &A[0], &A_1[0], rows);
  while (invOk == 0 && depth < 10)
  {
    depth++;
    indicesIn.clear();
    if (resl == 1)
      structStencil(dimPb, center, depth, ni, nj, nk, indicesIn);
    else if (resl == 2)
      NGONStencil(dimPb, center, depth, cVN, xtDnr, ytDnr, ztDnr, indicesIn);

    findRadius(dimPb, indicesIn, depth, order, pt,
               xtDnr, ytDnr, ztDnr, axis, radius, axisConst);
    indicesInEllipse(pt, radius, axis, xtDnr, ytDnr, ztDnr, indicesIn);
    // Attention: prendre en compte la valeur du cellN des points dans la sphere...
    eraseCellN(nature, cellNtDnr, indicesIn);
    K_INTERP::matrixA(indicesIn.size(), sizeBasis, indicesIn,
                      xtDnr, ytDnr, ztDnr,
                      order, dimPb,
                      pt, radius, axis, axisConst, &A[0]);
    invOk = K_LINEAR::inv(sizeBasis, &A[0], &A_1[0], rows);
  }

  if (depth >= 10) return -2;

  return 1;
}

//=============================================================================
/*
   Recherche du meilleur donneur
   On ne selectionne que les zones qui sont "proches" du point receveur et on compare leur volume
   On peut penaliser les cellules qui sont entourees de pt de cellN!=1

   IN: pt: point receveur
   IN: nDnrZones: nombre de zones donneuses
   IN: fields: champs des zones donneuses
   IN: vectOfKdTrees:  KDTree des points de cellN1 des zones donneuses
   IN: resl: 1 structure, 2 NGON
   IN: a2, a3, a4: //ni,nj,nk ou cnt en NS
   IN: vectOfCorres: correspondance entre les indices dans le KDT et les indices globaux
   IN: vectOfcVF: connectivites Vertex/Faces
   IN: vectOfcEV: connectivites Element/noeuds
   IN: vectOfcVN: connectivites Vertex/Vertex voisins
   IN: vectOfcFE: connectivites Faces/Elts
   IN: vectOfPosElt: tableaux de position des elements dans la connectivite
   IN: vectOfPosFace: tableaux de position des faces dans la connectivite
   IN: vectOfDimElt: tableaux de la dimension des elements
   IN: depth: profondeur du stencil donneur
   IN: penalty: 0 pas de penalisation
                1 on penalise les points qui ont beaucoup de voisins avec un cellN != 1
   OUT: noBest: indice du meilleur donneur
   OUT: indBestDnr: indice du point le plus proche de pt dans la meilleur zone donneuse
*/
//=============================================================================
void K_INTERP::getBestDonor(
  E_Int dimPb, E_Float *pt, E_Int nDnrZones, vector<FldArrayF*>& fields,
  vector<K_SEARCH::KdTree<FldArrayF>*>& vectOfKdTrees,
  vector<E_Int>& resl,
  vector<void*>& a2, vector<void*>& a3, vector<void*>& a4,
  vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs,
  vector<E_Int>& poscs, vector<FldArrayI*>& vectOfCorres,
  vector<vector<vector<E_Int> > >& vectOfcVF,
  vector<vector<vector<E_Int> > >& vectOfcEV,
  vector<vector<vector<E_Int> > >& vectOfcVN,
  vector<FldArrayI>& vectOfcFE,
  vector<FldArrayI>& vectOfPosElt,
  vector<FldArrayI>& vectOfPosFace,
  vector<FldArrayI>& vectOfDimElt,
  E_Int depth, E_Int penalty,
  E_Int& noBest, E_Int& indBestDnr, E_Float& dBestDnr)
{

  E_Int axisConst[3]; // =1 si l'axe est une direction constante, 0 sinon
  // Taille du probleme et direction constante
  // si axisConst=1, direction par prise en compte dans la base
  axisConst[0] = axisConst[1] = axisConst[2] = 0;

  if (dimPb == 3)        // pas de direction constante et probleme 3D
  { axisConst[0] = 0; axisConst[1] = 0; axisConst[2] = 0;}
  else if (dimPb == 2)   // pas de direction constante mais probleme 2D
  {
    axisConst[2] = 1; // pas de z
  }
  else if (dimPb == 1)
  {
    axisConst[1] = 1; axisConst[2] = 1; // pas de y/z
  }

  E_Float depth2s2 = depth*depth*1.41;
  E_Float volMin = K_CONST::E_MAX_FLOAT;
  noBest = -1;
  indBestDnr = -1;

  vector<E_Int> indCandidates;
  vector<E_Int> blkCandidates;
  vector<E_Float> volCandidates;
  E_Int pen;
  E_Int inddummy = -1;  // a laisser absolument a -1 pour K_METRIC::compVolOfStructCell3D


  for (E_Int no = 0; no < nDnrZones; no++)
  {
    // donnees du donneur
    E_Float* xtDnr = fields[no]->begin(posxs[no]);
    E_Float* ytDnr = fields[no]->begin(posys[no]);
    E_Float* ztDnr = fields[no]->begin(poszs[no]);
    E_Float* cellN = fields[no]->begin(poscs[no]);

    // distance entre le point le plus proche et le receveur
    E_Float dDnr;
    E_Int indDnr, temp;
    E_Int* corres = vectOfCorres[no]->begin();
    indDnr = vectOfKdTrees[no]->getClosest(pt, dDnr);
    indDnr = corres[indDnr];
    temp = indDnr;   // la fonction K_METRIC::compMeanLengthOfStructCell modifie l'indice de la cellule...
    //(pas K_METRIC::compVolOfStructCell3D)

    E_Float vol = 0.; //volume de la cellule donneuse eventuellement penalise
    E_Float dRef=0.; //distance maximale pour laquelle on compare les volumes
    E_Float realVol=0;// volume de la cellule donneuse effectif
    if (resl[no] == 1) //cas structure
    {
      // dimensions du maillage
      E_Int ni = *(E_Int*)a2[no];
      E_Int nj = *(E_Int*)a3[no];
      E_Int nk = *(E_Int*)a4[no];
      // Calcul du volume de la cellule donneuse
      if ((ni==1)||(nj==1)||(nk==1))   // 2D
        K_METRIC::compSurfOfStructCell(
          ni, nj, nk, temp,
          xtDnr, ytDnr, ztDnr, vol);
      else if ((ni>1)&&(nj>1)&&(nk>1)) //3D
        K_METRIC::compVolOfStructCell3D(
          ni, nj, nk, inddummy, temp,
          xtDnr, ytDnr, ztDnr, vol);
      realVol = vol;
      // On penalise la cellule si elle a des voisins de cellN != 1 et si elle est sur le bord
      if (penalty == 1)
      {
        pen = penalizeBorderStruct(dimPb, indDnr, depth, cellN, ni, nj, nk);
        vol += pen*1.e3;
      }

      // taille de la cellule du donneur
      dRef = compMaxLengthOf1stStructNbr(dimPb, indDnr, xtDnr, ytDnr, ztDnr, ni, nj, nk);
    }
    else if (resl[no] == 2) //cas NGON
    {
      // connectivite
      FldArrayI& cNG = *(FldArrayI*)a2[no];

      // Calcul du volume de la cellule donneuse
      // on moyenne les volumes des elements qui contiennent le point
      vector<E_Int> ENbrs; //elements qui contiennent le point
      connectNG2ENbrs(indDnr, vectOfcVF[no], vectOfcFE[no], ENbrs);
      for (size_t i = 0; i < ENbrs.size(); i++)
      {
        E_Float volElt = 0.; 
        K_METRIC::compNGonVolOfElement(
          xtDnr, ytDnr, ztDnr, cNG,
          ENbrs[i], vectOfcEV[no], vectOfPosElt[no],
          vectOfPosFace[no], vectOfDimElt[no],
          volElt);
        vol += volElt;
      }
      vol /= float(ENbrs.size());
      realVol = vol;
      // On penalise la cellule si elle a des voisins de cellN != 1 et si elle est sur le bord
      if (penalty == 1)
      {
        pen = penalizeBorderNGON(dimPb, indDnr, depth, cellN,
                                 xtDnr, ytDnr, ztDnr, vectOfcVN[no],
                                 vectOfcVF[no], vectOfcFE[no]);
        vol += pen*1.e3;
      }
      dRef = compMaxLengthOf1stNgonNbr(indDnr, xtDnr, ytDnr, ztDnr, vectOfcVN[no]);
    }

    // Selection du meilleur donneur
    if (dDnr < depth2s2*dRef)
    {
      E_Int order = 3; E_Int nature = 1;
      E_Int sizeBasis = K_FUNC::fact(dimPb+order-1)*1./(K_FUNC::fact(dimPb)*K_FUNC::fact(order-1));
      E_Int ni = 0, nj = 0, nk = 0;
      if (resl[no] == 1)
      {
        ni=*(E_Int*)a2[no];
        nj=*(E_Int*)a3[no];
        nk=*(E_Int*)a4[no];
      }
      if (dDnr < K_CONST::E_GEOM_CUTOFF)//pt coincident
      {
        indCandidates.push_back(indDnr);
        blkCandidates.push_back(no);
        volCandidates.push_back(realVol);
        dBestDnr = dDnr;
        indBestDnr = indDnr;
        noBest = no;
        volMin = vol;
      }
      else 
      {
        E_Float radius[3]; // longueurs des axes de l'ellipse
        radius[0] = radius[1] = radius[2] = 1.;
        E_Float axis[9];  // Axes de l'ellipse
        vector<E_Int> indicesIn;
        E_Int ok = buildStencil(dimPb, order, nature, sizeBasis,
                                indDnr, resl[no], ni, nj, nk, vectOfcVN[no],
                                xtDnr, ytDnr, ztDnr, cellN,
                                pt, radius, axis, axisConst, indicesIn);

        if (vol < volMin && ok == 1)
        {
          volMin = vol;    // volume minimal
          noBest = no;     // numero de la zone donneuse
          indBestDnr = indDnr; // indice du meilleur donneur
          dBestDnr = dDnr;
        }
      }
    }
  }

  // si points coincidents, on definit celui de plus petit volume comme donneur
  volMin = K_CONST::E_MAX_FLOAT;
  for (size_t i = 0; i < indCandidates.size(); i++)
  {
    if (volCandidates[i] < volMin)
    {
      volMin = volCandidates[i];
      indBestDnr = indCandidates[i];
      noBest = blkCandidates[i];
      dBestDnr = 0.;
    }
  }
  dBestDnr = sqrt(dBestDnr);
  return;
}

//=============================================================================
/*
  Determine si le sommet d'indice indV est sur le bord pour un maillage NGON conforme

  IN: indV: indice du point considere
  IN: cVF: connectivite Vertex/Faces
  IN: cFE: connectivite Faces/Elts
  OUT: onEdge: nombre de fois ou on est sur un bord
*/
//=============================================================================
E_Int K_INTERP::isOnEgdeNGON(E_Int indV, vector< vector<E_Int> >& cVF, FldArrayI& cFE)
{
  E_Int onEdge = 0;

  // cFE est de la forme (nfaces, 2)
  E_Int* facep1 = cFE.begin(1);  //element d'un cote de la face
  E_Int* facep2 = cFE.begin(2);  //element de l'autre cote de la face (0 si sur le bord)
  E_Int face;

  for (unsigned i = 0; i < cVF[indV].size(); i++)
  {
    face = cVF[indV][i]-1;
    if (facep1[face]==0 || facep2[face]==0) onEdge += 1;
  }

  return onEdge;
}

//=============================================================================
/*
  On penalise les cellules sur le bord ou qui ont des voisins de cellN!=1
  (maillage NGON conforme)

  IN: indV: indice du point considere
  IN: depth: profondeur du stencil du donneur (demi-largeur du voisinage)
  IN: cellN: cellN des points du maillage
  IN: cVN, cVF, cFE: connectivite Vertex/Vertex voisins, Vertex/Faces, Faces/Elts
  OUT: penalty: nombre de voisins qui possedent un cellN different de 1
              + le nombre de fois ou on est sur un bord
*/
//=============================================================================
E_Int K_INTERP::penalizeBorderNGON(
  E_Int dimPb, E_Int indV, E_Int depth, E_Float *cellN,
  E_Float* xtDnr, E_Float* ytDnr, E_Float* ztDnr,
  vector< vector<E_Int> >& cVN, vector< vector<E_Int> >& cVF,
  FldArrayI& cFE)
{
  E_Int penalty = 0;

  // voisinage de largeur depth
  vector<E_Int> nbr;
  NGONStencil(dimPb, indV, depth, cVN, xtDnr, ytDnr, ztDnr, nbr);

  for (unsigned int i = 0; i < nbr.size(); i++)
  {
    // On penalise si on a des voisins de cellN !=1
    if (cellN[nbr[i]] != 1) penalty += 1;

    // On penalise si on a des voisins sur le bord
    penalty += isOnEgdeNGON(indV, cVF, cFE);
  }

  return penalty;
}

//=============================================================================
/*
  On penalise les cellules sur le bord ou qui ont des voisins de cellN!=1 (maillage Structure)

  IN: indV: indice du point considere
  IN: depth: taille du recouvrement (demi-largeur du voisinage)
  IN: cellN: cellN des points du maillage
  IN: ni, nj, nk: dimensions du maillage
  OUT: penalty: nombre de voisins qui possedent un cellN different de 1
              + le nombre de fois ou on est sur un bord
*/
//=============================================================================
E_Int K_INTERP::penalizeBorderStruct(
  E_Int dimPb, E_Int indV, E_Int depth, E_Float *cellN,
  E_Int ni, E_Int nj, E_Int nk)
{
  E_Int penalty = 0;

  // voisinage de largeur depth
  vector<E_Int> nbr;
  structStencil(dimPb, indV, depth, ni, nj, nk, nbr);

  for (unsigned int i = 0; i < nbr.size(); i++)
  {
    // On penalise si on a des voisins de cellN !=1
    if (cellN[nbr[i]] != 1) penalty += 1;

    // On penalise si on a des voisins sur le bord
    penalty += isOnEgdeStruct(indV, ni, nj, nk);
  }
  return penalty;
}

//=============================================================================
/*
  Determine si un point est sur un bord pour un maillage Structure
  IN: ind: indice du point considere
  IN: ni, nj, nk: dimensions du maillage
  OUT: onEdge: nombre de fois ou on est sur un bord
*/
//=============================================================================
E_Int K_INTERP::isOnEgdeStruct(E_Int ind, E_Int ni, E_Int nj, E_Int nk)
{
  E_Int onEdge=0;

  // ind = i_ + j_*ni + k_*ni*nj
  E_Int i_, j_, k_;
  E_Int nij = ni*nj;
  i_ = (ind % nij)%ni;
  j_ = (ind % nij)/ni;
  k_ = ind/nij;

  if (i_==0 || i_==ni) onEdge += 1;
  if (j_==0 || j_==nj) onEdge += 1;
  if (k_==0 || k_==nk) onEdge += 1;

  return onEdge;
}

//=============================================================================
/*
  Determine la distance maximale entre les points du stencil et le centre
  maillage NGon
  IN: ind: indice du centre du stencil
  IN: x, y, z: coordonnees x, y, z du maillage
  IN: cVN: connectivite Vertex/Vertex voisins
  OUT: retourne la distance max au carre
*/
//=============================================================================
E_Float K_INTERP::compMaxLengthOf1stNgonNbr(
  E_Int ind, E_Float* x, E_Float* y, E_Float* z,
  vector<vector<E_Int> >& cVN)
{
  E_Float dMax = -1.*K_CONST::E_MAX_FLOAT;
  E_Int indV;
  E_Float d, dx, dy, dz;

  for (unsigned int v = 0; v < cVN[ind].size(); v++)
  {
    indV = cVN[ind][v]-1;
    dx = x[ind]-x[indV]; dy = y[ind]-y[indV]; dz = z[ind]-z[indV];
    d = dx*dx+dy*dy+dz*dz;
    if (d > dMax) dMax = d;
  }
  return dMax;
}

//=============================================================================
/*
  Determine la distance maximale entre les points du stencil et le centre
  maillage structure
  IN: ind: indice du centre du stencil
  IN: x, y, z: coordonnees x, y, z du maillage
  IN: ni, nj, nk: dimensions du maillage
  OUT: retourne la distance max au carre
*/
//=============================================================================
E_Float K_INTERP::compMaxLengthOf1stStructNbr(
  E_Int dimPb, E_Int ind, E_Float* x, E_Float* y, E_Float* z,
  E_Int ni, E_Int nj, E_Int nk)
{
  E_Float dMax = -1.*K_CONST::E_MAX_FLOAT;
  E_Int indV; E_Float d, dx, dy, dz;

  // premier voisinage
  vector<E_Int> nbr;
  structStencil(dimPb, ind, 1, ni, nj, nk, nbr);

  for (unsigned int v = 0; v < nbr.size(); v++)
  {
    indV = nbr[v];
    dx = x[ind]-x[indV]; dy = y[ind]-y[indV]; dz = z[ind]-z[indV];
    d = dx*dx+dy*dy+dz*dz;
    if (d > dMax) dMax = d;
  }
  return dMax;
}

//=============================================================================
/*
  Calcule les elements contenants le sommet d'indice indV pour un NGON conforme.
  IN: indV: indice du point considere
  IN: cVF: connectivite Vertex/Faces
  IN: cFE: connectivite Faces/Elts
  OUT: ENbrs: indices des elements contenant le sommet indV
*/
//=============================================================================
void K_INTERP::connectNG2ENbrs(E_Int indV, vector< vector<E_Int> >& cVF, FldArrayI& cFE,
                               vector<E_Int>& ENbrs)
{
  // cFE est de la forme (nfaces, 2)
  E_Int* facep1 = cFE.begin(1);  //element d'un cote de la face
  E_Int* facep2 = cFE.begin(2);  //element de l'autre cote de la face (0 si sur le bord)

  for (size_t i = 0; i < cVF[indV].size(); i++)
  {
    E_Int face = cVF[indV][i]-1;
    if (facep1[face] != 0) ENbrs.push_back(facep1[face]-1);
    if (facep2[face] != 0) ENbrs.push_back(facep2[face]-1);
  }
  sort(ENbrs.begin(), ENbrs.end());
  ENbrs.erase(unique(ENbrs.begin(), ENbrs.end()), ENbrs.end());
}
