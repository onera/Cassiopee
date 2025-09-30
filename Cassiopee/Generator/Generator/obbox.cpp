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

# include "generator.h"
# include <iostream>     // for printing cout
# include <math.h>       // for the maths
# include <stdio.h>      // for printf function

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Returns the oriented bounding box of an array as an array */
//=============================================================================
PyObject* K_GENERATOR::obbox(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int Weighting;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &Weighting))return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, 
                                     ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "obbox: invalid array.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "obbox: cannot find coordinates in array.");
    return NULL;        
  }
  posx++; posy++; posz++;
  
  if (Weighting == 1) //Input zone must be a TRI mesh
  {
    if (res == 1) //structured
    {
      PyErr_SetString(PyExc_TypeError,"obbox: 1st arg must be a TRI zone.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
    else
    {
      if (strcmp(eltType,"TRI") != 0)
      {
	PyErr_SetString(PyExc_TypeError,"obbox: 1st arg must be a TRI zone.");
	RELEASESHAREDB(res, array, f, cn); return NULL;  
      }  
    }     
  }

  // Extracts coordinates list
  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);
  E_Int nt = f->getSize();
  E_Int api = f->getApi();
  
  // Initializes the mean point vector and covariance matrix
  FldArrayF mu(3); mu.setAllValuesAtNull();
  FldArrayF Cov(3,3); Cov.setAllValuesAtNull();
  if (Weighting == 0)// no weighting by triangle surfaces : cloud used
  {
    // METHOD 1: point cloud (Weighting=0)
    // Calculation of mu
    for (E_Int i=0; i < nt; i++)
    {
      mu[0] = mu[0] + xt[i];
      mu[1] = mu[1] + yt[i];
      mu[2] = mu[2] + zt[i];
    }
    mu[0] = mu[0]/nt;
    mu[1] = mu[1]/nt;
    mu[2] = mu[2]/nt;
    // Calculation of the covariance matrix
    for (E_Int k=1; k < 4; k++)
    {
      for (E_Int j=0; j < 3; j++)
      {
        for (E_Int i=0; i < nt; i++)
        {
          Cov(j,k) = Cov(j,k) + ((*f)(i,j+1)-mu[j])*((*f)(i,k)-mu[k-1]);
        }
      }
    }
    for (E_Int k=1; k < 4; k++)
      {
      for (E_Int j=0; j < 3; j++)
      {
        Cov(j,k) = Cov(j,k)/nt;
      }
    }
  }
  else
  {
  // METHOD 2: triangular dense sampling (for unstructured surfaces of TRI)
  E_Int nelts = cn->getSize();
  // Construction of vertex coordinates vectors: p q r
  // as well as the vectors qp(=q-p) and rp(=r-p)
  FldArrayF p(3,nelts);
  FldArrayF q(3,nelts);
  FldArrayF r(3,nelts);
  FldArrayF qp(3,nelts);
  FldArrayF rp(3,nelts);
  for (E_Int e=0; e < nelts; e++)
  {
      for (E_Int i=1; i < 4; i++)
      { 
        p(i-1,e+1) = (*f)((*cn)(e,1)-1,i);
        q(i-1,e+1) = (*f)((*cn)(e,2)-1,i);
        r(i-1,e+1) = (*f)((*cn)(e,3)-1,i);
        qp(i-1,e+1) = q(i-1,e+1)-p(i-1,e+1);
        rp(i-1,e+1) = r(i-1,e+1)-p(i-1,e+1);
      }
  } 
  // Calculation of the triangle area
  FldArrayF m(nelts);
  E_Float mTot = 0.;
  for (E_Int e=0; e < nelts; e++)
  {
      m[e] = 0.5*sqrt( pow(qp(1,e+1)*rp(2,e+1)-rp(1,e+1)*qp(2,e+1),2) +
                       pow(qp(0,e+1)*rp(2,e+1)-rp(0,e+1)*qp(2,e+1),2) +
                       pow(qp(0,e+1)*rp(1,e+1)-rp(0,e+1)*qp(1,e+1),2) );
      mTot = mTot + m[e];
  }
  // Calculation of mu
  for (E_Int e=0; e < nelts; e++)
  {
    mu[0] = mu[0] + (p(0,e+1)+q(0,e+1)+r(0,e+1))*m[e];
    mu[1] = mu[1] + (p(1,e+1)+q(1,e+1)+r(1,e+1))*m[e];
    mu[2] = mu[2] + (p(2,e+1)+q(2,e+1)+r(2,e+1))*m[e];
  }
  mu[0] = mu[0]/(3.*mTot);
  mu[1] = mu[1]/(3.*mTot);
  mu[2] = mu[2]/(3.*mTot);
  // Calculates the covariance matrix
  for (E_Int k=1; k < 4; k++)
  {
     for (E_Int j=0; j < 3; j++)
     {
       for (E_Int e=0; e < nelts; e++)
       {
        Cov(j,k) = Cov(j,k) + m[e]*((p(j,e+1)+q(j,e+1)+r(j,e+1) - 3.*mu[j])*
                                    (p(k-1,e+1)+q(k-1,e+1)+r(k-1,e+1)-3.*mu[k-1]) +
                 (p(j,e+1)-mu[j])*(p(k-1,e+1)-mu[k-1]) +
                 (q(j,e+1)-mu[j])*(q(k-1,e+1)-mu[k-1]) + 
                 (r(j,e+1)-mu[j])*(r(k-1,e+1)-mu[k-1]));
       }
     }
  }
  for (E_Int k=1; k < 4; k++)
  {
      for (E_Int j=0; j < 3; j++)
      {
        Cov(j,k) = Cov(j,k)/(24.0*nelts);
      }
   }
  }
  /*
  cout << "Mu calculated and is:" << endl;
  cout << mu[0] << ", " << mu[1] << ", " << mu[2] << endl;  
  cout << "Cov calculated and is:" << endl;
  cout << Cov(0,1) << ", " << Cov(1,2) << ", " << Cov(2,3) << ", ";
  cout << Cov(0,2) << ", " << Cov(0,3) << ", " << Cov(1,3) << endl;  
  */
/*
  printf("Cov = [%16g, %16g, %16g;\n %16g, %16g, %16g;\n %16g, %16g, %16g];\n",
                                              Cov(0,1),Cov(0,2),Cov(0,3),
                                              Cov(1,1),Cov(1,2),Cov(1,3),
                                              Cov(2,1),Cov(2,2),Cov(2,3));
*/
  // Calculates the main axes of the OBB
  FldArrayF v0(3);
  FldArrayF v1(3);
  FldArrayF v2(3);
  FldArrayF v0box(3);
  FldArrayF v1box(3);
  FldArrayF v2box(3);
  FldArrayF SmallEdges(2);
  FldArrayI SmallEdgesInd(2); SmallEdgesInd.setAllValuesAt(-1);
  E_Float lambda0;
  E_Float lambda1;
  E_Float lambda2;  
  
  /* eigen3 does not work for the moment. Using eigen3bis instead
  E_Int rank = K_LINEAR::eigen3(Cov(0,1), Cov(0,2), Cov(0,3),
               		        Cov(1,2), Cov(1,3), Cov(2,3),
               		        lambda0, lambda1, lambda2, 
               		        v0.begin(), v1.begin(), v2.begin());
  */

  E_Float A[3][3];
  E_Float V[3][3];
  E_Float d[3];
  for (E_Int k=1; k < 4; k++)
  {
    for (E_Int j=0; j < 3; j++)
    {
      A[j][k-1] = Cov(j,k);
    }
  }
  K_LINEAR::eigen3bis(A, V, d);

  for (E_Int i=0; i < 3; i++)
  {
    v0[i] = V[i][0];
    v1[i] = V[i][1];
    v2[i] = V[i][2];
  }
  lambda0 = d[0];
  lambda1 = d[1];
  lambda2 = d[2];
  //printf("lambda = [%g, %g, %g];\n", lambda0, lambda1, lambda2);
  /*
  cout << "Eigenvectors and Eigenvalues calculated" << endl;
  printf("lambda = [%g, %g, %g];\n", lambda0, lambda1, lambda2);
  printf("v0 = [%g, %g, %g];\n", v0[0], v0[1], v0[2]);
  printf("v1 = [%g, %g, %g];\n", v1[0], v1[1], v1[2]);
  printf("v2 = [%g, %g, %g];\n", v2[0], v2[1], v2[2]);  
  */
  
  E_Float VolumeAABB = 0;
  E_Int isBoxShaped = 0;
  FldArrayF vmax(3);
  FldArrayF vmin(3);
  FldArrayI imax(3); imax.setAllValuesAt(-1);
  FldArrayI imin(3); imin.setAllValuesAt(-1);
  vmin[0]=xt[0];
  vmax[0]=xt[0];  
  vmin[1]=yt[0];
  vmax[1]=yt[0]; 
  vmin[2]=zt[0];
  vmax[2]=zt[0];  
  E_Float Tol=1.e-12;
if ( (K_FUNC::E_abs(lambda0-lambda1)<Tol)
   || (K_FUNC::E_abs(lambda0-lambda2)<Tol)
   || (K_FUNC::E_abs(lambda1-lambda2)<Tol) )
  {
  isBoxShaped = 1;
  //cout << "IT IS A BOX!!" << endl;
  // It is a box. Recalculates main axes.
  for (E_Int i=0; i < nt; i++)
  {
     if ( (xt[i] >= vmax[0])
        && (i!=imax[0])
        && (i!=imax[1])
        && (i!=imax[2])
        && (i!=imin[0])
        && (i!=imin[1])
        && (i!=imin[2]) )  {vmax[0]=xt[i]; imax[0]=i;}
  }
  for (E_Int i=0; i < nt; i++)
  {
     if ( (xt[i] <= vmin[0])
        && (i!=imax[0])
        && (i!=imax[1])
        && (i!=imax[2])
        && (i!=imin[0])
        && (i!=imin[1])
        && (i!=imin[2]) )  {vmin[0]=xt[i]; imin[0]=i;}
  }
  for (E_Int i=0; i < nt; i++)
  {
     if ( (yt[i] >= vmax[1])
        && (i!=imax[0])
        && (i!=imax[1])
        && (i!=imax[2])
        && (i!=imin[0])
        && (i!=imin[1])
        && (i!=imin[2]) )  {vmax[1]=yt[i]; imax[1]=i;}
  }

  for (E_Int i=0; i < nt; i++)
  {
     if ( (yt[i] <= vmin[1])
        && (i!=imax[0])
        && (i!=imax[1])
        && (i!=imax[2])
        && (i!=imin[0])
        && (i!=imin[1])
        && (i!=imin[2]) ) {vmin[1]=yt[i]; imin[1]=i;}
  }

  for (E_Int i=0; i < nt; i++)
  {
     if ( (zt[i] >= vmax[2])
        && (i!=imax[0])
        && (i!=imax[1])
        && (i!=imax[2])
        && (i!=imin[0])
        && (i!=imin[1])
        && (i!=imin[2]) )  {vmax[2]=zt[i]; imax[2]=i;}
  }
  for (E_Int i=0; i < nt; i++)
  {
     if ( (zt[i] <= vmin[2])
        && (i!=imax[0])
        && (i!=imax[1])
        && (i!=imax[2])
        && (i!=imin[0])
        && (i!=imin[1])
        && (i!=imin[2]) )  {vmin[2]=zt[i]; imin[2]=i;}
  }
  //printf("imax = [%d; %d; %d];\n",imax[0],imax[1],imax[2]);
  //printf("imin = [%d; %d; %d];\n",imin[0],imin[1],imin[2]);
  
  // Computes the volume of the corresponding AABB
  VolumeAABB = (vmax[0]-vmin[0])*(vmax[1]-vmin[1])*(vmax[2]-vmin[2]);

  // based on the 6 bounding points, searches the orthogonal system
  FldArrayF Pts(6,3);
  for (E_Int i=0; i<3;i++)
  {
    Pts(i,1)=xt[imin[i]];
    Pts(i,2)=yt[imin[i]];
    Pts(i,3)=zt[imin[i]];
  }
  for (E_Int i=0; i<3;i++)
  {
    Pts(i+3,1)=xt[imax[i]];
    Pts(i+3,2)=yt[imax[i]];
    Pts(i+3,3)=zt[imax[i]];
  }
  //cout << "Bounding points:" << endl;
  //cout << Pts << endl;
  FldArrayF Edges(5);
  for (E_Int e=0; e<5;e++)
  {
     Edges[e] = sqrt(pow(Pts(e+1,1)-Pts(0,1),2) +
                     pow(Pts(e+1,2)-Pts(0,2),2) +
                     pow(Pts(e+1,3)-Pts(0,3),2));
  }

  // Looks for the 2 smallest edges
  SmallEdges[0] = Edges[0];
  SmallEdges[1] = Edges[0];
  // This is the smallest edge
  for (E_Int e=0; e < 5; e++)
  {
     if (Edges[e] <= SmallEdges[0])  
     {
        SmallEdges[0]=Edges[e];
        SmallEdgesInd[0]=e;
     }
  }
  //v0
  for (E_Int i=0; i<3;i++)
  {
     v0box[i]=Pts(SmallEdgesInd[0]+1,i+1)-Pts(0,i+1);
  }
  for (E_Int i=0; i<3;i++) {v0box[i]=v0box[i]/SmallEdges[0];}
  //cout << "v0box" << endl;
  //cout << v0box << endl;

  // This is the second smallest edge being orthogonal to the previous one
  E_Float v0v1;
  // Temporarily assigns SmallEdges[1] the maximum value of Edges for
  // consistency reasons
  SmallEdges[1] = Edges[0];
  for (E_Int e=0; e < 5; e++)
  {
     if (Edges[e] >= SmallEdges[1])  
     {
        SmallEdges[1]=Edges[e];
     }
  }
  //cout << "e   Edges[e]  SmallEdges[1]  SmallEdges[0]" << endl;
  for (E_Int e=0; e < 5; e++)
  {
  //printf("%d   %0.6f  %0.8f    %0.8f\n",e,Edges[e],SmallEdges[1],SmallEdges[0]);
     if ((Edges[e] <= SmallEdges[1]+Tol) && (e!=SmallEdgesInd[0]))
     {
        //v1box
        for (E_Int i=0; i<3;i++)
        {
           v1box[i]=Pts(e+1,i+1)-Pts(0,i+1);
        }
        for (E_Int i=0; i<3;i++) {v1box[i]=v1box[i]/Edges[e];}

        v0v1= v0box[0]*v1box[0] + v0box[1]*v1box[1] + v0box[2]*v1box[2];
        //printf("%d   v0v1= %g\n",e,v0v1);
        if (K_FUNC::E_abs(v0v1) < Tol)
        {
        SmallEdges[1]=Edges[e];
        SmallEdgesInd[1]=e;
        //cout << "Verfied at e=" << e << endl;
        }     
     }
  }
/*
  cout << "SmallEdges:" << endl;
  cout << SmallEdges << endl;
  cout << "SmallEdgesInd" << endl;
  cout << SmallEdgesInd << endl;
*/
  //v1box
  for (E_Int i=0; i<3;i++)
  {
     v1box[i]=Pts(SmallEdgesInd[1]+1,i+1)-Pts(0,i+1);
  }
  for (E_Int i=0; i<3;i++) {v1box[i]=v1box[i]/SmallEdges[1];}
  //cout << "v1box" << endl;
  //cout << v1box << endl;

  //v2box
  v2box[0]=v0box[1]*v1box[2]-v1box[1]*v0box[2];
  v2box[1]=-v0box[0]*v1box[2]+v1box[0]*v0box[2];
  v2box[2]=v0box[0]*v1box[1]-v1box[0]*v0box[1];
  }

  if (SmallEdgesInd[1]!=-1) // It is a box
  {
    for (E_Int i=0; i<3;i++)
    {
        v0[i]=v0box[i];
        v1[i]=v1box[i];
        v2[i]=v2box[i];
    }
  }


  // Constructs the transposed Eigenvectors matrix
  FldArrayF Eig(3,3);
  for (E_Int j=0; j < 3; j++)
  {
      Eig(0,j+1) = v0[j];
      Eig(1,j+1) = v1[j];
      Eig(2,j+1) = v2[j];
  }
/*
  printf("Eig = [%12g, %12g, %12g;\n %12g, %12g, %12g;\n %12g, %12g, %12g];\n",
                                              Eig(0,1),Eig(0,2),Eig(0,3),
                                              Eig(1,1),Eig(1,2),Eig(1,3),
                                              Eig(2,1),Eig(2,2),Eig(2,3));
*/
  /* Calculates the OBB dimensions */
  // Projects the mesh coordinates into the OBB main axes
  FldArrayF Projv0(nt);
  FldArrayF Projv1(nt);
  FldArrayF Projv2(nt);
  for (E_Int i=0; i < nt; i++)
  {
     Projv0[i] = xt[i]*v0[0] + yt[i]*v0[1] + zt[i]*v0[2];
     Projv1[i] = xt[i]*v1[0] + yt[i]*v1[1] + zt[i]*v1[2];
     Projv2[i] = xt[i]*v2[0] + yt[i]*v2[1] + zt[i]*v2[2];
  }
  // Calculates the OBB radii
  vmin[0]=Projv0[0];
  vmax[0]=Projv0[0];  
  vmin[1]=Projv1[0];
  vmax[1]=Projv1[0]; 
  vmin[2]=Projv2[0];
  vmax[2]=Projv2[0];
  for (E_Int i=0; i < nt; i++)
  {
     if (Projv0[i] >= vmax[0])  {vmax[0]=Projv0[i];}
     if (Projv0[i] <= vmin[0])  {vmin[0]=Projv0[i];}
     if (Projv1[i] >= vmax[1])  {vmax[1]=Projv1[i];}
     if (Projv1[i] <= vmin[1])  {vmin[1]=Projv1[i];}
     if (Projv2[i] >= vmax[2])  {vmax[2]=Projv2[i];}
     if (Projv2[i] <= vmin[2])  {vmin[2]=Projv2[i];}
  }
  
  FldArrayF fobb(8,3);
  E_Int MakeAABB = 0;
  // Computes the volume of the corresponding OBB
  E_Float VolumeOBB;
  VolumeOBB = (vmax[0]-vmin[0])*(vmax[1]-vmin[1])*(vmax[2]-vmin[2]);
  if (isBoxShaped == 1)
  {
    if (VolumeAABB < VolumeOBB)
    {
      MakeAABB = 1;
    }
  }

  if (MakeAABB == 1)
  {
    // Point 0
    fobb(0,1) = xt[imin[0]];
    fobb(0,2) = yt[imin[1]];
    fobb(0,3) = zt[imin[2]];
    // Point 1
    fobb(1,1) = xt[imax[0]];
    fobb(1,2) = yt[imin[1]];
    fobb(1,3) = zt[imin[2]];
    // Point 2
    fobb(2,1) = xt[imin[0]];
    fobb(2,2) = yt[imax[1]];
    fobb(2,3) = zt[imin[2]];
    // Point 3
    fobb(3,1) = xt[imax[0]];
    fobb(3,2) = yt[imax[1]];
    fobb(3,3) = zt[imin[2]];
    // Point 4
    fobb(4,1) = xt[imin[0]];
    fobb(4,2) = yt[imin[1]];
    fobb(4,3) = zt[imax[2]];
    // Point 5
    fobb(5,1) = xt[imax[0]];
    fobb(5,2) = yt[imin[1]];
    fobb(5,3) = zt[imax[2]];        
    // Point 6
    fobb(6,1) = xt[imin[0]];
    fobb(6,2) = yt[imax[1]];
    fobb(6,3) = zt[imax[2]];
    // Point 7
    fobb(7,1) = xt[imax[0]];
    fobb(7,2) = yt[imax[1]];
    fobb(7,3) = zt[imax[2]];
    //cout << "makes AABB" << endl;
  }
  else // makes OBB
  {
  //cout << "makes OBB" << endl;
  // Calculates the rotation matrix,
  // which is the inverse of the EigenvectorT matrix
  FldArrayF Rot(3,3); Rot.setAllValuesAtNull();
  /*
  if ( (K_FUNC::E_abs(Eig(0,1)) < Tol)
       && (K_FUNC::E_abs(Eig(1,2)) < Tol)
       && (K_FUNC::E_abs(Eig(2,3)) < Tol) )
  {
      Eig(0,1) = 1.;
      Eig(1,2) = 1.;
      Eig(2,3) = 1.;
      Rot(0,1) = 1.;
      Rot(1,2) = 1.;
      Rot(2,3) = 1.;      
  }
  */
  K_LINEAR::inv3(Eig.begin(), Rot.begin());

  //cout << "makes OBB" << endl;
  //Constructs the OBB array by going back to the absolute reference frame  
  FldArrayF VectorRel(3);
  FldArrayF VectorAbs(3);
  
  // Point 0
  VectorRel[0] = vmin[0];
  VectorRel[1] = vmin[1];
  VectorRel[2] = vmin[2];
  K_LINEAR::prodv(3, 3, Rot.begin(), VectorRel.begin(), VectorAbs.begin());
  fobb(0,1) = VectorAbs[0];
  fobb(0,2) = VectorAbs[1];
  fobb(0,3) = VectorAbs[2];
  // Point 1
  VectorRel[0] = vmax[0];
  VectorRel[1] = vmin[1];
  VectorRel[2] = vmin[2];
  K_LINEAR::prodv(3, 3, Rot.begin(), VectorRel.begin(), VectorAbs.begin());
  fobb(1,1) = VectorAbs[0];
  fobb(1,2) = VectorAbs[1];
  fobb(1,3) = VectorAbs[2];
  // Point 2
  VectorRel[0] = vmin[0];
  VectorRel[1] = vmax[1];
  VectorRel[2] = vmin[2];
  K_LINEAR::prodv(3, 3, Rot.begin(), VectorRel.begin(), VectorAbs.begin());
  fobb(2,1) = VectorAbs[0];
  fobb(2,2) = VectorAbs[1];
  fobb(2,3) = VectorAbs[2];
  // Point 3
  VectorRel[0] = vmax[0];
  VectorRel[1] = vmax[1];
  VectorRel[2] = vmin[2];
  K_LINEAR::prodv(3, 3, Rot.begin(), VectorRel.begin(), VectorAbs.begin());
  fobb(3,1) = VectorAbs[0];
  fobb(3,2) = VectorAbs[1];
  fobb(3,3) = VectorAbs[2];
  // Point 4
  VectorRel[0] = vmin[0];
  VectorRel[1] = vmin[1];
  VectorRel[2] = vmax[2];
  K_LINEAR::prodv(3, 3, Rot.begin(), VectorRel.begin(), VectorAbs.begin());
  fobb(4,1) = VectorAbs[0];
  fobb(4,2) = VectorAbs[1];
  fobb(4,3) = VectorAbs[2];
  // Point 5
  VectorRel[0] = vmax[0];
  VectorRel[1] = vmin[1];
  VectorRel[2] = vmax[2];
  K_LINEAR::prodv(3, 3, Rot.begin(), VectorRel.begin(), VectorAbs.begin());
  fobb(5,1) = VectorAbs[0];
  fobb(5,2) = VectorAbs[1];
  fobb(5,3) = VectorAbs[2];
  // Point 6
  VectorRel[0] = vmin[0];
  VectorRel[1] = vmax[1];
  VectorRel[2] = vmax[2];
  K_LINEAR::prodv(3, 3, Rot.begin(), VectorRel.begin(), VectorAbs.begin());
  fobb(6,1) = VectorAbs[0];
  fobb(6,2) = VectorAbs[1];
  fobb(6,3) = VectorAbs[2];
  // Point 7
  VectorRel[0] = vmax[0];
  VectorRel[1] = vmax[1];
  VectorRel[2] = vmax[2];
  K_LINEAR::prodv(3, 3, Rot.begin(), VectorRel.begin(), VectorAbs.begin());
  fobb(7,1) = VectorAbs[0];
  fobb(7,2) = VectorAbs[1];
  fobb(7,3) = VectorAbs[2];
  }
/*
  printf("xt = [%12g, %12g, %12g, %12g, %12g, %12g, %12g, %12g];\n",
                xt[0],xt[1],xt[2],xt[3],xt[4],xt[5],xt[6],xt[7]);
  printf("yt = [%12g, %12g, %12g, %12g, %12g, %12g, %12g, %12g];\n",
                yt[0],yt[1],yt[2],yt[3],yt[4],yt[5],yt[6],yt[7]);
  printf("zt = [%12g, %12g, %12g, %12g, %12g, %12g, %12g, %12g];\n",
                zt[0],zt[1],zt[2],zt[3],zt[4],zt[5],zt[6],zt[7]);
*/
  
  // Constructs the Python array
  if (res == 1) 
  {
    PyObject* tpl = K_ARRAY::buildArray3(fobb, "x,y,z", 2, 2, 2, api);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else 
  {
    PyObject* tpl = K_ARRAY::buildArray3(fobb, "x,y,z", 2, 2, 2, api);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
}


