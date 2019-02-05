/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "pdm_poly_surf_gen.h"
#include "pdm_config.h"
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#define ABS(a)     ((a) <  0  ? -(a) : (a))
#define MIN(a,b)   ((a) > (b) ?  (b) : (a))

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


static double random01()
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / ABS(rsigna - rsignb);
  double resultat = sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}

void PDM_poly_surf_gen
(
PDM_MPI_Comm     localComm,
double       xmin,
double       xmax,
double       ymin,
double       ymax,
int          haveRandom,
int          initRandom,
PDM_g_num_t   nx,
PDM_g_num_t   ny,
PDM_g_num_t  *nGFace,
PDM_g_num_t  *nGVtx,
PDM_g_num_t  *nGEdge,
int         *dNVtx,
double     **dVtxCoord,
int         *dNFace,
int        **dFaceVtxIdx,
PDM_g_num_t **dFaceVtx,
PDM_g_num_t **dFaceEdge,
int         *dNEdge,
PDM_g_num_t **dEdgeVtx,
PDM_g_num_t **dEdgeFace,
int         *nEdgeGroup,
int        **dEdgeGroupIdx,
PDM_g_num_t **dEdgeGroup
)
{

  int nRank;
  PDM_MPI_Comm_size(localComm, &nRank);
  
  int localRank;
  PDM_MPI_Comm_rank(localComm, &localRank);
  
  const double coefRand = 0.3;
  srand(initRandom);
  
  /* Comptage des entites du maillage */
  /* -------------------------------- */

  /* nx et ny doivent etre pair et superieur a 4 */
  
  PDM_g_num_t nx1 = nx;
  PDM_g_num_t ny1 = ny;
  
  if (nx1 < 6)
    nx1 = 6;
  
  if (ny1 < 6)
    ny1 = 6;
  
  if (nx1 % 2 != 0)
    nx1 += 1;
  
  if (ny1 % 2 != 0)
    ny1 += 1;

  /* cotes d un octogone */
  
  double cote1;
  if (nx1 == 4)
    cote1 = (xmax - xmin) / (sqrt(2) + 1);
  else
    cote1 = (xmax - xmin) / ((nx1 - 2)/2 + ((nx1 - 2)/2 - 1) * sqrt(2) + sqrt(2));
  double cote2;
  if (ny1 == 4)
    cote2 = (ymax - ymin) / (sqrt(2) + 1);
  else
    cote2 = (ymax - ymin) / ((ny1 - 2)/2 + ((ny1 - 2)/2 - 1) * sqrt(2) + sqrt(2));
  
  /* Comptage des Sommets */
  
  PDM_g_num_t dNVtxTotal = 0;
  
  PDM_g_num_t cpt = 0;
  PDM_g_num_t cptMax;
  
  if (ny1 > 4)
    cptMax = ny1 + (ny1-4)/2;
  else
    cptMax = ny1;
  
  while(cpt < cptMax) {
    if (cpt % 3 == 1 || cpt % 3 == 2) {
      dNVtxTotal +=  nx1/2;
    }
    else if ((cpt % 3) == 0) {
      if ((cpt == (cptMax-1)) || (cpt == 0)) {
        dNVtxTotal++;
      }
      
      dNVtxTotal++;
      
      for (PDM_g_num_t ix = 2; ix < nx1-1; ix++) {
        dNVtxTotal++;
      }
      if ((cpt == (cptMax-1)) || (cpt == 0)) {
        dNVtxTotal++;
      }
    }
    cpt++;
  }
  
  /* Nombre de limites */

  *nEdgeGroup = 4;
  *dEdgeGroupIdx = (int *) malloc(sizeof(int) * (*nEdgeGroup + 1));
  (*dEdgeGroupIdx)[0] = 0;
  
  PDM_g_num_t *dNEdgeLim13Rank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (nRank+1));
  dNEdgeLim13Rank[0] = 0;
  for (int i = 0; i < nRank; i++)
    dNEdgeLim13Rank[i+1] = (PDM_g_num_t) (nx1 - 1) / nRank;

  PDM_g_num_t _reste = (nx1 - 1) % nRank;
  int reste =  (int) _reste;

  for (int i = 0; i < nRank; i++) {
    if (i < reste) {
      dNEdgeLim13Rank[i+1] += 1;
    }
  }

  
  PDM_g_num_t dNEdgeLim13 = dNEdgeLim13Rank[localRank + 1];

  for (int i = 0; i < nRank; i++) {
    dNEdgeLim13Rank[i+1] = dNEdgeLim13Rank[i+1] + dNEdgeLim13Rank[i];
  }

  PDM_g_num_t *dNEdgeLim24Rank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (nRank+1));
  dNEdgeLim24Rank[0] = 0;
  for (int i = 0; i < nRank; i++)
    dNEdgeLim24Rank[i+1] = (PDM_g_num_t) (ny1 - 1) / nRank;

  _reste = (ny1 - 1) % nRank;
  reste =  (int) (_reste);

  for (int i = 0; i < nRank; i++) {
    if (i < reste) {
      dNEdgeLim24Rank[i+1] += 1;
    }
  }
   
  PDM_g_num_t dNEdgeLim24 = dNEdgeLim24Rank[localRank + 1];

  for (int i = 0; i < nRank; i++) {
    dNEdgeLim24Rank[i+1] = dNEdgeLim24Rank[i+1] + dNEdgeLim24Rank[i];
  }

  *dEdgeGroup = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (2*dNEdgeLim13 + 2*dNEdgeLim24));

  (*dEdgeGroupIdx)[1] = (int) dNEdgeLim13;
  (*dEdgeGroupIdx)[2] = (int) dNEdgeLim24;
  (*dEdgeGroupIdx)[3] = (int) dNEdgeLim13;
  (*dEdgeGroupIdx)[4] = (int) dNEdgeLim24;

  for (int i= 0; i < 4; i++)
    (*dEdgeGroupIdx)[i+1] = (*dEdgeGroupIdx)[i] + (*dEdgeGroupIdx)[i+1];


  /* Construction des sommets */
  /* ------------------------ */
  
  PDM_g_num_t *dNVtxRank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (nRank+1));
  
  dNVtxRank[0] = 0;
  for (int i = 0; i < nRank; i++)
    dNVtxRank[i+1] = (PDM_g_num_t) dNVtxTotal / nRank;

  _reste = dNVtxTotal % nRank;
  reste =  (int) (_reste);
  
  for (int i = 0; i < nRank; i++) {
    if (i < reste) {
      dNVtxRank[i+1] += 1;
    }
  }
  
  *dNVtx = (int) dNVtxRank[localRank + 1];
  
  for (int i = 0; i < nRank; i++) {
    dNVtxRank[i+1] = dNVtxRank[i+1] + dNVtxRank[i];
  }

  *dVtxCoord = (double*) malloc (sizeof(double) * 3 * (*dNVtx));

  /* Construction des elements et des aretes */
  /* --------------------------------------- */

  /* Comptage des elements et des aretes */

  PDM_g_num_t  dNFaceTotal = 0;

  /* Triangles */

  /* -- Premiere ligne */
  dNFaceTotal += nx1/2;
  
  /* -- Autres triangles (un a gauche un a droite */
  PDM_g_num_t nbLi = (ny1-4)/2;
  dNFaceTotal += 2*nbLi;

  /* -- Derniere ligne */
  dNFaceTotal += nx1/2;

  /* Quadrangles */
  PDM_g_num_t nxQuad = (nx1-4)/2;
  PDM_g_num_t nyQuad = (ny1-4)/2;
  dNFaceTotal += nyQuad * nxQuad;

  /* Polygones */
  PDM_g_num_t nxPoly = (nx1-2)/2;
  PDM_g_num_t nyPoly = (ny1-2)/2;
  dNFaceTotal += nxPoly * nyPoly;
  
  PDM_g_num_t dNEdgeTotal = (6 * nxPoly + 1) * nyPoly + nxPoly; /* Aretes touchant un polygone */
  dNEdgeTotal += nx1 + ny1; /* Aretes des cotes ne touchant pa sun polygone */
  
  PDM_g_num_t nPoly = nxPoly * nyPoly;
  PDM_g_num_t nQuad = nxQuad * nyQuad;
  PDM_g_num_t nTri  = dNFaceTotal - nQuad - nPoly;
  
  /* Allocation des elts */
  
  PDM_g_num_t *dNFaceRank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (nRank+1));
  PDM_g_num_t *dNEdgeRank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (nRank+1));
  
  dNFaceRank[0] = 0;
  for (int i = 0; i < nRank; i++)
    dNFaceRank[i+1] = dNFaceTotal / nRank ;

  _reste = dNFaceTotal % nRank;
  reste =  (int) _reste ;
  for (int i = 0; i < nRank; i++) {
    if (i < reste) {
      dNFaceRank[i+1] += 1;
    }
  }
  
  dNEdgeRank[0] = 0;
  for (int i = 0; i < nRank; i++)
    dNEdgeRank[i+1] = dNEdgeTotal / nRank ;
  
  _reste = dNEdgeTotal % nRank;
  reste =  (int) _reste;
  for (int i = 0; i < nRank; i++) {
    if (i < reste) {
      dNEdgeRank[i+1] += 1;
    }
  }
  
  *dNFace = (int) dNFaceRank[localRank+1];
  *dNEdge = (int) dNEdgeRank[localRank+1];
  
  for (int i = 0; i < nRank; i++) {
    dNFaceRank[i+1] = dNFaceRank[i+1] + dNFaceRank[i];
    dNEdgeRank[i+1] = dNEdgeRank[i+1] + dNEdgeRank[i];
  }
  
  *dFaceVtxIdx = (int *) malloc(sizeof(int) * ((*dNFace)+1));
  *dFaceVtx    = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 8 * (*dNFace));
  *dFaceEdge      = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 8 * (*dNFace));
  *dEdgeVtx    = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 2 * (*dNEdge));
  *dEdgeFace      = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 2 * (*dNEdge));

  /* Calcul des entites */
  /* ------------------ */
  
  /* Calcul des coordonnees */
  
  PDM_g_num_t dNVtxTmp = 0;

  cpt = 0;
  if (ny1 > 4)
    cptMax = ny1 + (ny1-4)/2;
  else
    cptMax = ny1;
  
  double ycourant = ymin;
  const double eps = 1e-5;
  
  while (cpt < cptMax) {
    if (cpt % 3 == 1 || cpt % 3 == 2) {
      for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
        if ((dNVtxRank[localRank] <= dNVtxTmp) && (dNVtxRank[localRank+1] > dNVtxTmp)) {
          PDM_g_num_t _localIdx = dNVtxTmp - dNVtxRank[localRank];
          int localIdx = (int) _localIdx;
          (*dVtxCoord)[3*localIdx]   = xmin + ix * (1+sqrt(2)) * cote1;
          (*dVtxCoord)[3*localIdx+1] = ycourant;
          (*dVtxCoord)[3*localIdx+2] = 0.;
        }
        dNVtxTmp++;
      }
    }
    else if ((cpt % 3) == 0) {
      if ((cpt == (cptMax-1)) || (cpt == 0)) {
        if ((dNVtxRank[localRank] <= dNVtxTmp) && (dNVtxRank[localRank+1] > dNVtxTmp)) {
          PDM_g_num_t _localIdx = dNVtxTmp - dNVtxRank[localRank];
          int localIdx = (int) _localIdx;
          (*dVtxCoord)[3*localIdx]   = xmin;
          (*dVtxCoord)[3*localIdx+1] = ycourant;
          (*dVtxCoord)[3*localIdx+2] = 0.;
        }
        dNVtxTmp++;
      }
      
      if ((dNVtxRank[localRank] <= dNVtxTmp) && (dNVtxRank[localRank+1] > dNVtxTmp)) {
        PDM_g_num_t _localIdx = dNVtxTmp - dNVtxRank[localRank];
        int localIdx = (int) _localIdx;
        (*dVtxCoord)[3*localIdx]   = xmin + cote1/sqrt(2);
        (*dVtxCoord)[3*localIdx+1] = ycourant;
        (*dVtxCoord)[3*localIdx+2] = 0.;
      }
      dNVtxTmp++;
      
      double xbase = xmin + cote1/sqrt(2);

      PDM_g_num_t nx11 = (nx1-2)/2;
      for (PDM_g_num_t ix = 0; ix < nx11; ix++) {
        if ((dNVtxRank[localRank] <= dNVtxTmp) && (dNVtxRank[localRank+1] > dNVtxTmp)) {
          PDM_g_num_t _localIdx = dNVtxTmp - dNVtxRank[localRank];
          int localIdx = (int) _localIdx;
          (*dVtxCoord)[3*localIdx] = xbase + (ix+1) * cote1 + ix * cote1*sqrt(2);
          (*dVtxCoord)[3*localIdx+1] = ycourant;
          (*dVtxCoord)[3*localIdx+2] = 0.;
        }
        dNVtxTmp++;
        if (ix < (nx11 - 1)) {
          if ((dNVtxRank[localRank] <= dNVtxTmp) && (dNVtxRank[localRank+1] > dNVtxTmp)) {
            PDM_g_num_t _localIdx = dNVtxTmp - dNVtxRank[localRank];
            int localIdx = (int) _localIdx;
            (*dVtxCoord)[3*localIdx] = xbase + (ix+1) * cote1 + (ix+1) * cote1*sqrt(2);
            (*dVtxCoord)[3*localIdx+1] = ycourant;
            (*dVtxCoord)[3*localIdx+2] = 0.;
          }
          dNVtxTmp++;
        }
      }
      
      if ((cpt == (cptMax-1)) || (cpt == 0)) {
        if ((dNVtxRank[localRank] <= dNVtxTmp) && (dNVtxRank[localRank+1] > dNVtxTmp)) {
          PDM_g_num_t _localIdx = dNVtxTmp - dNVtxRank[localRank];
          int localIdx = (int) _localIdx;
          (*dVtxCoord)[3*localIdx]   = xmax;
          (*dVtxCoord)[3*localIdx+1] = ycourant;
          (*dVtxCoord)[3*localIdx+2] = 0.;
        }
        dNVtxTmp++;
      }
    }
    cpt++;
    if ((cpt % 3 == 1) || (cpt % 3 == 0))
      ycourant += cote2/sqrt(2);
    else
      ycourant += cote2;
  }

  *nGVtx = dNVtxTmp;

  /* Perturbation des coordonnees + Creation d'une courbure en Z */
  
  for (PDM_g_num_t ix = 0; ix <(*dNVtx) ; ix++) {
    if (ABS(xmin-(*dVtxCoord)[3*ix]) > eps &&
        ABS(xmax-(*dVtxCoord)[3*ix]) > eps &&
        ABS(ymax-(*dVtxCoord)[3*ix+1]) > eps &&
        ABS(ymin-(*dVtxCoord)[3*ix+1]) > eps) {
      if (haveRandom != 0) {
        (*dVtxCoord)[3*ix]   += random01() * coefRand * cote1;
        (*dVtxCoord)[3*ix+1] += random01() * coefRand * cote2;
      }
    }
    //dVtxCoord[3*ix+2]=2*sin(3*dVtxCoord[3*ix+1])*sin(3*dVtxCoord[3*ix]);
    //dVtxCoord[3*ix+2]=20*sin(dVtxCoord[3*ix+1]/5.)*sin(dVtxCoord[3*ix]/5.);
    (*dVtxCoord)[3*ix+2] = 0.;
  }
  
  /* Construction simultanée des connectivites des elements et des aretes internes */
  /* ----------------------------------------------------------------------------- */

  int dNFaceTmp = 0;
  int dNEdgeTmp = 0;
  int dNFaceAbs = 0;
  int dNEdgeAbs = 0;
  
  PDM_g_num_t n1;
  PDM_g_num_t n2;
  PDM_g_num_t n3;
  PDM_g_num_t n4;
  PDM_g_num_t n5;
  PDM_g_num_t n6;
  PDM_g_num_t n7;
  PDM_g_num_t n8;
  
  /* Triangle */
    
  /* -- Premiere ligne */
  
  n1 = 1;
  n2 = n1+nx1;
  (*dFaceVtxIdx)[0] = 0;
  
  for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
    PDM_g_num_t ix1 = 2 * ix;

    if ((dNFaceRank[localRank] <= dNFaceAbs) && (dNFaceRank[localRank+1] > dNFaceAbs)) {
      int ideb = (*dFaceVtxIdx)[dNFaceTmp] ;
      (*dFaceVtx)[ideb]   = n1+ix1;
      (*dFaceVtx)[ideb+1] = n1+ix1+1;
      (*dFaceVtx)[ideb+2] = n2+ix;
      if (ix == 0)
        (*dFaceEdge)[ideb]   = 1;
      else
        (*dFaceEdge)[ideb]   = nTri + 4 + 6*(ix-1) + 2;
      if (ix == 0)
        (*dFaceEdge)[ideb+1]   = dNEdgeAbs + 2;
      else
        (*dFaceEdge)[ideb+1]   = dNEdgeAbs + 1;
      if (ix == (nx1/2 - 1))
        (*dFaceEdge)[ideb+2]   = dNEdgeAbs + 2;
      else if (ix == ((nx1/2 - 1) - 1))
        (*dFaceEdge)[ideb+2]   = nTri + 4 + 6*ix + 7;
      else
        (*dFaceEdge)[ideb+2]   = nTri + 4 + 6*ix + 6;
      dNFaceTmp += 1;
      (*dFaceVtxIdx)[dNFaceTmp] = ideb + 3;
    }
    dNFaceAbs += 1;

    if (ix == 0) {
      if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
        (*dEdgeVtx)[2*dNEdgeTmp]     = n2+ix;
        (*dEdgeVtx)[2*dNEdgeTmp + 1] = n1+ix1;
        (*dEdgeFace)[2*dNEdgeTmp]       = dNFaceAbs;
        (*dEdgeFace)[2*dNEdgeTmp + 1]   = 0;
        dNEdgeTmp += 1;
      }
      dNEdgeAbs += 1;
    }

    if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
      (*dEdgeVtx)[2*dNEdgeTmp]     = n1+ix1;
      (*dEdgeVtx)[2*dNEdgeTmp + 1] = n1+ix1+1;
      (*dEdgeFace)[2*dNEdgeTmp]       = dNFaceAbs;
      (*dEdgeFace)[2*dNEdgeTmp + 1]   = 0;
      dNEdgeTmp += 1;
    }
    dNEdgeAbs += 1;

    if (ix == (nx1/2 - 1)) {
      if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
        (*dEdgeVtx)[2*dNEdgeTmp]     = n1+ix1+1;
        (*dEdgeVtx)[2*dNEdgeTmp + 1] = n2+ix;
        (*dEdgeFace)[2*dNEdgeTmp]       = dNFaceAbs;
        (*dEdgeFace)[2*dNEdgeTmp + 1]   = 0;
        dNEdgeTmp += 1;
      }
      dNEdgeAbs += 1;
    }

  }
  
  /* -- Autres triangles (un a gauche un a droite */
  
  n1 = 1 + nx1 + nx1/2;
  n2 = n1 + nx1/2;
  n3 = n2 + nx1-2;
  
  for (PDM_g_num_t itri = 0; itri < nbLi; itri++) {
    n4 = n1 + nx1/2 - 1;
    n5 = n2 + nx1-2 - 1;
    n6 = n3 + nx1/2 - 1;

    if ((dNFaceRank[localRank] <= dNFaceAbs) && (dNFaceRank[localRank+1] > dNFaceAbs)) {
      int ideb = (*dFaceVtxIdx)[dNFaceTmp];
      (*dFaceVtx)[ideb]   = n1;
      (*dFaceVtx)[ideb+1] = n2;
      (*dFaceVtx)[ideb+2] = n3;
      (*dFaceEdge)[ideb]     = dNEdgeAbs + 1;
      (*dFaceEdge)[ideb+1]   = nTri + 4 + (6 * nxPoly + 1) * itri + 4;
      if (itri == nbLi - 1) 
        (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (itri+1) + 7;
      else
        (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (itri+1) + 6;
      dNFaceTmp += 1;
      (*dFaceVtxIdx)[dNFaceTmp] = ideb + 3;
    }
    dNFaceAbs += 1;

    if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
      (*dEdgeVtx)[2*dNEdgeTmp]     = n3;
      (*dEdgeVtx)[2*dNEdgeTmp + 1] = n1;
      (*dEdgeFace)[2*dNEdgeTmp]       = dNFaceAbs;
      (*dEdgeFace)[2*dNEdgeTmp + 1]   = 0;
      dNEdgeTmp += 1;
    }
    dNEdgeAbs += 1;

    if ((dNFaceRank[localRank] <= dNFaceAbs) && (dNFaceRank[localRank+1] > dNFaceAbs)) {
      int ideb = (*dFaceVtxIdx)[dNFaceTmp];
      (*dFaceVtx)[ideb]   = n4;
      (*dFaceVtx)[ideb+1] = n6;
      (*dFaceVtx)[ideb+2] = n5;
      (*dFaceEdge)[ideb]     = dNEdgeAbs + 1;
      if (itri == nbLi - 1) 
        (*dFaceEdge)[ideb+1]   = nTri + 4 + (6 * nxPoly + 1) * (itri + 1) + (7 * nxPoly + 1) - 6;
      else
        (*dFaceEdge)[ideb+1]   = nTri + 4 + (6 * nxPoly + 1) * (itri+2) + 2;
      (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (itri + 1) - 3;
      dNFaceTmp += 1;
      (*dFaceVtxIdx)[dNFaceTmp] = ideb + 3;
    }
    dNFaceAbs += 1;

    if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
      (*dEdgeVtx)[2*dNEdgeTmp]     = n4;
      (*dEdgeVtx)[2*dNEdgeTmp + 1] = n6;
      (*dEdgeFace)[2*dNEdgeTmp]       = dNFaceAbs;
      (*dEdgeFace)[2*dNEdgeTmp + 1]   = 0;
      dNEdgeTmp += 1;
    }
    dNEdgeAbs += 1;

    n1 = n3 + nx1/2;
    n2 = n1 + nx1/2;
    n3 = n2 + nx1-2;
  }
  
  /* -- Derniere ligne */
  
  n2 = n1 + nx1/2;
  for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
    PDM_g_num_t ix1 = 2 * ix;
    if ((dNFaceRank[localRank] <= dNFaceAbs) && (dNFaceRank[localRank+1] > dNFaceAbs)) {
      int ideb = (*dFaceVtxIdx)[dNFaceTmp] ;
      (*dFaceVtx)[ideb]   = n1 + ix;
      (*dFaceVtx)[ideb+1] = n2 + ix1 + 1;
      (*dFaceVtx)[ideb+2] = n2 + ix1;
      if (ix == 0)
        (*dFaceEdge)[ideb+1]   = dNEdgeAbs + 2;
      else
        (*dFaceEdge)[ideb+1]   = dNEdgeAbs + 1;
      if (ix == 0)
        (*dFaceEdge)[ideb]   = dNEdgeAbs + 1;
      else if (ix == (nx1/2 - 1))
        (*dFaceEdge)[ideb]   = nTri + 4 + (6 * nxPoly + 1) * (nyPoly - 1) + 6*(ix-1) + 5;
      else
        (*dFaceEdge)[ideb]   = nTri + 4 + (6 * nxPoly + 1) * (nyPoly - 1) + 6*(ix-1) + 3;
      if (ix == (nx1/2 - 1))
        (*dFaceEdge)[ideb+2]   = dNEdgeAbs + 2;
      else if (ix == ((nx1/2 - 1) - 1))
        (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (nyPoly - 1) + 6*(nxPoly - 1) + 7;
      else
        (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (nyPoly - 1) + 6*ix + 5;
      dNFaceTmp += 1;
      (*dFaceVtxIdx)[dNFaceTmp] = ideb + 3;
    }
    dNFaceAbs += 1;

    if (ix == 0) {
      if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
        (*dEdgeVtx)[2*dNEdgeTmp]     = n2 + ix1;
        (*dEdgeVtx)[2*dNEdgeTmp + 1] = n1 + ix;
        (*dEdgeFace)[2*dNEdgeTmp]       = dNFaceAbs;
        (*dEdgeFace)[2*dNEdgeTmp + 1]   = 0;
        dNEdgeTmp += 1;
      }
      dNEdgeAbs += 1;
    }

    if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
      (*dEdgeVtx)[2*dNEdgeTmp]     = n2 + ix1;
      (*dEdgeVtx)[2*dNEdgeTmp + 1] = n2 + ix1 + 1;
      (*dEdgeFace)[2*dNEdgeTmp]       = dNFaceAbs;
      (*dEdgeFace)[2*dNEdgeTmp + 1]   = 0;
      dNEdgeTmp += 1;
    }
    dNEdgeAbs += 1;

    if (ix == (nx1/2 - 1)) {
      if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
        (*dEdgeVtx)[2*dNEdgeTmp]     = n1 + ix;
        (*dEdgeVtx)[2*dNEdgeTmp + 1] = n2 + ix1 + 1;
        (*dEdgeFace)[2*dNEdgeTmp]       = dNFaceAbs;
        (*dEdgeFace)[2*dNEdgeTmp + 1]   = 0;
        dNEdgeTmp += 1;
      }
      dNEdgeAbs += 1;
    }
  }
  
  /* Quadrangle */
  
  for (PDM_g_num_t iy = 0; iy < nyQuad; iy++) {
    for (PDM_g_num_t ix = 0; ix < nxQuad; ix++) {
      n1 = iy*(2*nx1-2) + nx1 + nx1/2 + 1 + ix + 1;
      n2 = iy*(2*nx1-2) + 2*nx1 + 1 + 2*ix + 1;
      n3 = n2 + 1;
      n4 = iy*(2*nx1-2) + 3*nx1 - 2 + 1 + ix + 1;
      if ((dNFaceRank[localRank] <= dNFaceAbs) && (dNFaceRank[localRank+1] > dNFaceAbs)) {
        int ideb = (*dFaceVtxIdx)[dNFaceTmp];
        (*dFaceVtx)[ideb]   = n1;
        (*dFaceVtx)[ideb+1] = n3;
        (*dFaceVtx)[ideb+2] = n4;
        (*dFaceVtx)[ideb+3] = n2;
        (*dFaceEdge)[ideb]     = nTri + 4 + (6 * nxPoly + 1) * iy     + 6*ix     + 3;
        if (ix == nxQuad - 1)
          (*dFaceEdge)[ideb+1]   = nTri + 4 + (6 * nxPoly + 1) * iy     + 6*(ix+1) + 5;
        else
          (*dFaceEdge)[ideb+1]   = nTri + 4 + (6 * nxPoly + 1) * iy     + 6*(ix+1) + 4;
        if (ix == nyQuad - 1)
          if (ix == nxQuad - 1)
            (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (iy+1) + 7*(ix+1) + 8;
          else
            (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (iy+1) + 7*(ix+1) + 7;
        else
          if (ix == nxQuad - 1)
            (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (iy+1) + 6*(ix+1) + 7;
          else
            (*dFaceEdge)[ideb+2]   = nTri + 4 + (6 * nxPoly + 1) * (iy+1) + 6*(ix+1) + 6;
        (*dFaceEdge)[ideb+3]   = nTri + 4 + (6 * nxPoly + 1) * (iy+1) + 6*ix     + 2;
        dNFaceTmp += 1;
        (*dFaceVtxIdx)[dNFaceTmp] = ideb + 4;
      }
      dNFaceAbs += 1;
    }
  }
  
  /* Polygones */
  
  PDM_g_num_t delta = 0;
  PDM_g_num_t ipoLy = 0;
  
  for (PDM_g_num_t iy = 0; iy < nyPoly; iy++) {
    if (iy == 1)
      delta += 2*nx1;
    else if (iy != 0)
      delta += 2*nx1-2;
    for (PDM_g_num_t ix = 0; ix < nxPoly; ix++) {
      ipoLy += 1;
      if (iy == 0)
        n1 = delta + 1 + 2*ix +1;
      else
        n1 = delta + 1 + 2*ix ;
      n2 = n1 + 1;
      n3 = iy*(2*nx1-2) + 1 + nx1 + ix;
      n4 = n3 + 1;
      n5 = iy*(2*nx1-2) + 1 + nx1 + nx1/2 + ix;
      n6 = n5 + 1;
      n7 = iy*(2*nx1-2) + 1 + 2*nx1 + 2*ix;
      if (iy == (nyPoly - 1))
        n7 = iy*(2*nx1-2) + 1 + 2*nx1 + 2*ix + 1;
      n8 = n7 + 1;

      PDM_g_num_t connecPoly[8];
      connecPoly[0]      = n1;
      connecPoly[1]      = n2;
      connecPoly[2]      = n4;
      connecPoly[3]      = n6;
      connecPoly[4]      = n8;
      connecPoly[5]      = n7;
      connecPoly[6]      = n5;
      connecPoly[7]      = n3;

      if ((dNFaceRank[localRank] <= dNFaceAbs) && (dNFaceRank[localRank+1] > dNFaceAbs)) {
        int ideb = (*dFaceVtxIdx)[dNFaceTmp] ;
        for (int k = 0; k < 8; k++)
          (*dFaceVtx)[ideb+k] = connecPoly[k];

        int id = 0;
        (*dFaceEdge)[ideb]     = dNEdgeAbs + (++id);
        (*dFaceEdge)[ideb+1]   = dNEdgeAbs + (++id);
        if (ix == (nxPoly - 1))
          (*dFaceEdge)[ideb+2]   = dNEdgeAbs + (++id);
        else 
          (*dFaceEdge)[ideb+2]   = dNEdgeAbs + (6 * nxPoly + 1) + 5;
        (*dFaceEdge)[ideb+3]   = dNEdgeAbs + (++id);
        if (iy == (nyPoly - 1))
          (*dFaceEdge)[ideb+4]   = dNEdgeAbs + (++id);
        else
          (*dFaceEdge)[ideb+4]   = dNEdgeAbs + (6 * nxPoly + 1) + 1;
        (*dFaceEdge)[ideb+5]   = dNEdgeAbs + (++id);
        (*dFaceEdge)[ideb+6]   = dNEdgeAbs + (++id);
        (*dFaceEdge)[ideb+7]   = dNEdgeAbs + (++id);
        dNFaceTmp += 1;
        (*dFaceVtxIdx)[dNFaceTmp] = ideb + 8;
      }
      dNFaceAbs += 1;
      
      /* Definition de toutes les aretes internes */
      
      for (int k = 0; k < 8; k++) {
        if (!((k == 2) && (ix != (nxPoly - 1))) &&
            !((k == 4) && (iy != (nyPoly - 1)))) {


          if ((dNEdgeRank[localRank] <= dNEdgeAbs) && (dNEdgeRank[localRank+1] > dNEdgeAbs)) {
            (*dEdgeVtx)[2*dNEdgeTmp]     = connecPoly[k];
            (*dEdgeVtx)[2*dNEdgeTmp + 1] = connecPoly[(k + 1)%8];
            (*dEdgeFace)[2*dNEdgeTmp]       = nTri + nQuad + ipoLy;
            if (k == 0) {
              if (iy == 0)
                (*dEdgeFace)[2*dNEdgeTmp + 1] = 0;
              else
                (*dEdgeFace)[2*dNEdgeTmp + 1] = nTri + nQuad + ipoLy - nxPoly;
            }
            else if (k == 1) {
              if (iy == 0)
                (*dEdgeFace)[2*dNEdgeTmp + 1] = ix + 2;
              else {
                if (ix == (nxPoly - 1))
                  (*dEdgeFace)[2*dNEdgeTmp + 1] = nx1/2 + 2*(iy-1) + 2;
                else
                  (*dEdgeFace)[2*dNEdgeTmp + 1] = nTri + (iy - 1)*nxQuad + ix + 1;
              }
            }
            else if (k == 2) {
              (*dEdgeFace)[2*dNEdgeTmp + 1] = 0;
            }
            else if (k == 3) {
              if (iy == (nyPoly - 1))
                (*dEdgeFace)[2*dNEdgeTmp + 1] = nx1/2 + 2*nbLi + ix +2;
              else {
                if (ix == (nxPoly - 1))
                  (*dEdgeFace)[2*dNEdgeTmp + 1] = nx1/2 + 2*(iy+1);
                else
                  (*dEdgeFace)[2*dNEdgeTmp + 1] = nTri + iy*nxQuad + ix + 1;
              }
            }
            else if (k == 4) {
              (*dEdgeFace)[2*dNEdgeTmp + 1] = 0;
            }
            else if (k == 5) {
              if (iy == (nyPoly - 1))
                (*dEdgeFace)[2*dNEdgeTmp + 1] = nx1/2 + 2*nbLi + ix+1;
              else {
                if (ix == 0)
                  (*dEdgeFace)[2*dNEdgeTmp + 1] = nx1/2 + 2*(iy+1) - 1;
                else
                  (*dEdgeFace)[2*dNEdgeTmp + 1] = nTri + iy*nxQuad + ix + 1 - 1;
              }
            }
            else if (k == 6) {
              if (ix == 0)
                (*dEdgeFace)[2*dNEdgeTmp + 1] = 0;
              else
                (*dEdgeFace)[2*dNEdgeTmp + 1] = nTri + nQuad + ipoLy - 1;
            }
            else if (k == 7) {
              if (iy == 0)
                (*dEdgeFace)[2*dNEdgeTmp + 1] = ix + 1;
              else {
                if (ix == 0)
                  (*dEdgeFace)[2*dNEdgeTmp + 1] = nx1/2 + 2*(iy-1) + 1;
                else
                  (*dEdgeFace)[2*dNEdgeTmp + 1] = nTri + (iy-1)*nxQuad + ix + 1 - 1;
              }
            }
            dNEdgeTmp += 1;
          }
          dNEdgeAbs += 1;
        }            
      }
    }
  }

  *nGEdge = dNEdgeAbs;
  *nGFace = dNFaceAbs;

  /* Definition des limites */
  /* ---------------------- */

  int idxDEdgeGroup = 0;

  /* - lim 1 - (triangles ligne bas + octogones bas) */

  PDM_g_num_t dNEdgeGroupAbs = 0;
  for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
    if ((dNEdgeLim13Rank[localRank] <= dNEdgeGroupAbs) && (dNEdgeLim13Rank[localRank+1] > dNEdgeGroupAbs)) {
      (*dEdgeGroup)[idxDEdgeGroup++]  = 2 + ix; 
    }
    ++dNEdgeGroupAbs;
  }

  for (PDM_g_num_t ix = 0; ix < nxPoly; ix++) {
    if ((dNEdgeLim13Rank[localRank] <= dNEdgeGroupAbs) && (dNEdgeLim13Rank[localRank+1] > dNEdgeGroupAbs)) {
      (*dEdgeGroup)[idxDEdgeGroup++]  = nx1 + ny1 + (6 * ix) + 1; 
    }
    ++dNEdgeGroupAbs;
  }

  /* - lim 2 - */

  dNEdgeGroupAbs = 0;
  for (PDM_g_num_t iy = 0; iy < ny1/2; iy++) {
    if ((dNEdgeLim24Rank[localRank] <= dNEdgeGroupAbs) && (dNEdgeLim24Rank[localRank+1] > dNEdgeGroupAbs)) {
      if (iy == (ny1/2 - 1))
        (*dEdgeGroup)[idxDEdgeGroup++]  = nx1 + ny1;
      else
        (*dEdgeGroup)[idxDEdgeGroup++]  = nx1/2 + 2*iy + 2; 
    }
    ++dNEdgeGroupAbs;
  }

  for (PDM_g_num_t iy = 0; iy < nyPoly; iy++) {
    if ((dNEdgeLim24Rank[localRank] <= dNEdgeGroupAbs) && (dNEdgeLim24Rank[localRank+1] > dNEdgeGroupAbs)) {
      if (iy == nyPoly - 1)
        (*dEdgeGroup)[idxDEdgeGroup++]  =  nx1 + ny1 + ((6 * nxPoly) + 1) * (nyPoly - 1) + 7 * (nxPoly - 1) + 3;
      else
        (*dEdgeGroup)[idxDEdgeGroup++]  =  nx1 + ny1 + ((6 * nxPoly) + 1) * (iy+1) - 4;
    }
    ++dNEdgeGroupAbs;
  }

  /* - lim 3 - */

  dNEdgeGroupAbs = 0;
  for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
    if ((dNEdgeLim13Rank[localRank] <= dNEdgeGroupAbs) && (dNEdgeLim13Rank[localRank+1] > dNEdgeGroupAbs)) {
      (*dEdgeGroup)[idxDEdgeGroup++]  = nx1 + ny1 - nx1/2 + ix; 
    }
    ++dNEdgeGroupAbs;
  }

  for (PDM_g_num_t ix = 0; ix < nxPoly; ix++) {
    if ((dNEdgeLim13Rank[localRank] <= dNEdgeGroupAbs) && (dNEdgeLim13Rank[localRank+1] > dNEdgeGroupAbs)) {
      if (ix == (nxPoly - 1))
        (*dEdgeGroup)[idxDEdgeGroup++]  = nx1 + ny1 + ((6 * nxPoly) + 1) * (nyPoly - 1) + 7 * ix + 5;
      else
        (*dEdgeGroup)[idxDEdgeGroup++]  = nx1 + ny1 + ((6 * nxPoly) + 1) * (nyPoly - 1) + 7 * ix + 4;
    }
    ++dNEdgeGroupAbs;
  }

  /* - lim 4 - */

  dNEdgeGroupAbs = 0;
  for (PDM_g_num_t iy = 0; iy < ny1/2; iy++) {
    if ((dNEdgeLim24Rank[localRank] <= dNEdgeGroupAbs) && (dNEdgeLim24Rank[localRank+1] > dNEdgeGroupAbs)) {
      if (iy == 0)
        (*dEdgeGroup)[idxDEdgeGroup++]  = 1; 
      else
        (*dEdgeGroup)[idxDEdgeGroup++]  = 2 + nx1/2 + (2*(iy-1)) + 1;
    }
    ++dNEdgeGroupAbs;
  }

  for (PDM_g_num_t iy = 0; iy < nyPoly; iy++) {
    if ((dNEdgeLim24Rank[localRank] <= dNEdgeGroupAbs) && (dNEdgeLim24Rank[localRank+1] > dNEdgeGroupAbs)) {
      if (iy == nyPoly - 1)
        (*dEdgeGroup)[idxDEdgeGroup++]  =  nx1 + ny1 + ((6 * nxPoly) + 1) * iy + 6;
      else
        (*dEdgeGroup)[idxDEdgeGroup++]  =  nx1 + ny1 + ((6 * nxPoly) + 1) * iy + 5;
    }
    ++dNEdgeGroupAbs;
  }

  free(dNEdgeLim13Rank); 
  free(dNEdgeLim24Rank); 
  free(dNVtxRank);
  free(dNFaceRank);
  free(dNEdgeRank);
}


/* void PROCF(creemaillagepolygone2d_f, CREEMAILLAGEPOLYGONE2D_F) */
/* ( */
/* PDM_MPI_Fint   *localFComm, */
/* double     *xmin, */
/* double     *xmax, */
/* double     *ymin, */
/* double     *ymax, */
/* int        *initRandom, */
/* PDM_g_num_t *nx, */
/* PDM_g_num_t *ny, */
/* int        *dNVtx, */
/* double     *dVtxCoord_f, */
/* int        *dNFace, */
/* int        *dFaceVtxIdx_f, */
/* PDM_g_num_t *dFaceVtx_f, */
/* PDM_g_num_t *dFaceEdge_f,    */
/* int        *dNEdge, */
/* PDM_g_num_t *dEdgeVtx_f, */
/* PDM_g_num_t *dEdgeFace_f */
/* ) */
/* { */
/*   PDM_MPI_Comm localComm = PDM_MPI_Comm_f2c(*localFComm); */
/*   int dNVtx_f = *dNVtx; */
/*   int dNFace_f = *dNFace; */

/*   double *dVtxCoord = NULL; */
/*   int    *dFaceVtxIdx = NULL; */
/*   int    *dFaceVtx = NULL; */

/*   creeMaillagePolygone2D(*order, */
/*                          localComm, */
/*                          *xmin, */
/*                          *xmax, */
/*                          *ymin, */
/*                          *ymax, */
/*                          *initRandom, */
/*                          *nx, */
/*                          *ny, */
/*                          dNVtx, */
/*                          &dVtxCoord, */
/*                          dNFace, */
/*                          &dFaceVtxIdx, */
/*                          &dFaceVtx); */

/*   if (dNVtx_f < *dNVtx) { */
/*     PDM_printf("Augmenter le nombre de sommets Fortran a : %i \n", *dNVtx); */
/*     exit(1); */
/*   } */

/*   if (dNFace_f < *dNFace) { */
/*     PDM_printf("Augmenter le nombre d'elements a : %i \n", *dNFace); */
/*     exit(1); */
/*   } */

/*   if (*lEltsConnecPointer_f < dFaceVtxIdx[*dNFace]) { */
/*     PDM_printf("Augmenter la taille du tableau de connectivite a : %i \n", dFaceVtxIdx[*dNFace]); */
/*     exit(1); */
/*   } */

/*   for(int i = 0; i < 3*(*dNVtx); i++) */
/*     dVtxCoord_f[i] = dVtxCoord[i]; */

/*   for(int i = 0; i < *dNFace + 1; i++) */
/*     dFaceVtxIdx_f[i] = dFaceVtxIdx[i]; */

/*   for(int i = 0; i < dFaceVtxIdx[*dNFace]; i++) */
/*     dFaceVtx_f[i] = dFaceVtx[i]; */

/*   free(dVtxCoord); */
/*   free(dFaceVtxIdx); */
/*   free(dFaceVtx); */
/* } */

#ifdef __cplusplus
}
#endif /* __cplusplus */
