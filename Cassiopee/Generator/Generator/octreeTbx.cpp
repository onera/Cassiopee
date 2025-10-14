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

// Outils pour l'octree

# include "generator.h"
# include "Nuga/include/BbTree.h"
# include "Nuga/include/ArrayAccessor.h"

//=============================================================================
/* Recuperation des elements voisins par arete, avec prise en compte de la non
   conformite des elements 
   IN: npts: nb de points dans l'octree
   IN: xt, yt, zt: coordonnees des points du maillage
   IN: cEV: connectivite element/vertex
   OUT: nghbrs: connectivite elements/elements voisins non conformes
   IN: type: 0: voisin ayant une facette en commun,
             1: voisin ayant une facette ou un point en commun.
   IN: mindh: pas de la plus petite cellule (uniquement si type=1) 
   Retourne -1 si pas possible
*/
//=============================================================================
E_Int K_GENERATOR::getNeighbourElts(
  E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
  K_FLD::FldArrayI& cEV, 
  std::vector< std::vector<E_Int> >& nghbrs, E_Int type, E_Float mindh)
{
  E_Int nvert = cEV.getNfld();
  if (nvert != 4 && nvert != 8) return -1;

  if (type == 0)
  {
    if (nvert == 4) return getNeighbourQuads(npts, xt, yt, zt, cEV, nghbrs);
    else return getNeighbourHexas(npts, xt, yt, zt, cEV, nghbrs);
  }
  else
  {
    if (nvert == 4) return getNeighbourQuads2(npts, xt, yt, zt, cEV, nghbrs,
                                              mindh);
    else return getNeighbourHexas2(npts, xt, yt, zt, cEV, nghbrs, mindh);
  }
}

//=============================================================================
// Retourne les voisins par elt ayant une facette en commun
//=============================================================================
E_Int K_GENERATOR::getNeighbourQuads(
  E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
  K_FLD::FldArrayI& cEV, 
  std::vector< std::vector<E_Int> >& nghbrs)
{
  E_Float eps = 1.e-10;
  E_Int nelts = cEV.getSize();
  E_Int* cn1 = cEV.begin(1); E_Int* cn2 = cEV.begin(2);
  E_Int* cn3 = cEV.begin(3); E_Int* cn4 = cEV.begin(4);

  /* 1 - creation du bbtree */
  typedef K_SEARCH::BoundingBox<3> BBox3DType; 
  std::vector<BBox3DType*> boxes(nelts);// a detruire a la fin
  K_FLD::FldArrayF bbox(nelts,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cEV, xt, yt, zt, bbox);
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);

#pragma omp parallel default(shared)
  {
     E_Float minB0[3];  E_Float maxB0[3];
#pragma omp for
     for (E_Int et = 0; et < nelts; et++)
     {
       minB0[0] = xminp[et]-eps; minB0[1] = yminp[et]-eps; minB0[2] = zminp[et]-eps;
       maxB0[0] = xmaxp[et]+eps; maxB0[1] = ymaxp[et]+eps; maxB0[2] = zmaxp[et]+eps; 
       boxes[et] = new BBox3DType(minB0, maxB0);
     }
  }
  K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);

#pragma omp parallel default(shared)
  {
    E_Int indA1, indB1, indC1, indD1, indA2, indB2, indC2, indD2, et2;
  E_Float xA1, xB1, xC1, xD1, yA1, yB1, yC1, yD1;
  E_Float xA2, xB2, xC2, xD2, yA2, yB2, yC2, yD2;
  E_Float minB[3];  E_Float maxB[3];
#pragma omp for
  for (E_Int et1 = 0; et1 < nelts; et1++)
  { 
    indA1 = cn1[et1]-1; indB1 = cn2[et1]-1; indC1 = cn3[et1]-1; indD1 = cn4[et1]-1;
    xA1 = xt[indA1]; xB1 = xt[indB1]; xC1 = xt[indC1]; xD1 = xt[indD1];
    yA1 = yt[indA1]; yB1 = yt[indB1]; yC1 = yt[indC1]; yD1 = yt[indD1];
    //zA1 = zt[indA1]; zB1 = zt[indB1]; zC1 = zt[indC1]; zD1 = zt[indD1];
    std::vector<E_Int> voisins;
 
    // recherche des candidats 
    minB[0] = xminp[et1]-eps; minB[1] = yminp[et1]-eps; minB[2] = zminp[et1]-eps;
    maxB[0] = xmaxp[et1]+eps; maxB[1] = ymaxp[et1]+eps; maxB[2] = zmaxp[et1]+eps; 
    std::vector<E_Int> indicesBB;
    bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
    for (size_t noi = 0; noi < indicesBB.size(); noi++)
    {
      et2 = indicesBB[noi];
      if (et2 == et1) goto nextet2; 

      indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1; indD2 = cn4[et2]-1;
      xA2 = xt[indA2]; xB2 = xt[indB2]; xC2 = xt[indC2]; xD2 = xt[indD2];
      yA2 = yt[indA2]; yB2 = yt[indB2]; yC2 = yt[indC2]; yD2 = yt[indD2];
      //zA2 = zt[indA2]; zB2 = zt[indB2]; zC2 = zt[indC2]; zD2 = zt[indD2];
      if (xC2 < xA1-eps) goto nextet2;
      if (xA2 > xC1+eps) goto nextet2;
      if (yC2 < yA1-eps) goto nextet2;
      if (yA2 > yC1+eps) goto nextet2;
      
      //test face A1B1 avec D2C2
      if (K_FUNC::fEqualZero(yA1-yD2, eps) == true)
      {
        if (xC2 > xA1+eps && xD2 < xB1-eps) voisins.push_back(et2); 
        goto nextet2;
      }
      //test face D1C1 avec A2B2
      if (K_FUNC::fEqualZero(yC1-yA2, eps) == true)
      {
        if (xC1 > xA2+eps && xD1 < xB2-eps) voisins.push_back(et2);
        goto nextet2;
      }
      //test A1D1 avec B2C2
      if (K_FUNC::fEqualZero(xA1-xB2, eps) == true)
      {
        if (yC2 > yA1+eps && yB2 < yD1-eps) voisins.push_back(et2);
        goto nextet2;
      }
      //test B1C1 avec A2D2
      if (K_FUNC::fEqualZero(xB1-xA2, eps) == true)
      {
        if (yC1 > yA2+eps && yB1 < yD2-eps) voisins.push_back(et2);
        goto nextet2;
      }
      nextet2:;
    }// fin et2
    nghbrs[et1] = voisins; // copie
  }// fin et1
  }

  delete bbtree;
  for (E_Int i = 0; i < nelts; i++) delete boxes[i];
  return 1;
}
//=============================================================================
// Retourne les voisins par elt ayant une facette en commun
//=============================================================================
E_Int K_GENERATOR::getNeighbourHexas(
  E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
  K_FLD::FldArrayI& cEV, 
  std::vector< std::vector<E_Int> >& nghbrs)
{
  E_Float eps = 1.e-10;
  E_Int nelts = cEV.getSize();
  E_Int* cn1 = cEV.begin(1); //E_Int* cn2 = cEV.begin(2);
  //E_Int* cn3 = cEV.begin(3); E_Int* cn4 = cEV.begin(4);
  //E_Int* cn5 = cEV.begin(5); E_Int* cn6 = cEV.begin(6);
  E_Int* cn7 = cEV.begin(7); //E_Int* cn8 = cEV.begin(8);
  /* 1 - creation du bbtree */
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  std::vector<BBox3DType*> boxes(nelts);// a detruire a la fin
  K_FLD::FldArrayF bbox(nelts,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cEV, xt, yt, zt, bbox);
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);

#pragma omp parallel default(shared)
  {
    E_Float minB0[3];  E_Float maxB0[3];
#pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      minB0[0] = xminp[et]-eps; minB0[1] = yminp[et]-eps; minB0[2] = zminp[et]-eps;
      maxB0[0] = xmaxp[et]+eps; maxB0[1] = ymaxp[et]+eps; maxB0[2] = zmaxp[et]+eps; 
      boxes[et] = new BBox3DType(minB0, maxB0);
    }
  }
  K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);

#pragma omp parallel default(shared)
  {
  E_Int indA1, indG1, indA2, indG2, et2;
  E_Float xmin1, ymin1, zmin1, xmax1, ymax1, zmax1;
  E_Float xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
  E_Float minB[3];  E_Float maxB[3];
#pragma omp for
  for (E_Int et1 = 0; et1 < nelts; et1++)
  {
    indA1 = cn1[et1]-1; indG1 = cn7[et1]-1;
    xmin1 = xt[indA1]; ymin1 = yt[indA1]; zmin1 = zt[indA1]; 
    xmax1 = xt[indG1]; ymax1 = yt[indG1]; zmax1 = zt[indG1]; 
 
    minB[0] = xminp[et1]-eps; minB[1] = yminp[et1]-eps; minB[2] = zminp[et1]-eps;
    maxB[0] = xmaxp[et1]+eps; maxB[1] = ymaxp[et1]+eps; maxB[2] = zmaxp[et1]+eps;
    std::vector<E_Int> indicesBB;
    bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
    // recherche des candidats
    for (size_t noi = 0; noi < indicesBB.size(); noi++)
    {
      et2 = indicesBB[noi];
      indA2 = cn1[et2]-1; indG2 = cn7[et2]-1;
      xmin2 = xt[indA2]; ymin2 = yt[indA2]; zmin2 = zt[indA2]; 
      xmax2 = xt[indG2]; ymax2 = yt[indG2]; zmax2 = zt[indG2]; 
      if (et2 == et1) goto nextet2; 
      if (xmax2 < xmin1-eps) goto nextet2;
      if (xmin2 > xmax2+eps) goto nextet2;
      if (ymax2 < ymin1-eps) goto nextet2;
      if (ymin2 > ymax2+eps) goto nextet2;
      if (zmax2 < zmin1-eps) goto nextet2;
      if (zmin2 > zmax2+eps) goto nextet2;

      //test face ABCD avec EFGH z = cste
      if (K_FUNC::fEqualZero(zmin1-zmax2,eps) == true || K_FUNC::fEqualZero(zmax1-zmin2,eps) == true)
      {
        if (xmin2 < xmax1-eps && ymin2 < ymax1-eps && xmax2 > xmin1+eps && ymax2 > ymin1+eps)
          nghbrs[et1].push_back(et2); 
        goto nextet2;
      }
      //test face ADHE avec BCGF x = cste
      if (K_FUNC::fEqualZero(xmin1-xmax2,eps) == true || K_FUNC::fEqualZero(xmax1-xmin2,eps) == true)
      {
        if ( zmin2 < zmax1-eps && ymin2 < ymax1-eps && zmax2 > zmin1+eps &&  ymax2 > ymin1+eps)
          nghbrs[et1].push_back(et2); 
        goto nextet2;
      }
      //test face ABFE avec DCGH y = cste
      if (K_FUNC::fEqualZero(ymin1-ymax2,eps) == true || K_FUNC::fEqualZero(ymax1-ymin2,eps) == true)
      {
        if (zmin2 < zmax1-eps && xmin2 < xmax1-eps && zmax2 > zmin1+eps && xmax2 > xmin1+eps)
          nghbrs[et1].push_back(et2); 
        goto nextet2;
      }
      nextet2:;
    }// fin et2
  }// fin et1
  } // OMP
  delete bbtree;
  for (E_Int i = 0; i < nelts; i++) delete boxes[i];
  return 1;
}

//=============================================================================
// Retourne les voisins par elt ayant une facette ou un point en commun
// IN: npts: nbre de pts de l'octree
// IN: xt,yt,zt: coord de l'octree
// IN: cEV: connectivite elt/vertex de l'octree
// IN: mindh: pas de la plus petite cellule de l'octree
// OUT: nghbrs: pour chaque elt, cellules voisines.
//=============================================================================
E_Int K_GENERATOR::getNeighbourQuads2(
  E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
  K_FLD::FldArrayI& cEV,
  std::vector< std::vector<E_Int> >& nghbrs, E_Float mindh )
{
  E_Float eps = 0.1*mindh;
  E_Int nelts = cEV.getSize();
  /* 1 - creation du bbtree */
  typedef K_SEARCH::BoundingBox<2>  BBox2DType; 
  std::vector<BBox2DType*> boxes(nelts);
  K_FLD::FldArrayF bbox(nelts,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cEV, xt, yt, zt, bbox); 
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);

#pragma omp parallel default(shared)
  {
    E_Float minB0[2];  E_Float maxB0[2];
#pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      minB0[0] = xminp[et]-eps; minB0[1] = yminp[et]-eps;
      maxB0[0] = xmaxp[et]+eps; maxB0[1] = ymaxp[et]+eps;
      boxes[et] = new BBox2DType(minB0, maxB0);
    }
  }
  K_SEARCH::BbTree2D* bbtree = new K_SEARCH::BbTree2D(boxes);

#pragma omp parallel default(shared)
  {
  E_Float minB[2];  E_Float maxB[2];
  E_Int et2;
#pragma omp for
  for (E_Int et1 = 0; et1 < nelts; et1++)
  {
    minB[0] = xminp[et1]-eps; minB[1] = yminp[et1]-eps; 
    maxB[0] = xmaxp[et1]+eps; maxB[1] = ymaxp[et1]+eps;
    std::vector<E_Int> indicesBB;
    bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
    // recherche des candidats
    for (size_t noi = 0; noi < indicesBB.size(); noi++)
    {
      et2 = indicesBB[noi];
      if (et2 != et1) nghbrs[et1].push_back(et2);
    }// fin et2
  }// fin et1
  }

  for (E_Int et = 0; et < nelts; et++) delete boxes[et];
  delete bbtree;
  return 1;
}
//=============================================================================
// Retourne les voisins par elt ayant une facette ou un point en commun
// IN: npts: nbre de pts de l'octree
// IN: xt,yt,zt: coord de l'octree
// IN: cEV: connectivite elt/vertex de l'octree
// IN: mindh: pas de la plus petite cellule de l'octree
// OUT: nghbrs: pour chaque elt, cellules voisines.
//=============================================================================
E_Int K_GENERATOR::getNeighbourHexas2(
  E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
  K_FLD::FldArrayI& cEV,
  std::vector< std::vector<E_Int> >& nghbrs, E_Float mindh)
{
  E_Float eps = 0.1*mindh;
  E_Int nelts = cEV.getSize();
  /* 1 - creation du bbtree */
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  std::vector<BBox3DType*> boxes(nelts);
  K_FLD::FldArrayF bbox(nelts,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cEV, xt, yt, zt, bbox); 
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);

#pragma omp parallel default(shared)
  {
    E_Float minB0[3]; E_Float maxB0[3];
#pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      minB0[0] = xminp[et]-eps; minB0[1] = yminp[et]-eps; minB0[2] = zminp[et]-eps;
      maxB0[0] = xmaxp[et]+eps; maxB0[1] = ymaxp[et]+eps; maxB0[2] = zmaxp[et]+eps; 
      boxes[et] = new BBox3DType(minB0, maxB0);
    }
  }
  K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);

#pragma omp parallel default(shared)
  {
  E_Float minB[3]; E_Float maxB[3];
  E_Int et2;
#pragma omp for
  for (E_Int et1 = 0; et1 < nelts; et1++)
  {
    minB[0] = xminp[et1]-eps; minB[1] = yminp[et1]-eps; minB[2] = zminp[et1]-eps;
    maxB[0] = xmaxp[et1]+eps; maxB[1] = ymaxp[et1]+eps; maxB[2] = zmaxp[et1]+eps;
    std::vector<E_Int> indicesBB;
    bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
    // recherche des candidats
    for (size_t noi = 0; noi < indicesBB.size(); noi++)
    {
      et2 = indicesBB[noi];
      if (et2 != et1) nghbrs[et1].push_back(et2);
    }// fin et2
  }// fin et1
  }

  for (E_Int et = 0; et < nelts; et++) delete boxes[et];
  delete bbtree;
  return 1;
}
//=============================================================================
/* Modification de l'indicateur de raffinement/deraffinement afin de respecter 
   un rapport de taille de maille 2:1 entre 2 elements adjacents 
   Priorite donnee au niveau fin */
//=============================================================================
void K_GENERATOR::modifyIndicator(
  E_Int nelts, K_FLD::FldArrayF& dht, std::vector< std::vector<E_Int> >& cEEN,
  E_Int posi, E_Float* indict)
{
  E_Float eps = 1.e-10;
  E_Float indic1, indic2, dh1, dh2;
  E_Int et2; E_Int ok = 0;
  E_Float* dhtp = dht.begin(1);
  while (ok == 0)
  {
    ok = 1; 
    for (E_Int et1 = 0; et1 < nelts; et1++)
    {
      dh1 = dhtp[et1];
      std::vector<E_Int>& voisins = cEEN[et1]; E_Int nvoisins = voisins.size();
      for (E_Int noet = 0; noet < nvoisins; noet++)
      {
        et2 = voisins[noet]; dh2 = dhtp[et2];
        indic1 = indict[et1]; indic2 = indict[et2]; 

        if (K_FUNC::fEqual(indic1, indic2) == false) // indic1 != indic2
        {       
          if (indic1 > 0. && K_FUNC::fEqualZero(indic2) == true) 
          {
            if (K_FUNC::fEqualZero(2.*dh1-dh2, eps) == true) 
            { indict[et2] = 1.; ok = 0; }
          }
          else if (K_FUNC::fEqualZero(indic1) == true && indic2 > 0.)
          {
            if (K_FUNC::fEqualZero(dh1-2.*dh2, eps) == true) 
            { indict[et1] = 1.; ok = 0; }
          }
          else if (K_FUNC::fEqualZero(indic1) == true && indic2 < 0.)
          {
            if (K_FUNC::fEqualZero(2.*dh1-dh2, eps) == true) 
            { indict[et2] = 0.; ok = 0; }
          }
          else if (indic1 < 0. && K_FUNC::fEqualZero(indic2) == true)
          {
            if (K_FUNC::fEqualZero(dh1-2.*dh2, eps) == true) 
            { indict[et1] = 0.; ok = 0; }
          }
          else if (indic1 > 0. && indic2 < 0.) 
          {
            if (K_FUNC::fEqualZero(dh1-dh2, eps) == true) 
            { indict[et2] = 0.; ok = 0; }
            else if (K_FUNC::fEqualZero(2.*dh1-dh2, eps) == true) 
            { indict[et2] = 1.; ok = 0; }
          }
          else if (indic1 < 0. && indic2 > 0.)
          {
            if (K_FUNC::fEqualZero(dh1-dh2, eps) == true) 
            { indict[et1] = 0.; ok = 0; }
            else if (K_FUNC::fEqualZero(dh1-2.*dh2, eps) == true) 
            { indict[et1] = 1.; ok = 0; }
          }
        }
      }
    }
  }
}

//=============================================================================
/* Recherche des voisins valides pour la fusion 
   IN: xs, ys, zs: point min de l'octree global
   IN: xe, ye, ze: point max de l'octree global */
//=============================================================================
void K_GENERATOR::getValidNgbrsForMerge(
  E_Int et, E_Float* indict, E_Float* dht, 
  E_Float xs, E_Float ys, E_Float zs,
  E_Float xe, E_Float ye, E_Float ze,
  E_Float* xt, E_Float* yt, E_Float* zt,
  K_FLD::FldArrayIS& dejaVu, K_FLD::FldArrayI& cn,
  std::vector<E_Int>& candidats)
{
  E_Float eps = 1.e-10;
  E_Int nelts = cn.getSize(); E_Int nvert = cn.getNfld();
  E_Float dh = dht[et];
  E_Int* cn1 = cn.begin(1);

  //et de niveau l (fin) determination du niveau l-1 (plus grossier)
  E_Float dh2 = 2*dh;
  E_Int indA = cn1[et]-1; 
  E_Float xmin = xt[indA]; E_Float ymin = yt[indA]; E_Float zmin = zt[indA];
  E_Int i = E_Int((xmin-xs)/dh2+0.1); 
  E_Int j = E_Int((ymin-ys)/dh2+0.1); 
  E_Int k = E_Int((zmin-zs)/dh2+0.1);

  // bbox du niveau l-1 contenant et et ses voisins a fusionner
  xmin = xs + i*dh2; ymin = ys + j*dh2; zmin = zs + k*dh2;
  E_Float xmax = xmin + dh2; E_Float ymax = ymin + dh2; E_Float zmax = zmin + dh2;
  if (nvert == 4) {zmin = 0.; zmax = 1.;}
  // determination des voisins valides de meme niveau et contenus dans le meme niveau -1
  E_Int indA2;
  E_Float xA2, yA2, zA2;
  for (E_Int et2 = 0; et2 < nelts; et2++)
  {
    if (et2 != et && K_FUNC::fEqual(indict[et2], -1.) == true && dejaVu[et2] == 0 && K_FUNC::fEqualZero(dht[et2]-dh,eps) == true)
    {
      indA2 = cn1[et2]-1; xA2 = xt[indA2]; yA2 = yt[indA2]; zA2 = zt[indA2];
      if (xA2 > xmin-eps && xA2 < xmax-eps && 
          yA2 > ymin-eps && yA2 < ymax-eps &&
          zA2 > zmin-eps && zA2 < zmax-eps)
        candidats.push_back(et2);
    }    
  }
}
//=============================================================================
/* Fusionne un element avec 3/8 voisins valides */
//=============================================================================
E_Int
K_GENERATOR::mergeOctreeElement(E_Int et, E_Int npts, E_Float indic,
                                K_FLD::FldArrayI& cn, 
                                E_Float xs, E_Float ys, E_Float zs,
                                E_Float xe, E_Float ye, E_Float ze,
                                E_Float* xt, E_Float* yt, E_Float* zt, 
                                E_Float* dht, E_Float* indict, 
                                K_FLD::FldArrayI& cno, K_FLD::FldArrayF& fo, 
                                K_FLD::FldArrayF& indicout,
                                E_Int& no, E_Int& eto, 
                                K_FLD::FldArrayIS& dejaVu)
{
  E_Float eps = 1.e-10;
  // Recherche des elements voisins de meme niveau
  E_Int nvert = cno.getNfld(); E_Int nelts = cno.getSize();
  E_Int nptso = fo.getSize();
  if (no >= nptso-9*nvert) 
  {nptso = nptso+npts; fo.reAllocMat(nptso,3);}
  if (eto >= cno.getSize()-17) 
  {nelts = cno.getSize()+nelts; cno.reAllocMat(nelts,nvert); indicout.reAlloc(nelts); }

  E_Float* xo = fo.begin(1); E_Float* yo = fo.begin(2); E_Float* zo = fo.begin(3);
  E_Float* indico = indicout.begin();
  E_Float dh = dht[et];
  E_Float xmin, ymin, xmax, ymax, zmin, zmax;

  std::vector<E_Int> candidats;
  getValidNgbrsForMerge(et, indict, dht, xs, ys, zs, xe, ye, ze, 
                        xt, yt, zt, dejaVu, cn, candidats);
  E_Int ncandidats = candidats.size();
  // 1- verifier si on a au moins un nb de candidats suffisant (3 pour le 2D).
  // 2- sinon: 3 QUAD candidats: forment-ils un carre ? Immediat 
  // 3- plus de nmax candidats: impossible, sinon ca sent le bug qqpart...
  E_Int nmax = nvert-1;
  E_Int indA1, indB1, indC1, indE1;
  E_Int indA2, indB2, indC2, indE2;
  
  if (ncandidats < nmax) 
  {
    return -1; //cas 1
  }
  else if (ncandidats == nmax) //cas 2
  {
    if (nvert == 4) // QUAD
    {
      E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3); //E_Int* cn4 = cn.begin(4);
      indA1 = cn1[et]-1; indB1 = cn2[et]-1; indC1 = cn3[et]-1;
      xmin = xt[indA1]; xmax = xt[indB1]; 
      ymin = yt[indA1]; ymax = yt[indC1];
      zmin = zt[indA1]; zmax = zmin;
      for (E_Int v = 0; v < ncandidats; v++)
      {
        E_Int et2 = candidats[v]; 
        indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1; //indD2 = cn4[et2]-1;
        xmin = K_FUNC::E_min(xmin, xt[indA2]); xmax = K_FUNC::E_max(xmax, xt[indB2]); 
        ymin = K_FUNC::E_min(ymin, yt[indA2]); ymax = K_FUNC::E_max(ymax, yt[indC2]); 
      }
      if (K_FUNC::fEqualZero(2*dh-(xmax-xmin),eps) == false || 
          K_FUNC::fEqualZero(2*dh-(ymax-ymin),eps) == false ) return -2;

      // construction du quad
      xo[no] = xmin; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,1) = no; 
      xo[no] = xmax; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,2) = no; 
      xo[no] = xmax; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,3) = no; 
      xo[no] = xmin; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,4) = no; 
      indico[eto] = K_FUNC::E_min(0.,indic+1);
      eto++;
      for (E_Int v = 0; v < ncandidats; v++)
      { E_Int et2 = candidats[v]; dejaVu[et2] = 1; }
      dejaVu[et] = 1;
      return 1;
    }
    else // HEXA
    {
      E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2); 
      E_Int* cn3 = cn.begin(3); //E_Int* cn4 = cn.begin(4);
      E_Int* cn5 = cn.begin(5); //E_Int* cn6 = cn.begin(6); 
      //E_Int* cn7 = cn.begin(7); E_Int* cn8 = cn.begin(8);
      indA1 = cn1[et]-1; indB1 = cn2[et]-1; indC1 = cn3[et]-1; //indD1 = cn4[et]-1;
      indE1 = cn5[et]-1; //indF1 = cn6[et]-1; indG1 = cn7[et]-1; indH1 = cn8[et]-1;
      xmin = xt[indA1]; xmax = xt[indB1]; 
      ymin = yt[indA1]; ymax = yt[indC1];
      zmin = zt[indA1]; zmax = zt[indE1];
    
      for (E_Int v = 0; v < ncandidats; v++)
      {
        E_Int et2 = candidats[v]; 
        indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1; indE2 = cn5[et2]-1;
        xmin = K_FUNC::E_min(xmin, xt[indA2]); xmax = K_FUNC::E_max(xmax, xt[indB2]); 
        ymin = K_FUNC::E_min(ymin, yt[indA2]); ymax = K_FUNC::E_max(ymax, yt[indC2]);
        zmin = K_FUNC::E_min(zmin, zt[indA2]); zmax = K_FUNC::E_max(zmax, zt[indE2]);
      }
      
      if (K_FUNC::fEqualZero(2*dh-(xmax-xmin),eps) == false || 
          K_FUNC::fEqualZero(2*dh-(ymax-ymin),eps) == false ||
          K_FUNC::fEqualZero(2*dh-(zmax-zmin),eps) == false) return -3;
  
      // construction de l hexa
      xo[no] = xmin; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,1) = no; 
      xo[no] = xmax; yo[no] = ymin; zo[no] = zmin; no++; cno(eto,2) = no; 
      xo[no] = xmax; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,3) = no; 
      xo[no] = xmin; yo[no] = ymax; zo[no] = zmin; no++; cno(eto,4) = no; 
      xo[no] = xmin; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,5) = no; 
      xo[no] = xmax; yo[no] = ymin; zo[no] = zmax; no++; cno(eto,6) = no; 
      xo[no] = xmax; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,7) = no; 
      xo[no] = xmin; yo[no] = ymax; zo[no] = zmax; no++; cno(eto,8) = no; 
      indico[eto] = K_FUNC::E_min(0.,indic+1);
      eto++;
      for (E_Int v = 0; v < ncandidats; v++)
      { E_Int et2 = candidats[v]; dejaVu[et2] = 1; }
      dejaVu[et] = 1;
      return 1;
    }
  }
  printf("Error: mergeElement: too many candidates (" SF_D_ ") to merge with element " SF_D_ ". Check the octree.\n",
         ncandidats, et+1);
  
  return -1;
}
