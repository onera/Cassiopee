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
// coarsen elements of a surface triangular mesh by edge contraction

# include <math.h>
# include <string.h>
# include <stdio.h>
# include "post.h"

using namespace std;
using namespace K_FLD;


//===========================================================================
/* Fusion des elements d'un maillage surfacique TRI en fonction d'un critere
   de deraffinement.
   Retourne un maillage surfacique TRI modifie */
//===========================================================================
PyObject* K_POST::coarsen(PyObject* self, PyObject* args)
{
  // surf: maillage a deraffiner (x,y,z+sol)
  // indic: indicateur de deraffinement: vaut 0 ou 1
  PyObject* surf; PyObject* aindic;
  E_Float eps, argqual;
  
  if (!PYPARSETUPLE_(args, OO_ RR_,
                    &surf, &aindic, &argqual, &eps))
  {
    return NULL;
  }
  
  // check argqual: between 0 and 0.5
  if (argqual > 0.5 or argqual < 0.)
  {
    printf("Warning: coarsen: argqual must be between 0 and 0.5. Set to default value: 0.25.\n");
    argqual = 0.25;
  }
  /*-----------------------------------------------*/
  /* Extraction des donnees du maillage surfacique */ 
  /*-----------------------------------------------*/
  char* varString0; char* eltType0;
  FldArrayF* f; FldArrayI* cn;
  E_Int nil, njl, nkl;
  E_Int res = 
    K_ARRAY::getFromArray3(surf, varString0, f, nil, njl, nkl, cn, eltType0);
  
  if (res != 2 || strcmp(eltType0, "TRI") != 0)
  {
    RELEASESHAREDB(res, surf, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "coarsen: array is invalid.");
    return NULL;
  }
  
  // Check size of array
  E_Int posxu = K_ARRAY::isCoordinateXPresent(varString0);
  E_Int posyu = K_ARRAY::isCoordinateYPresent(varString0);
  E_Int poszu = K_ARRAY::isCoordinateZPresent(varString0);

  if (posxu == -1 || posyu == -1 || poszu == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "coarsen: array must contain coordinates.");
    RELEASESHAREDU(surf, f, cn);
    return NULL;
  }
  posxu++; posyu++; poszu++;
  
  /*-------------------------------------------*/
  /* Extraction de l'indicateur de raffinement */ 
  /*-------------------------------------------*/
  char* varString1;
  FldArrayF* findic;
  E_Int ni1, nj1, nk1;
  FldArrayI* cn1 = NULL;
  char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray3(aindic, varString1, findic, 
                                      ni1, nj1, nk1, cn1, eltType1);
  E_Int nelts = cn->getSize();
  if (res1 == 1 )
  {
    if (ni1*nj1*nk1 != nelts)
    {
      RELEASESHAREDU(surf, f, cn);
      RELEASESHAREDS(aindic, findic);
      PyErr_SetString(
        PyExc_TypeError,
        "coarsen: dimension of refinement indicator array must be equal to the number of elements.");
      return NULL;
    }
  }
  else if (res1 == 2)
  {
    if ( findic->getSize() != nelts )
    {
      RELEASESHAREDU(surf, f, cn);
      RELEASESHAREDU(aindic, findic, cn1);
      PyErr_SetString(
        PyExc_TypeError,
        "coarsen: dimension of refinement indicator array must be equal to the number of elements.");
      return NULL;
    }
    delete cn1;
  }
  else
  {
    RELEASESHAREDU(surf, f, cn);  
    PyErr_SetString(PyExc_TypeError,
                    "coarsen: refinement indicator array is invalid.");
    return NULL;
  }

  E_Int api = f->getApi();
  FldArrayIS indic(nelts);
  short* indicp = indic.begin();
  E_Float* findicp = findic->begin(); 
  for (E_Int i = 0; i < nelts; i++) indicp[i] = short(findicp[i]);
  delete findic;

  /* fusion des elements */  
  mergeElements(*cn, *f, posxu ,posyu, poszu, argqual, indic, eps);

  /* retour */ 
  PyObject* t = K_ARRAY::buildArray3(*f, varString0, *cn, "TRI", api);
  RELEASESHAREDU(surf, f, cn);
    
  return t;
}

//==========================================================================
/* Fusion des elements */
//==========================================================================
void K_POST::mergeElements(FldArrayI& connect, FldArrayF& field,
                           E_Int posx, E_Int posy, E_Int posz,
                           E_Float argqual,
                           FldArrayIS& indic, E_Float& eps)
{
  // 0-determination des pts externes
  E_Int npts = field.getSize();
  FldArrayI extNodes(npts);
  getExternNodes(connect, field, extNodes);

  E_Float* xt = field.begin(posx);
  E_Float* yt = field.begin(posy);
  E_Float* zt = field.begin(posz);

  // 1-calcul de la connectivite vertex->elements
  vector< vector<E_Int> > cVE(npts); 
  K_CONNECT::connectEV2VE(connect, cVE);

  // 2-distance de projection 
  E_Int nelts = connect.getSize();
  E_Int* cn1 = connect.begin(1);
  E_Int* cn2 = connect.begin(2);
  E_Int* cn3 = connect.begin(3);

  FldArrayI selectedVertices(6, 2);
  
  for (E_Int et = 0; et < nelts; et++)
  { 
    selectedVertices.setAllValuesAt(-1);
    E_Int c = 0;
    if (isDegenerated(et, connect, field.begin(posx), 
                      field.begin(posy), field.begin(posz)) == 0 && 
         indic[et] == 1) 
    {
      E_Int test, indo1, indo2;
      // determination des combinaisons de fusion
      E_Int ind1 = cn1[et]-1; E_Int ind2 = cn2[et]-1; E_Int ind3 = cn3[et]-1;
            
      if (extNodes[ind1] == 0) 
      {
        test = testFusion(et, ind1, ind2, connect, indic, 
                          xt, yt, zt, cVE, eps, 
                          indo1, indo2);
        if (test == 1) {selectedVertices(c,1)=indo1; selectedVertices(c,2)=indo2; c++;}
        test = testFusion(et, ind1, ind3, connect, indic, 
                          xt, yt, zt, cVE, eps,
                          indo1, indo2); 
        if (test == 1) {selectedVertices(c,1)=indo1; selectedVertices(c,2)=indo2; c++;}
      }
      if (extNodes[ind2] == 0) 
      {
        test = testFusion(et, ind2, ind1, connect, indic, 
                          xt, yt, zt, cVE, eps, 
                          indo1, indo2); 
        if (test == 1) {selectedVertices(c,1)=indo1; selectedVertices(c,2)=indo2; c++;}
        test = testFusion(et, ind2, ind3, connect, indic, 
                          xt, yt, zt, cVE, eps,
                          indo1, indo2); 
        if (test == 1) {selectedVertices(c,1)=indo1; selectedVertices(c,2)=indo2; c++;}
      }
      if (extNodes[ind3] == 0) 
      {
        test = testFusion(et, ind3, ind1, connect, indic, 
                          xt, yt, zt, cVE, eps,
                          indo1, indo2);  
        if (test == 1) {selectedVertices(c,1)=indo1; selectedVertices(c,2)=indo2; c++;}
        test = testFusion(et, ind3, ind2, connect, indic, 
                          xt, yt, zt, cVE, eps,
                          indo1, indo2); 
        if (test == 1) {selectedVertices(c,1)=indo1; selectedVertices(c,2)=indo2; c++;}
      }
      

      if ( c > 0 )
      {
        if (c == 1) { indo1 = selectedVertices(0,1); indo2 = selectedVertices(0,2);}
        else 
          selectBestCandidatesForMerge(et, selectedVertices, connect, 
                                       argqual, xt, yt, zt, cVE, indo1, indo2);

        if (indo1 != -1 && indo2 != -1) 
        {
          xt[indo1] = xt[indo2]; yt[indo1] = yt[indo2]; zt[indo1] = zt[indo2];
 
          for (E_Int et2 = 0; et2 < nelts; et2++)
          {
            if (cn1[et2] == indo1+1) {cn1[et2] = indo2+1; indic[et2] = 0;}
            if (cn2[et2] == indo1+1) {cn2[et2] = indo2+1; indic[et2] = 0;}
            if (cn3[et2] == indo1+1) {cn3[et2] = indo2+1; indic[et2] = 0;}
          }
        }
      }
    }//test non degenere
  }//parcours de ts les elts
  
  // 3- close
  E_Float tolc = 1.e-6; 
  K_CONNECT::cleanConnectivity(posx, posy, posz, tolc, 
                               "TRI", field, connect);
}

//=============================================================================
// fusion des pts
//=============================================================================
E_Int K_POST::testFusion(E_Int et, E_Int ind1, E_Int ind2, FldArrayI& connect,
                         FldArrayIS& indic,
                         E_Float* xt, E_Float* yt, E_Float* zt, 
                         vector< vector<E_Int> >& cVE, E_Float eps,
                         E_Int& indo1, E_Int& indo2)
{
  E_Float eps2 = eps*eps;
  indo1 = -1; indo2 = -1;

  E_Int* cn1 = connect.begin(1);
  E_Int* cn2 = connect.begin(2);
  E_Int* cn3 = connect.begin(3);
  E_Float xp, yp, zp, dist2;
  E_Float sigma0, sigma1;
  E_Float p1[3]; E_Float p2[3]; E_Float p3[3]; E_Float p[3];
  vector<E_Int>& cVE1 = cVE[ind1];
  E_Int neti = cVE1.size();

  for (E_Int noeti = 0; noeti < neti; noeti++)
  {
    E_Int eti = cVE1[noeti];
    if (eti != et) 
    {
      if (indic[eti] != 1) return 0;

      E_Int indt1 = cn1[eti]-1; 
      E_Int indt2 = cn2[eti]-1; 
      E_Int indt3 = cn3[eti]-1;
      if (indt1 == ind1) indt1 = ind2;
      if (indt2 == ind1) indt2 = ind2;
      if (indt3 == ind1) indt3 = ind2;
      
      p1[0] = xt[indt1]; p1[1] = yt[indt1]; p1[2] = zt[indt1];
      p2[0] = xt[indt2]; p2[1] = yt[indt2]; p2[2] = zt[indt2];
      p3[0] = xt[indt3]; p3[1] = yt[indt3]; p3[2] = zt[indt3];
      p[0] = xt[ind1]; p[1] = yt[ind1]; p[2] = zt[ind1]; 
      
      // 1 - verification critere de distance
      E_Bool in = false;
      K_COMPGEOM::distanceToTriangle(
        p1, p2, p3, p, 1,
        dist2, in, xp, yp, zp, sigma0, sigma1);
      if (dist2 < eps2)   
      {
        //2 - verification de la convexite de la galette
        E_Int convex = isConvex(ind1, connect, cVE1, xt, yt, zt);
        if (convex == 1) {indo1 = ind1; indo2 = ind2;}
        else return 0;
      }
      else return 0;
    }// parcours des voisins du pt 1
  }
  return 1;
}
//=========================================================================
/* Calcul de la convexite de la galette */
//=========================================================================
E_Int K_POST::isConvex(E_Int indA, FldArrayI& connect, vector<E_Int> & cVE,
                       E_Float* xt, E_Float* yt, E_Float* zt)
{
  E_Float eps  = 1e-4;//tolerance sur le sinus de l angle teta1+teta2

  E_Int ntri = cVE.size();// nb de triangles ds la galette
  FldArrayIS tag(ntri);
  tag.setAllValuesAtNull();

  E_Int* cn1 = connect.begin(1);
  E_Int* cn2 = connect.begin(2);
  E_Int* cn3 = connect.begin(3);
  
  for (E_Int noet1 = 0; noet1 < ntri; noet1++)
  {
    if ( tag[noet1] < 2 ) 
    {
      E_Int et1 = cVE[noet1];
      E_Int indB1 = -1; E_Int indC1 = -1;

      if ( cn1[et1]-1 == indA) 
      {indB1 = cn2[et1]-1; indC1 = cn3[et1]-1;}
      else if ( cn2[et1]-1 == indA) 
      {indB1 = cn1[et1]-1; indC1 = cn3[et1]-1;}
      else //if ( cn3[et1]-1 == indA)
      {indB1 = cn1[et1]-1; indC1 = cn2[et1]-1;}
      
      E_Int indB2 = -1; E_Int indC2 = -1;
      E_Int indB = -1; E_Int indC = -1; E_Int indD = -1;

      for (E_Int noet2 = 0; noet2 < ntri; noet2++)
      {
        if ( noet2 != noet1 && tag[noet2] < 2)
        {
          E_Int et2 = cVE[noet2];
          if ( cn1[et2]-1 == indA) 
          {indB2 = cn2[et2]-1; indC2 = cn3[et2]-1;}
          else if ( cn2[et2]-1 == indA) 
          {indB2 = cn1[et2]-1; indC2 = cn3[et2]-1;}
          else //if ( cn3[et2]-1 == indA)
          {indB2 = cn1[et2]-1; indC2 = cn2[et2]-1;}

          if ( indB1 == indB2)
          {indB = indB1; indC = indC1; indD = indC2;}
          else if ( indB1 == indC2)
          {indB = indB1; indC = indC1; indD = indB2;}
          else if ( indC1 == indB2)
          {indB = indC1; indC = indB1; indD = indC2;}
          else if ( indC1 == indC2)
          {indB = indC1; indC = indB1; indD = indB2;}
          else goto next2;
          //ABC : triangle et1 - ABD : triangle et2
          E_Float t1pt2 = computeAngle(indA, indB, indC, indD, xt, yt, zt);

          // t1+t2 > 180 deg ->concavite trouvee
          if (t1pt2 > 180 - eps) return 0;
          
          // maj du tag : nb de traitmts effectues par triangle
          tag[noet1] = tag[noet1]+1;
          tag[noet2] = tag[noet2]+1;
          if (tag[noet1] == 2) goto next1;//chq triangle est traite 2 fois maxi
          
        }
        next2:;//triangle et2 suivant
      }
    }                                           
    next1:;///triangle et1 suivant
  }
  return 1;
}
//===========================================================================
/* Calcul de teta1+teta2: teta1 = (BC,BA), teta2 = (BA,BD) */
//===========================================================================
E_Float K_POST::computeAngle(E_Int indA, E_Int indB, E_Int indC, E_Int indD, 
                             E_Float* xt, E_Float* yt, E_Float* zt)
{
  E_Float rad2deg = K_CONST::E_PI_DEG/K_CONST::E_PI;
  E_Float xBC = xt[indC]-xt[indB];
  E_Float yBC = yt[indC]-yt[indB];
  E_Float zBC = zt[indC]-zt[indB];

  E_Float xBA = xt[indA]-xt[indB];
  E_Float yBA = yt[indA]-yt[indB];
  E_Float zBA = zt[indA]-zt[indB];

  E_Float xBD = xt[indD]-xt[indB];
  E_Float yBD = yt[indD]-yt[indB];
  E_Float zBD = zt[indD]-zt[indB];

  E_Float dBC = xBC*xBC+yBC*yBC+zBC*zBC;
  E_Float dBA = xBA*xBA+yBA*yBA+zBA*zBA;
  E_Float dBD = xBD*xBD+yBD*yBD+zBD*zBD;
  E_Float dBCdBA = sqrt(dBC*dBA);
  E_Float dBAdBD = sqrt(dBA*dBD);

  E_Float ps1 = xBC*xBA + yBC*yBA + zBC*zBA;//<BC,BA>
  E_Float ps2 = xBA*xBD + yBA*yBD + zBA*zBD;//<BA,BD>

  E_Float teta1 = acos(ps1/dBCdBA)*rad2deg;
  E_Float teta2 = acos(ps2/dBAdBD)*rad2deg;
  return teta1+teta2;
}

//=========================================================================
/* Retourne 0 si triangle non degenere, 1 sinon */
//=========================================================================
E_Int K_POST::isDegenerated(E_Int noet, FldArrayI& cn, 
                            E_Float* xt, E_Float* yt, E_Float* zt)
{
  E_Int ind1 = cn(noet,1)-1; 
  E_Int ind2 = cn(noet,2)-1; 
  E_Int ind3 = cn(noet,3)-1; 

  E_Float dx = xt[ind2]-xt[ind1];
  E_Float dy = yt[ind2]-yt[ind1];
  E_Float dz = zt[ind2]-zt[ind1];

  if ( K_FUNC::fEqualZero(dx) == true &&
       K_FUNC::fEqualZero(dy) == true && 
       K_FUNC::fEqualZero(dz) == true ) return 1;
  
  dx = xt[ind3]-xt[ind1];
  dy = yt[ind3]-yt[ind1];
  dz = zt[ind3]-zt[ind1];
  if ( K_FUNC::fEqualZero(dx) == true &&
       K_FUNC::fEqualZero(dy) == true && 
       K_FUNC::fEqualZero(dz) == true ) return 1;

  dx = xt[ind2]-xt[ind3];
  dy = yt[ind2]-yt[ind3];
  dz = zt[ind2]-zt[ind3];
  if ( K_FUNC::fEqualZero(dx) == true &&
       K_FUNC::fEqualZero(dy) == true && 
       K_FUNC::fEqualZero(dz) == true ) return 1;
  
  return 0;
}
//==========================================================================
/* Determination des sommets externes. 
   extNodes vaut 1 si sommet externe, 0 sinon */
//==========================================================================
void K_POST::getExternNodes(FldArrayI& cn, FldArrayF& coord,
                            FldArrayI& extNodes)
{
  E_Int nfaces = 3; E_Int nvert = 2;
  FldArrayF fext; 
  FldArrayI cnext;
  exteriorFacesBasic(nfaces, nvert, coord, cn, fext, cnext); 
  K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-6, "TRI", fext, cnext);
  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);
    
  E_Float* xb = fext.begin(1);
  E_Float* yb = fext.begin(2);
  E_Float* zb = fext.begin(3);
  
  E_Float dx, dy, dz;
  extNodes.setAllValuesAtNull();
  
  for (E_Int ind1 = 0; ind1 < coord.getSize(); ind1++)
  {
    for (E_Int inde = 0; inde < fext.getSize(); inde++)
    {
      dx = xt[ind1]-xb[inde];
      dy = yt[ind1]-yb[inde];
      dz = zt[ind1]-zb[inde];
      if (K_FUNC::fEqualZero(dx) == true && 
          K_FUNC::fEqualZero(dy) == true &&
          K_FUNC::fEqualZero(dz) == true)
      {
        extNodes[ind1] = 1;
        goto next;
      }
    }
    next:;
  }
}

//=============================================================================
/* Selection de la meilleure galette au sens ou les triangles obtenus sont 
   de la meilleure qualite parmi les candidats*/
//=============================================================================
void K_POST::selectBestCandidatesForMerge(
  E_Int et0, FldArrayI& selectedVertices,
  FldArrayI& connect, E_Float argqual,
  E_Float* xt, E_Float* yt, E_Float* zt, 
  vector< vector<E_Int> >& cVE, 
  E_Int& indo1, E_Int& indo2)
{
  indo1 = -1; indo2 = -1;
  E_Float p1[3]; E_Float p2[3]; E_Float p3[3];

  E_Int* cn1 = connect.begin(1);
  E_Int* cn2 = connect.begin(2);
  E_Int* cn3 = connect.begin(3);
  E_Float best = argqual;
  
  // parcours des possibilites de merges possibles
  for (E_Int c = 0; c < selectedVertices.getSize(); c++)
  {
    E_Int ind1 = selectedVertices(c,1);
    E_Int ind2 = selectedVertices(c,2);
    if ( ind1 != -1 && ind2 != -1 )
    {
      vector<E_Int>& cVE1 = cVE[ind1];
      E_Int neti = cVE1.size();

      E_Float minratio = 1.;
      //parcours des triangles de la galette de ind1
      for (E_Int noeti = 0; noeti < neti; noeti++)
      {
        E_Int eti = cVE1[noeti];
        if ( eti != et0 ) 
        {
          E_Int indt1 = cn1[eti]-1; 
          E_Int indt2 = cn2[eti]-1; 
          E_Int indt3 = cn3[eti]-1;
          if ( indt1 == ind1 ) indt1 = ind2;
          if ( indt2 == ind1 ) indt2 = ind2;
          if ( indt3 == ind1 ) indt3 = ind2;
          
          p1[0] = xt[indt1]; p1[1] = yt[indt1]; p1[2] = zt[indt1];
          p2[0] = xt[indt2]; p2[1] = yt[indt2]; p2[2] = zt[indt2];
          p3[0] = xt[indt3]; p3[1] = yt[indt3]; p3[2] = zt[indt3];
          E_Float r1 = K_COMPGEOM::inscribedCircleRadius(p1, p2, p3);
          E_Float r2 = K_COMPGEOM::circumCircleRadius(p1[0], p1[1], p1[2], 
                                                      p2[0], p2[1], p2[2],
                                                      p3[0], p3[1], p3[2]);
          if (K_FUNC::fEqualZero(r2) == false)
          { E_Float crit = r1/r2; if (crit < minratio) minratio = crit; }
        }
      }//parcours de ts les triangles de la triangulation candidate
      if (minratio > best)
      {
        best = minratio; indo1 = ind1; indo2 = ind2;
      }
    }
  }
}
