/*    
    Copyright 2013-2019 Onera.

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

#include "BlkInterp.h"

# include <list>
# include <stack>

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

struct stackData
{
  K_KINTERP::BlkIntTreeNode* current;
  E_Float xxmax,yymax,zzmax,ppmax,qqmax,rrmax;
  E_Float xxmin,yymin,zzmin,ppmin,qqmin,rrmin;
  E_Int coupure;
};

//=============================================================================
/* Recherche de la liste des cellules candidates. Retourne la taille de 
   listOfCandidateCells */
//=============================================================================
E_Int K_KINTERP::BlkInterpAdt::
getListOfCandidateCells(E_Float x, E_Float y, E_Float z,
                        list<E_Int>& listOfCandidateCells)
{
  if (x < _xmin) return 0;
  if (y < _ymin) return 0;
  if (z < _zmin) return 0;
  if (x > _xmax) return 0;
  if (y > _ymax) return 0;
  if (z > _zmax) return 0;
  E_Float a1, a2, a3, a4, a5, a6;
  E_Float b1, b2, b3, b4, b5, b6;
  E_Float xmax, ymax, zmax, xmin, ymin, zmin;
  stack<stackData> stack;
  stackData dataForStack;
  
  a1 = _xmin;
  a2 = _ymin;
  a3 = _zmin;
  a4 = x;
  a5 = y;
  a6 = z;
   
  b1 = x;
  b2 = y;
  b3 = z;
  b4 = _xmax;
  b5 = _ymax;
  b6 = _zmax;

  K_KINTERP::BlkIntTreeNode* current;
  current = _tree;

  /* Get down the tree */
  E_Int coupure = 0;
  E_Int coupure2;
   
  E_Float xxmax, yymax, zzmax, xxmin, yymin, zzmin;
  E_Float ppmax, qqmax, rrmax, ppmin, qqmin, rrmin;

  xxmin = _xmin;    // boite de recherche
  yymin = _ymin;
  zzmin = _zmin;
  ppmin = _xmin;
  qqmin = _ymin;
  rrmin = _zmin;
  xxmax = _xmax;
  yymax = _ymax;
  zzmax = _zmax;
  ppmax = _xmax;
  qqmax = _ymax;
  rrmax = _zmax;

  dataForStack.current = current;
  dataForStack.xxmax = xxmax;
  dataForStack.yymax = yymax;
  dataForStack.zzmax = zzmax;
  dataForStack.ppmax = ppmax;
  dataForStack.qqmax = qqmax;
  dataForStack.rrmax = rrmax;
  dataForStack.xxmin = xxmin;
  dataForStack.yymin = yymin;
  dataForStack.zzmin = zzmin;
  dataForStack.ppmin = ppmin;
  dataForStack.qqmin = qqmin;
  dataForStack.rrmin = rrmin;
  dataForStack.coupure = coupure;
  
  stack.push(dataForStack);
   
  while (stack.size() != 0)
  {
    dataForStack = stack.top();
    stack.pop();

    current = dataForStack.current;
    xxmin = dataForStack.xxmin;
    yymin = dataForStack.yymin;
    zzmin = dataForStack.zzmin;
    ppmin = dataForStack.ppmin;
    qqmin = dataForStack.qqmin;
    rrmin = dataForStack.rrmin;
    xxmax = dataForStack.xxmax;
    yymax = dataForStack.yymax;
    zzmax = dataForStack.zzmax;
    ppmax = dataForStack.ppmax;
    qqmax = dataForStack.qqmax;
    rrmax = dataForStack.rrmax;
    coupure = dataForStack.coupure;

    coupure2 = coupure+1;
    if (coupure2 > 5) coupure2 = 0;
    
    /* examine node */
    current->getCellBB(xmax, ymax, zmax, xmin, ymin, zmin);
     
    if ( xmin>=a1 && xmin<=b1 &&
         ymin>=a2 && ymin<=b2 &&
         zmin>=a3 && zmin<=b3 &&
         xmax>=a4 && xmax<=b4 &&
         ymax>=a5 && ymax<=b5 &&
         zmax>=a6 && zmax<=b6 )
    {
      listOfCandidateCells.push_back(current->_ind);
    }
      
    switch (coupure)
    {
      // dichotomie direction 1
      case 0:
        xxmax = K_CONST::ONE_HALF*(xxmax+xxmin);
          
        if (a1 <= xxmax)
        {
          dataForStack.current = current->_left;
          dataForStack.xxmax = xxmax;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
            
        xxmax = 2*xxmax-xxmin;
        xxmin = K_CONST::ONE_HALF*(xxmax+xxmin);
        if (b1 >= xxmin)
        {
          dataForStack.current = current->_right;
          dataForStack.xxmax = xxmax;
          dataForStack.xxmin = xxmin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
            
        break;
            
        // dichotomie direction 2  
      case 1:
          
        yymax = K_CONST::ONE_HALF*(yymax+yymin);
          
        if (a2 <= yymax)
        {
          dataForStack.current = current->_left;
          dataForStack.yymax = yymax;
           dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        yymax = 2*yymax-yymin;
        yymin = K_CONST::ONE_HALF*(yymax+yymin);
        if (b2 >= yymin)
        {
          dataForStack.current = current->_right;
          dataForStack.yymax = yymax;
          dataForStack.yymin = yymin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
            
        break;

        // dichotomie direction 3     
      case 2:
            
        zzmax = K_CONST::ONE_HALF*(zzmax+zzmin);
          
        if (a3 <= zzmax)
        {
          dataForStack.current = current->_left;
          dataForStack.zzmax = zzmax;
           dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        zzmax = 2*zzmax-zzmin;
        zzmin = K_CONST::ONE_HALF*(zzmax+zzmin);
        if (b3 >= zzmin)
        {
          dataForStack.current = current->_right;
          dataForStack.zzmax = zzmax;
          dataForStack.zzmin = zzmin;
           dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        break;
          
        // dichotomie direction 4     
      case 3:
          
        ppmax = K_CONST::ONE_HALF*(ppmax+ppmin);
          
        if (a4 <= ppmax)
        {
          dataForStack.current = current->_left;
          dataForStack.ppmax = ppmax;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        ppmax = 2*ppmax-ppmin;
        ppmin = K_CONST::ONE_HALF*(ppmax+ppmin);
        if (b4 >= ppmin)
        {
          dataForStack.current = current->_right;
          dataForStack.ppmax = ppmax;
          dataForStack.ppmin = ppmin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        break;
          
        // dichotomie direction 5     
      case 4:
          
        qqmax = K_CONST::ONE_HALF*(qqmax+qqmin);
          
        if (a5 <= qqmax)
        {
          dataForStack.current = current->_left;
          dataForStack.qqmax = qqmax;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        qqmax = 2*qqmax-qqmin;
        qqmin = K_CONST::ONE_HALF*(qqmax+qqmin);
        if (b5 >= qqmin)
        {
          dataForStack.current = current->_right;
          dataForStack.qqmax = qqmax;
          dataForStack.qqmin = qqmin;
          dataForStack.rrmin = rrmin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        break;
          
        // dichotomie direction 6     
      case 5:
          
        rrmax = K_CONST::ONE_HALF*(rrmax+rrmin);
          
        if (a6 <= rrmax)
        {
          dataForStack.current = current->_left;
          dataForStack.rrmax = rrmax;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        rrmax = 2*rrmax-rrmin;
        rrmin = K_CONST::ONE_HALF*(rrmax+rrmin);
        if (b6 >= rrmin)
        {
          dataForStack.current = current->_right;
          dataForStack.rrmax = rrmax;
          dataForStack.rrmin = rrmin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        break;
    } /* switch */      
  }   
  E_Int size = listOfCandidateCells.size();
  return size;
}

// ============================================================================
/* Find the interpolation cell for point (x,y,z) */
// ============================================================================
short K_KINTERP::BlkInterpAdt::searchInterpolationCell(
  E_Float x, E_Float y, E_Float z,
  E_Int &ic, E_Int &jc, E_Int &kc,
  FldArrayF& cf)
{
  ic = 0; jc = 0; kc = 0;
  K_KINTERP::KMesh& mesh = const_cast<K_KINTERP::KMesh&>(_mesh);

  // recherche de la liste des cellules candidates
  list<E_Int> listOfCandidateCells; 

  E_Int found = getListOfCandidateCells(x, y, z, listOfCandidateCells);
 
 /* Find the right cell among candidate cells */
  if (found == 0) return 0; //listOfCandidateCells vide
  
  list<E_Int>::iterator itr;

  E_Boolean JUMP;

  E_Int ind;
  E_Int i;
  E_Int isomm;
  E_Int is, js, ks;

  E_Float xi, yi, zi;
  E_Float xt[15], yt[15], zt[15];

  /* Get informations on the KMesh on which is built the InterpCartGrid */
  E_Float* xl = mesh.getXVector();
  E_Float* yl = mesh.getYVector();
  E_Float* zl = mesh.getZVector();

  E_Int ni = mesh.getIm();
  E_Int nj = mesh.getJm();
  E_Int nk = mesh.getKm();
  E_Int nij = ni*nj;
  
  /* For 8 cells, apply the technique of jump to find interpolation cell */
  JUMP = false;
  itr = listOfCandidateCells.begin();
  
  ind = *itr ;

  for (i = 0; i < 8; i++)
  { 
    coordHexa(ind, ni, nj,
              xl, yl, zl,
              ic, jc, kc,
              xt, yt, zt);
    
    if (getCellJump(x, y, z,
                    xt, yt, zt,              
                    isomm,
                    xi, yi, zi) == false)
    {
      if (i == 7)
      {
        JUMP =  false;
        break;
      }
      /* Apply a technique of jump */
      is = E_min(8, E_Int((xi+E_sign(xi))*K_CONST::ONE_HALF));
      is = E_max(-8, is);
      js = E_min(8, E_Int((yi+E_sign(yi))*K_CONST::ONE_HALF));
      js = E_max(-8, js);
      ks = E_min(8, E_Int((zi+E_sign(zi))*K_CONST::ONE_HALF));
      ks = E_max(-8, ks);
      ind = ind+is+js*ni+ks*nij;

      kc = ind/nij;
      E_Int kcnij = kc*nij;
      jc = (ind-kcnij)/ni;
      E_Int jcni = jc*ni;
      ic = ind-jcni-kcnij;

      if (ic<0 || ic>ni-2 || jc<0 || jc>nj-2 || kc<0 || kc>nk-2)
      {
        JUMP =  false;
        break;
      }
    }
    else
    {        
      for (itr = listOfCandidateCells.begin();
           itr != listOfCandidateCells.end();
           itr++)
      {
        if (ind == *itr)
        {
          JUMP = true;
          goto saut;
        }
      }
      JUMP =  false;
      break;
    }
  }

  /* If the technique of jump succeeds, find the interpolation coefficients
     in the cell by cut it in 24 tetrahedras */
  saut:
  if (JUMP)
  {
    if (getCoeffInterpHexa(x, y, z,
                           isomm,
                           xi, yi, zi, 
                           xt, yt, zt,
                           cf) == true)
    {
      return 1;
    }
    else
    {
      JUMP =  false;
      listOfCandidateCells.erase(itr);
    }
  }
 
  /* If the technique of jump fails, we test all candidate cells */
  if (JUMP ==  false) 
  {    
    for (itr = listOfCandidateCells.begin();
         itr != listOfCandidateCells.end();
         itr++)
    {
      ind = (*itr);
          
      coordHexa(ind, ni, nj,
                xl, yl, zl,
                ic, jc, kc,
                xt, yt, zt);
      if (coeffInterpHexa(x, y, z,
                          xt, yt, zt,
                          cf) == true)
        return 1;
    }
  }
  return 0;
}
// ============================================================================
/* Find the extrapolation cell for point (x,y,z) for structured meshes */
// ============================================================================
short K_KINTERP::BlkInterpAdt::searchExtrapolationCell(
  E_Float x, E_Float y, E_Float z,
  E_Int& ic, E_Int& jc, E_Int& kc,
  FldArrayF& cf,
  E_Int order, E_Float cfMax,
  const FldArrayF& cellNatureField)
{
  // interpolation mesh
  //K_KINTERP::KMesh& mesh = _mesh;

  //E_Float* xl = mesh.getXVector();
  //E_Float* yl = mesh.getYVector();
  //E_Float* zl = mesh.getZVector();

  E_Int i,j,k;
  E_Int icsav = 0;
  E_Int jcsav = 0;
  E_Int kcsav = 0;
  E_Int foundInDomain = 0;
  E_Int ncf = 8;
  FldArrayF cf_sav(ncf);

  // Init cf
  cf.setAllValuesAtNull();

  // search list of candidate cells
  list<E_Int> listOfCandidateCells; 
  E_Int found = getListOfCandidateCells(x, y, z, listOfCandidateCells);
  
  if (found == 0) return 0; //listOfCandidateCells empty

  /* Find the right cell among candidate cells */
  list<E_Int>::iterator itr;
  E_Float sum_coef;
  E_Float saved_sum_coef = K_CONST::E_MAX_FLOAT;

  //const E_Float* cellNp = cellNatureField.begin(); // 0: blanked 1: other 
  E_Int cellNsize = cellNatureField.getSize();
  FldArrayI cellN(cellNsize);

  E_Int ind;
  // parameters for extrapolation routine
  E_Int is, js, ks, testNature;

  E_Int ni = _mesh.getIm();
  E_Int nj = _mesh.getJm();
  E_Int nij, kmnij, jmni;

  // Loop on all candidate cells : computation of coefficients cf
  for (itr = listOfCandidateCells.begin();
       itr != listOfCandidateCells.end();
       itr++)
  {
    // 1D-index of interpolation cell
    ind = *itr;
 
    // (i,j,k)-indices of interpolation cell
    nij = ni*nj;
    k = ind/nij+1; kmnij = (k-1)*nij;
    j = (ind-kmnij)/ni+1; jmni = (j-1)*ni;
    i = ind-jmni-kmnij+1;
    
    // Index of interpolation cell (in center)
    //E_Int ind0 = (i-1)+(j-1)*ni+(k-1)*ni*nj;
    //E_Int ind1 = ind0+1;
    //E_Int ind2 = ind0+ni;
    //E_Int ind3 = ind2+1;
    //E_Int ind4 = ind0+ni*nj;
    //E_Int ind5 = ind4+1;
    //E_Int ind6 = ind4+ni;
    //E_Int ind7 = ind6+1;
    
    // cellN for points of interpolation cell
    //E_Float cellNm0 = cellNp[ind0];
    //E_Float cellNm1 = cellNp[ind1];
    //E_Float cellNm2 = cellNp[ind2];
    //E_Float cellNm3 = cellNp[ind3];
    //E_Float cellNm4 = cellNp[ind4];
    //E_Float cellNm5 = cellNp[ind5];
    //E_Float cellNm6 = cellNp[ind6];
    //E_Float cellNm7 = cellNp[ind7];

    // compute interpolation coefficients for hexahedra
    // - the interpType is always O2CF in this routine
    K_KINTERP::BlkInterpData::InterpolationType interpType = K_KINTERP::BlkInterpData::O2CF;
    // - the interpolation mesh type is always NODE in this routine
    K_KINTERP::BlkInterpData::InterpMeshType interpMeshType = K_KINTERP::BlkInterpData::NODE;
    // - cellNatureField has to be an integer array
    for (E_Int ii = 0; ii < cellNsize; ii++)
      cellN[ii]=(E_Int)cellNatureField[ii];
    // - testNature=1: no blanked cells (criterion for the validity of an interpolation cell)
    testNature = 1;
    // - order: order of the extrapolation
    // (is, js, ks): neighbour cell (not used here)
    K_KINTERP::BlkInterpWithKMesh::getExtrapolationCoeffForCell(x, y, z, i, j, k, cf, cellN, testNature, order, 
                                                                is, js, ks, cfMax, interpType, interpMeshType);
    
    // Keep extrapolation cell with minimum distance (~ minimum sum of absolute coefficients) 
    // from the interpolated point
    sum_coef = E_abs(cf[0])+E_abs(cf[1])+E_abs(cf[2])+E_abs(cf[3])+E_abs(cf[4])+E_abs(cf[5])+E_abs(cf[6])+E_abs(cf[7]);
    if (sum_coef < saved_sum_coef && sum_coef < cfMax)
    {
      foundInDomain = 1;
      saved_sum_coef = sum_coef;
      icsav = i;
      jcsav = j;
      kcsav = k;
      cf_sav=cf;    
    }
  }

  ic = icsav; jc = jcsav; kc = kcsav;
  cf = cf_sav;
  return foundInDomain;
}
// ============================================================================
/* Find the extrapolation cell for point (x,y,z) for unstructured meshes */
// ============================================================================
short K_KINTERP::BlkInterpAdt::searchExtrapolationCell(
  E_Float x, E_Float y, E_Float z,
  E_Int& noelt,
  FldArrayF& cf,
  FldArrayI& indi,
  E_Int order,
  const FldArrayF& cellNatureField)
{
  // interpolation mesh
  K_KINTERP::KMesh& mesh = _mesh;

  E_Float* xt = mesh.getXVector();
  E_Float* yt = mesh.getYVector();
  E_Float* zt = mesh.getZVector();

  FldArrayI& connectEV = _mesh.getConnectivity();

  // datas for interpolation cell (tetrahedra)
  E_Float xp, yp, zp;
  E_Float xq, yq, zq;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Int indp, indq, indr, inds;
  E_Float xi, yi, zi;

  // threshold of sul of extrapolation coefs for order 1
  E_Float threshold = 1.e5;

  E_Int foundInDomain = 0;
  E_Int noeltsave = 0;
  E_Int ncf = 8;
  FldArrayF cf_sav(ncf);

  // Init cf
  cf.setAllValuesAtNull();

  // search list of candidate cells
  list<E_Int> listOfCandidateCells; 
  E_Int found = getListOfCandidateCells(x, y, z, listOfCandidateCells);
  
  if (found == 0) return 0; //listOfCandidateCells empty

  /* Find the right cell among candidate cells */
  list<E_Int>::iterator itr;
  E_Float sum;
  E_Float sum_coef;
  E_Float saved_sum_coef = K_CONST::E_MAX_FLOAT;
  E_Float max_coef;
  const E_Float* cellNp = cellNatureField.begin(); // 0: blanked 1: other  
  E_Int* cn1 = connectEV.begin(1);
  E_Int* cn2 = connectEV.begin(2);
  E_Int* cn3 = connectEV.begin(3);
  E_Int* cn4 = connectEV.begin(4);

  // Loop on all candidate cells : computation of coefficients cf
  for (itr = listOfCandidateCells.begin();
       itr != listOfCandidateCells.end();
       itr++)
  {
    // index of candidate cell
    noelt = *itr;
    
    // compute interpolation coefficients for tetrahedra
    indp = cn1[noelt]-1;
    indq = cn2[noelt]-1;
    indr = cn3[noelt]-1;
    inds = cn4[noelt]-1;
    
    xp = xt[indp]; yp = yt[indp]; zp = zt[indp];
    xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
    xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
    xs = xt[inds]; ys = yt[inds]; zs = zt[inds];
    coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                     xr, yr, zr, xs, ys, zs, xi, yi, zi);
    
    // cellN for points of interpolation cell
    E_Float cellNmp = cellNp[indp];
    E_Float cellNmq = cellNp[indq];
    E_Float cellNmr = cellNp[indr];
    E_Float cellNms = cellNp[inds];
    
    // extrapolation - order 0
    if (order == 0)
    {
      for (E_Int i = 0; i<cf.getSize();i++)
        cf[i]=0.;
      // test will point of the tetrahedra will be used for extrapolation
      max_coef = E_max(xi,yi); max_coef = E_max(max_coef,zi); 
      if (max_coef <= 0.5) 
        cf[0]=1.; // point P
      else if (K_FUNC::fEqualZero(max_coef,xi))
        cf[1]=1.; // point Q
      else if (K_FUNC::fEqualZero(max_coef,yi))
        cf[2]=1.; // point R
      else
        cf[3]=1.; // point S         
    }
    else if (order == 1) // extrapolation - order 1
    {
      cf[0] = 1-xi-yi-zi ; 
      cf[1] = xi;
      cf[2] = yi;
      cf[3] = zi;
      sum_coef = E_abs(cf[1])+E_abs(cf[2])+E_abs(cf[3]);
      // if sum_coef exceeds threshold, degenerate to order 0
      if (sum_coef >= threshold)
      {
        for (E_Int i = 0; i<cf.getSize();i++)
          cf[i]=0.;
        // test will point of the tetrahedra will be used for extrapolation
        max_coef = E_max(xi,yi); max_coef = E_max(max_coef,zi); 
        if (max_coef <= 0.5) 
          cf[0]=1.; // point P
        else if (K_FUNC::fEqualZero(max_coef,xi))
          cf[1]=1.; // point Q
        else if (K_FUNC::fEqualZero(max_coef,yi))
          cf[2]=1.; // point R
        else
          cf[3]=1.; // point S         
      }
    }
    // Keep only valid extrapolation cell (no blanked cell in tetrahedra)
    sum = cellNmp*cellNmq*cellNmr*cellNms;
    if (sum > K_CONST::E_CUTOFF)
    {
      // Keep extrapolation cell with minimum distance (~ minimum sum of absolute coefficients) 
      // from the interpolated point
      sum_coef = E_abs(cf[1])+E_abs(cf[2])+E_abs(cf[3]);
      if (sum_coef < saved_sum_coef)
      {
        foundInDomain = 1;
        saved_sum_coef = sum_coef;
        noeltsave = noelt;
        cf_sav=cf;    
      }
    }
  }

  noelt = noeltsave;
  cf = cf_sav;
  return foundInDomain;
}
// ============================================================================
/* 
   Find the interpolation cell for point (x,y,z) for a tetra kmesh
   le numero de l'elt trouve est retourne dans noelt
*/
// ============================================================================
short K_KINTERP::BlkInterpAdt::searchInterpolationCell(
  E_Float x, E_Float y, E_Float z,
  E_Int& noelt,
  FldArrayF& cf)
{ 
  K_KINTERP::KMesh& mesh = const_cast<K_KINTERP::KMesh&>(_mesh);

  cf.setAllValuesAtNull();

  // recherche de la liste des cellules candidates
  list<E_Int> listOfCandidateCells; 
  E_Int found = getListOfCandidateCells(x, y, z, listOfCandidateCells);

  if (found == 0) return 0; //listOfCandidateCells vide

  const E_Float EPS = _EPS_TETRA;

  E_Float* xt = mesh.getXVector();
  E_Float* yt = mesh.getYVector();
  E_Float* zt = mesh.getZVector();

  E_Float xp, yp, zp;
  E_Float xq, yq, zq;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Int indp, indq, indr, inds;
  E_Float xi, yi, zi, sum;
  E_Int et;

  list<E_Int>::iterator itr;
  FldArrayI& connectEV = mesh.getConnectivity();
  E_Int* cn1 = connectEV.begin(1);
  E_Int* cn2 = connectEV.begin(2);
  E_Int* cn3 = connectEV.begin(3);
  E_Int* cn4 = connectEV.begin(4);

  for (itr = listOfCandidateCells.begin(); 
       itr != listOfCandidateCells.end();
       itr++)
  {
    et = *itr;
    indp = cn1[et]-1;
    indq = cn2[et]-1;
    indr = cn3[et]-1;
    inds = cn4[et]-1;

    xp = xt[indp];
    yp = yt[indp];
    zp = zt[indp];
    
    xq = xt[indq];
    yq = yt[indq];
    zq = zt[indq];
    
    xr = xt[indr];
    yr = yt[indr];
    zr = zt[indr];
    
    xs = xt[inds];
    ys = yt[inds];
    zs = zt[inds];

    coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                     xr, yr, zr, xs, ys, zs, xi, yi, zi);

    sum = xi+yi+zi;
    if (xi > -EPS &&  yi > -EPS &&  zi > -EPS && sum < K_CONST::ONE+3*EPS)
    {
      cf[0] = 1-sum; 
      cf[1] = xi;
      cf[2] = yi;
      cf[3] = zi;
      noelt = et;
      return 1;
    }
  }
  noelt = -1;
  return 0;
}

//=============================================================================
short K_KINTERP::BlkInterpAdt::searchInterpolationCellByJump(
  E_Float x, E_Float y, E_Float z,
  E_Int& ic, E_Int& jc, E_Int& kc,
  FldArrayF& cf)
{
  return K_KINTERP::BlkInterpWithKMesh::
    searchInterpolationCellByJump(x, y, z, ic, jc, kc, cf);
}
 
//=============================================================================
void K_KINTERP::BlkInterpAdt::searchInterpolationCellByJumpv(
  FldArrayF& coord,E_Int istart, E_Int iend,
  FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
  FldArrayF& cf, FldArrayIS& found )
{
   K_KINTERP::BlkInterpWithKMesh::
     searchInterpolationCellByJumpv(coord, istart, iend,
                                    ic, jc, kc,cf, found);
}

// ============================================================================
/**/
// ============================================================================
void K_KINTERP::BlkInterpAdt::searchInterpolationCellv(
  FldArrayF& coord, E_Int istart, E_Int iend,
  FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
  FldArrayF& cf, FldArrayIS& found)
{
  FldArrayF cfloc(8); cfloc.setAllValuesAtNull();
  
  for (E_Int i = istart; i < iend; i++)
  {
    found[i] = searchInterpolationCell(coord(i,1), coord(i,2), coord(i,3), 
                                       ic[i], jc[i], kc[i],cfloc);
    for (E_Int j = 0; j < 8; j++) cf(i,j+1) = cfloc[j];
  }
}
//=============================================================================
/* Search the extrapolation cell on the list of meshes (the original one and the
   duplicated one(s))
   It starts from the current mesh (the original one if the previous
   interpolation cell was found on it) */
//=============================================================================
short 
K_KINTERP::BlkInterpAdt::getExtrapolationCell(E_Float x, E_Float y, E_Float z,
                                              E_Int& ic, E_Int& jc, E_Int& kc,
                                              FldArrayF& cf,
                                              const FldArrayI& cellNatureField,
                                              E_Int testNature,
                                              E_Float& test,
                                              K_KINTERP::BlkInterpData::InterpolationType interpType,
                                              K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
                                              E_Float cfMax)
{  
  return  K_KINTERP::BlkInterpWithKMesh::getExtrapolationCell( 
    x, y, z, ic, jc, kc,
    cf, cellNatureField, testNature, test, interpType, interpMeshType, cfMax);
}

// ============================================================================
short K_KINTERP::BlkInterpAdt::getExtrapolationCellStruct(E_Float x, E_Float y, E_Float z,
                                                         FldArrayI& indi,
                                                         FldArrayF& cf,
                                                         E_Int order, E_Float cfMax,
                                                         const FldArrayF& cellNatureField,
                                                         InterpolationType interpType)
{
  return  K_KINTERP::BlkInterpWithKMesh::getExtrapolationCellStruct( 
    x, y, z, indi, cf, order, cfMax, cellNatureField, interpType);
}

// ============================================================================
short K_KINTERP::BlkInterpAdt::getExtrapolationCellUnstr(E_Float x, E_Float y, E_Float z,
                                                        E_Int& noelt,
                                                        FldArrayF& cf,
                                                        FldArrayI& indi,
                                                        E_Int order,
                                                        const FldArrayF& cellNatureField)
{
  return  K_KINTERP::BlkInterpWithKMesh::getExtrapolationCellUnstr( 
    x, y, z, noelt, cf, indi, order, cellNatureField);
}

//=============================================================================
/* Compute the extrapolation coefficients for point (x,y,z) inside cell
   (ic,jc,kc) (in extended centers). It is based on the coefficients in the
   best tetrahedra. In this version (default), the mesh in centers is supposed
   to be known.*/
//=============================================================================
short K_KINTERP::BlkInterpAdt::getExtrapolationCoeffForCell(
  E_Float x, E_Float y, E_Float z,
  E_Int ic, E_Int jc, E_Int kc,
  FldArrayF& cf, E_Float cfMax,
  K_KINTERP::BlkInterpData::InterpolationType interpType,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType)
{
  return  K_KINTERP::BlkInterpWithKMesh::getExtrapolationCoeffForCell( 
    x, y, z, ic, jc, kc, cf, cfMax, interpType, interpMeshType);
}

//=============================================================================
short K_KINTERP::BlkInterpAdt::getExtrapolationCoeffForCell(
  E_Float x, E_Float y, E_Float z,
  E_Int ic, E_Int jc, E_Int kc,
  FldArrayF& cf,
  const FldArrayI& cellNatureField,
  E_Int testNature, E_Int order,
  E_Int& is, E_Int& js, E_Int& ks, E_Float cfMax,
  K_KINTERP::BlkInterpData::InterpolationType interpType,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType)
{
  return  K_KINTERP::BlkInterpWithKMesh::getExtrapolationCoeffForCell( 
    x, y, z, ic, jc, kc, cf,
    cellNatureField, testNature, order,
    is, js, ks, cfMax, interpType, interpMeshType);
}
// =====================Interp/BlkInterpAdt_getCell.cpp ==================
