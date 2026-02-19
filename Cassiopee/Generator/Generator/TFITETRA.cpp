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

// TFI generator - TETRA

# include "generator.h"

using namespace std;
using namespace K_FLD;

//===========================================================================
/* TFI TETRA */
//===========================================================================
PyObject* K_GENERATOR::TFITETRA(PyObject* arrays)
{
  // Extract infos from arrays
  vector<E_Int> res;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; 
  vector<E_Int> njt; 
  vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false; 
  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int ns0 = structF.size(); E_Int nu0 = unstrF.size();
  //Verification de la nature des arrays
  if ( unstrF.size() != 4 || isOk != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "TFITETRA: only 4 TRI arrays must be defined.");
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }
  // Verification des TRI
  for (E_Int v = 0; v < 4; v++)
  {
    if ( strcmp(eltType[v], "TRI") != 0 ) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "TFITETRA: unstructured arrays must be TRI-type.");    
      for (E_Int nos = 0; nos < 4; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
  }
  
  // verification de la taille des triangulations
  E_Int nelts= cnt[0]->getSize();
  E_Int nptsTRI = unstrF[0]->getSize();
  for (E_Int i = 1; i < 4; i++)
  {
    if (cnt[i]->getSize() != nelts) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "TFITETRA: all the TRI arrays must have the same number of elements.");    
      for (E_Int nos = 0; nos < 4; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
    if ( unstrF[i]->getSize() != nptsTRI)
    {
      PyErr_SetString(PyExc_TypeError,
                      "TFITETRA: all the TRI arrays must have the same number of points.");    
      for (E_Int nos = 0; nos < 4; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
  }
  vector<E_Int> posxt(4); vector<E_Int> posyt(4); vector<E_Int> poszt(4);
  for (E_Int i = 0; i < 4; i++)//coordinates already checked in getFromArrays
  {
    posxt[i] = K_ARRAY::isCoordinateXPresent(unstrVarString[i])+1;
    posyt[i] = K_ARRAY::isCoordinateYPresent(unstrVarString[i])+1;
    poszt[i] = K_ARRAY::isCoordinateZPresent(unstrVarString[i])+1;
  }
  E_Float* xt1 = unstrF[0]->begin(posxt[0]);
  E_Float* yt1 = unstrF[0]->begin(posyt[0]);
  E_Float* zt1 = unstrF[0]->begin(poszt[0]);
  
  E_Float* xt2 = unstrF[1]->begin(posxt[1]);
  E_Float* yt2 = unstrF[1]->begin(posyt[1]);
  E_Float* zt2 = unstrF[1]->begin(poszt[1]);

  E_Float* xt3 = unstrF[2]->begin(posxt[2]);
  E_Float* yt3 = unstrF[2]->begin(posyt[2]);
  E_Float* zt3 = unstrF[2]->begin(poszt[2]);

  E_Float* xt4 = unstrF[3]->begin(posxt[3]);
  E_Float* yt4 = unstrF[3]->begin(posyt[3]);
  E_Float* zt4 = unstrF[3]->begin(poszt[3]);

  E_Int n = nptsTRI;// nb de points : n = ni*(ni+1)/2
  E_Float delta= sqrt(1+8*n);
  E_Int ni = (E_Int)((-1.+delta)/2.);// nb de points sur une arete de triangle
  E_Float npts0 = ni*(ni*ni/6.+ni/2.+1./3);
  E_Int npts = E_Int(npts0+1.e-8);

  short isok = checkContinuousBndTETRA(ni, xt1, yt1, zt1, xt2, yt2, zt2, xt3, yt3, zt3, xt4, yt4, zt4);
  if ( isok == 0 )
  {
    for (E_Int nos = 0; nos < 4; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

    PyErr_SetString(PyExc_TypeError,
                    "TFI: boundaries do not match.");
    return NULL;
  }

  //connectivite
  E_Int ni1 = ni-1;
  E_Int ntetra = ni1*ni1*ni1;
  char* varString = unstrVarString[0];
  E_Int api = unstrF[0]->getApi();
  PyObject* tpl = K_ARRAY::buildArray3(3, varString, npts, ntetra, "TETRA", false, api);
  FldArrayF* coord; FldArrayI* cn;
  K_ARRAY::getFromArray3(tpl, coord, cn);
  cn->setAllValuesAt(1);
 
  E_Float* xt = coord->begin(1);
  E_Float* yt = coord->begin(2);
  E_Float* zt = coord->begin(3);

  E_Int ind = 0;
  E_Int invni1 = E_Int(1./E_Float(ni1));
  E_Float l1, l2, l3, l4;
  
  for (E_Int k = 0; k < ni; k++)
    for (E_Int j = 0; j < ni-k; j++)
      for (E_Int i = 0; i < ni-j-k; i++)
      {
        l2 = i*invni1;
        l3 = j*invni1;
        l4 = k*invni1;
        l1 = 1.-l2-l3-l4; 
        // Triangle t1: P1P2,P1P3 (lambda4=0)
        // Triangle t2: P2P3,P2P4 (lambda1=0)
        // Triangle t3: P3P4,P3P1 (lambda2=0)
        // Triangle t4: P4P1,P4P2 (lambda3=0)

        xt[ind] = 
          l1*xt1[i+j*ni-j*(j-1)/2] + l2*xt1[i+k+j*ni-j*(j-1)/2] + l3*xt1[i+(j+k)*ni-(j+k)*(j+k-1)/2]+

          l2*xt2[j+k*ni-k*(k-1)/2] + l3*xt2[ni1-i-k+k*ni-k*(k-1)/2] + l4*xt2[j+(ni1-i-j)*ni-(ni1-i-j)*(ni1-i-j-1)/2]+ 

          l3*xt3[k+(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2]+l4*xt3[i+k+(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2] + l1*xt3[k+(ni1-j-k)*ni-(ni1-j-k)*(ni1-j-k-1)/2] +

          l4*xt4[ni1-i-j-k+i*ni-i*(i-1)/2] + l1*xt4[ni1-i-k+i*ni-i*(i-1)/2] + l2*xt4[ni1-i-j-k+(i+j)*ni-(i+j)*(i+j-1)/2] 

          -l1*xt1[i] - l2*xt1[i+j+k] - l2*xt2[j] - l3*xt2[ni1-i]

          -l3*xt3[(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2] - l1*xt3[(ni1-j)*ni-(ni1-j)*(ni1-j-1)/2] - l1*xt4[ni1-k] -  l4*xt3[i+j+k + (ni1-i-j-k)*ni - (ni1-i-j-k)*(ni1-i-j-k-1)/2]

          -l2*xt2[k*ni-k*(k-1)/2]-l4*xt2[(ni1-i)*ni-(ni1-i)*(ni1-i-1)/2]-l3*xt2[ni1-k+k*ni-k*(k-1)/2] -l4*xt2[j+(ni1-j)*ni-(ni1-j)*(ni1-j-1)/2]


          + l1 * xt1[0] + l2 * xt2[0] + l3*xt3[0] +l4*xt4[0];


        yt[ind] = 
          l1*yt1[i+j*ni-j*(j-1)/2] + l2*yt1[i+k+j*ni-j*(j-1)/2] + l3*yt1[i+(j+k)*ni-(j+k)*(j+k-1)/2]+

          l2*yt2[j+k*ni-k*(k-1)/2] + l3*yt2[ni1-i-k+k*ni-k*(k-1)/2] + l4*yt2[j+(ni1-i-j)*ni-(ni1-i-j)*(ni1-i-j-1)/2]+ 

          l3*yt3[k+(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2]+l4*yt3[i+k+(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2] + l1*yt3[k+(ni1-j-k)*ni-(ni1-j-k)*(ni1-j-k-1)/2] +

          l4*yt4[ni1-i-j-k+i*ni-i*(i-1)/2] + l1*yt4[ni1-i-k+i*ni-i*(i-1)/2] + l2*yt4[ni1-i-j-k+(i+j)*ni-(i+j)*(i+j-1)/2] 

          -l1*yt1[i] - l2*yt1[i+j+k] - l2*yt2[j] - l3*yt2[ni1-i]

          -l3*yt3[(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2] - l1*yt3[(ni1-j)*ni-(ni1-j)*(ni1-j-1)/2] - l1*yt4[ni1-k] -  l4*yt3[i+j+k + (ni1-i-j-k)*ni - (ni1-i-j-k)*(ni1-i-j-k-1)/2]

          -l2*yt2[k*ni-k*(k-1)/2]-l4*yt2[(ni1-i)*ni-(ni1-i)*(ni1-i-1)/2]-l3*yt2[ni1-k+k*ni-k*(k-1)/2] -l4*yt2[j+(ni1-j)*ni-(ni1-j)*(ni1-j-1)/2]


          + l1 * yt1[0] + l2 * yt2[0] + l3*yt3[0] +l4*yt4[0];


        zt[ind] = 
          l1*zt1[i+j*ni-j*(j-1)/2] + l2*zt1[i+k+j*ni-j*(j-1)/2] + l3*zt1[i+(j+k)*ni-(j+k)*(j+k-1)/2]+

          l2*zt2[j+k*ni-k*(k-1)/2] + l3*zt2[ni1-i-k+k*ni-k*(k-1)/2] + l4*zt2[j+(ni1-i-j)*ni-(ni1-i-j)*(ni1-i-j-1)/2]+ 

          l3*zt3[k+(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2]+l4*zt3[i+k+(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2] + l1*zt3[k+(ni1-j-k)*ni-(ni1-j-k)*(ni1-j-k-1)/2] +

          l4*zt4[ni1-i-j-k+i*ni-i*(i-1)/2] + l1*zt4[ni1-i-k+i*ni-i*(i-1)/2] + l2*zt4[ni1-i-j-k+(i+j)*ni-(i+j)*(i+j-1)/2] 

          -l1*zt1[i] - l2*zt1[i+j+k] - l2*zt2[j] - l3*zt2[ni1-i]

          -l3*zt3[(ni1-i-j-k)*ni-(ni1-i-j-k)*(ni1-i-j-k-1)/2] - l1*zt3[(ni1-j)*ni-(ni1-j)*(ni1-j-1)/2] - l1*zt4[ni1-k] -  l4*zt3[i+j+k + (ni1-i-j-k)*ni - (ni1-i-j-k)*(ni1-i-j-k-1)/2]

          -l2*zt2[k*ni-k*(k-1)/2]-l4*zt2[(ni1-i)*ni-(ni1-i)*(ni1-i-1)/2]-l3*zt2[ni1-k+k*ni-k*(k-1)/2] -l4*zt2[j+(ni1-j)*ni-(ni1-j)*(ni1-j-1)/2]


          + l1 * zt1[0] + l2 * zt2[0] + l3*zt3[0] +l4*zt4[0];

        ind++;
      }

  //connectivity
  E_Int et = 0;
  E_Int Nj, Nk;
  E_Int offi=0, offj=0;
  E_Int offi2=0, offj2=0;
  
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);
  E_Int* cn4 = cn->begin(4);

  for (E_Int k = 0; k < ni1; k++)
  { 
    Nk = (ni-k+1)*(ni-k)/2;
    for (E_Int j = 0; j < ni1-k; j++)
    {
      offj = Nk-j;
      offj2 = Nk-(j+1);
      Nj = ni-j-k;
      for (E_Int i = 0; i < ni1-j-k-1; i++)
      { 
        cn1[et] = offi + 1; 
        cn2[et] = offi + 2;
        cn3[et] = offi + Nj + 1;
        cn4[et] = offi + offj + 1;
        et++;     

        cn1[et] = offi + 2;
        cn2[et] = offi + Nj + 2;
        cn3[et] = offi + Nj + 1;
        cn4[et] = offi + offj + 1;
        et++;

        cn1[et] = offi + 2;
        cn2[et] = offi + Nj + 2;
        cn3[et] = offi + offj + 1;
        cn4[et] = offi + offj + 2;
        et++;

        // deuxieme bande
        offi2 =  offi + Nj;
        cn1[et] = offi2 + 1; 
        cn2[et] = offi2 + 2; 
        cn3[et] = offi2 + offj2 + 1;  
        cn4[et] = offi + offj + 1; 
        et++;
        
        cn1[et] = offi2 + 2;
        cn2[et] = offi + offj + 1;
        cn3[et] = offi + offj + 2;        
        cn4[et] = offi2 + offj2 + 1;  
        et++;
          
        if ( i <  ni1-j-k-2)
        {
          cn1[et] = offi2 + 2;
          cn2[et] = offi + offj + 2;
          cn3[et] = offi2 + offj2 + 2;        
          cn4[et] = offi2 + offj2 + 1;  
          et++;
        }
        // fin deuxieme bande
        offi++;
      }


      cn1[et] = offi + 1; 
      cn2[et] = offi + 2;
      cn3[et] = offi + Nj + 1;
      cn4[et] = offi + offj + 1;
      et++;
    
      offi++;
      // pour aller au j suivant
      offi++;   
    }
    offi++;
  }

  RELEASESHAREDU(tpl, coord, cn);
  for (E_Int nos = 0; nos < 4; nos++)
    RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
  return tpl;
}
//=============================================================================
/**/
//=============================================================================
short K_GENERATOR::checkContinuousBndTETRA(E_Int ni,
                                           E_Float* xt1, E_Float* yt1, E_Float* zt1,
                                           E_Float* xt2, E_Float* yt2, E_Float* zt2,
                                           E_Float* xt3, E_Float* yt3, E_Float* zt3,
                                           E_Float* xt4, E_Float* yt4, E_Float* zt4)
{
  E_Float eps = 1.e-6;
  E_Float dx, dy, dz;
  E_Int ni1 = ni-1;
  // Triangle t1: P1P2,P1P3 (lambda4=0)
  // Triangle t2: P2P3,P2P4 (lambda1=0)
  // Triangle t3: P3P4,P3P1 (lambda2=0)
  // Triangle t4: P4P1,P4P2 (lambda3=0)
  // verification des sommets des 4 triangles
  // Sommet P1
  E_Int indtmax = ni*(ni+1)/2-1; 
  dx = xt1[0] - xt3[indtmax];
  dy = yt1[0] - yt3[indtmax];
  dz = zt1[0] - zt3[indtmax];
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
  {
    printf("TFITETRA: P1(tri1) must match P3(tri3)\n");
    return 0;
  }
  dx = xt1[0]-xt4[ni1];
  dy = yt1[0]-yt4[ni1];
  dz = zt1[0]-zt4[ni1];
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
  {
    printf("TFITETRA: P1(tri1) must match P2(tri4)\n");
    return 0;
  }
  //Sommet P2
  //triangles t1 et t2
  dx = xt2[0] - xt1[ni1];
  dy = yt2[0] - yt1[ni1];
  dz = zt2[0] - zt1[ni1];
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
  {
    printf("TFITETRA: P1(tri2) must match P2(tri1)\n");
    return 0;
  }
  //triangles t1 et t4
  dx = xt2[0] - xt4[indtmax];
  dy = yt2[0] - yt4[indtmax];
  dz = zt2[0] - zt4[indtmax];
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
  {
    printf("TFITETRA: P1(tri2) must match P3(tri4)\n");
    return 0;
  }
  // Sommet P3
  //triangle t3 et t1
  dx = xt3[0] - xt1[indtmax];
  dy = yt3[0] - yt1[indtmax];
  dz = zt3[0] - zt1[indtmax];
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
  {
    printf("TFITETRA: P1(tri3) must match P3(tri1)\n");
    return 0;
  }
  //triangles t3 et t2
  dx = xt3[0] - xt2[ni1];
  dy = yt3[0] - yt2[ni1];
  dz = zt3[0] - zt2[ni1];
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
  {
    printf("TFITETRA: P1(tri3) must match P2(tri2)\n");
    return 0;
  }
  // Sommet P4
  //triangles t4 et t2
  dx = xt4[0] - xt2[indtmax];
  dy = yt4[0] - yt2[indtmax];
  dz = zt4[0] - zt2[indtmax];
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
  {
    printf("TFITETRA: P1(tri4) must match P3(tri2)\n");
    return 0;
  }
  //triangles t4 et t3
  dx = xt4[0] - xt3[ni1];
  dy = yt4[0] - yt3[ni1];
  dz = zt4[0] - zt3[ni1];
  if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
  {
    printf("TFITETRA: P1(tri4) must match P2(tri3)\n");
    return 0;
  }
  return 1;
}
//=============================================================================
/* reorder des faces TRI TFI */
//=============================================================================
// void K_GENERATOR::reorderTFITETRA(E_Int ni, E_Int nptsTRI,
//                                   vector<E_Int>& posxt, vector<E_Int>& posyt, vector<E_Int>& poszt,
//                                   vector<FldArrayF*> fields, FldArrayIS& newOrder, E_Float eps)
// {
  // newOrder.setAllValuesAt(-1);
  // newOrder[0] = 0;//triangle t1 de reference
  // E_Int nfld = fields[0]->getNfld();
  // // triangle t1: faces B=P2(t1) et C=P3(t1) coincidentes avec un autre triangle
  // E_Int indBt1 = ni-1; E_Int indCt1 = ni*(ni+1)/2-1;
  // E_Int indBtopp=-1; E_Int indCtopp = -1; notopp = -1;// triangle t2 a trouver
  
  // E_Int posxt1 = posxt[0]; E_Int posyt1 = posyt[0]; E_Int poszt1 = poszt[0];  
  // E_Float* xt1 = fields[0]->begin(posxt1);
  // E_Float* yt1 = fields[0]->begin(posyt1);
  // E_Float* zt1 = fields[0]->begin(poszt1);
  // E_Int found = 0;
  // for (E_Int v = 1; v < 4; v++)
  // {
  //   E_Int posxtopp = posxt[v]; 
  //   E_Int posytopp = posyt[v]; 
  //   E_Int posztopp = poszt[v]; 
  //   E_Float* xtopp = fields[v]->begin(posxtopp);
  //   E_Float* ytopp = fields[v]->begin(posytopp);
  //   E_Float* ztopp = fields[v]->begin(posztopp);
  //   found = 0;
  //   for (E_Int indopp = 0; indopp < nptsTRI; indopp++)
  //   {
  //     if (K_FUNC::fEqualZero(xtopp[indopp]-xt1[indBt1],eps) == true ) 
  //     {
  //       indBtopp = indopp; found += 1;
  //     }
  //     if (K_FUNC::fEqualZero(xtopp[indopp]-xt1[indCt1],eps) == true ) 
  //     {
  //       indCtopp = indopp; found += 1;
  //     }
  //     if ( found == 2 ) {notopp = v; goto tri2;}
  //   }
  // }
  // tri2:;
  // //reorder eventuel ?
  // newOrder[1] = notopp;
  // FldArrayF coordsT2(nptsTRI,nfld);
  // if ( indBtopp == 0 && indCtopp == nptsTRI-1) coordsT2 = *fields[notopp];
  // // Reinit pour t3
// }
