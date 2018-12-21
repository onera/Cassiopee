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

// TFI generator - PENTA

# include "generator.h"

using namespace std;
using namespace K_FLD;

//===========================================================================
/* TFI PENTA */
//===========================================================================
PyObject* K_GENERATOR::TFIPENTA(PyObject* arrays)
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
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false; 
  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  E_Int ns0 = structF.size(); E_Int nu0 = unstrF.size();
  //Verification de la nature des arrays
  if ( unstrF.size() != 2 || structF.size() != 3 || isOk != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "TFI: only 2 TRI and 3 2D-structured arrays must be defined.");
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

    return NULL;
  }
  // verification des TRI
  E_Int nunstr = unstrF.size();
  E_Int nstruct = structF.size();

  for (E_Int v = 0; v < nunstr; v++)
  {
    if ( strcmp(eltType[v], "TRI") != 0 ) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "TFI: unstructured arrays must be TRI-type.");
      for (E_Int nos = 0; nos < ns0; nos++)
        RELEASESHAREDS(objs[nos], structF[nos]);
      for (E_Int nos = 0; nos < nu0; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
  }
  //verification des 2D struct
  for (E_Int v = 0; v < nstruct; v++)
  {
    E_Int ni = nit[v]; E_Int nj = njt[v]; E_Int nk = nkt[v];
    if ( ni < 2 || nj < 2 || nk != 1 )
    {
      for (E_Int nos = 0; nos < ns0; nos++)
        RELEASESHAREDS(objs[nos], structF[nos]);
      for (E_Int nos = 0; nos < nu0; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      PyErr_SetString(PyExc_TypeError,
                    "TFI: structured arrays must be 2D.");
      return NULL;
    }
  } 
  // coordonnees deja verifiees dans getFromArrays
  char* varString = structVarString[0];
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); 
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); 
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;
  
  FldArrayF* F1 = unstrF[0];
  FldArrayF* F2 = unstrF[1];
  FldArrayF* F3 = structF[0];
  FldArrayF* F4 = structF[1];
  FldArrayF* F5 = structF[2];
  FldArrayI* cn1 = cnt[0]; FldArrayI* cn2 = cnt[1]; 

  // determination du nb de noeuds sur une arete
  E_Int n11 = E_Int(sqrt(cn1->getSize()*1.)) + 1;
  if ( (n11 * n11 + n11)/2 != F1->getSize() ) 
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

    PyErr_SetString(PyExc_TypeError,
                    "TFI: first triangle must be a TFI triangle.");
    return NULL;
  }
  E_Int n21 = E_Int(sqrt(cn2->getSize()*1.)) + 1;
  if ( (n21 * n21 + n21)/2 != F2->getSize() ) 
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

    PyErr_SetString(PyExc_TypeError,
                    "TFI: second triangle must be a TFI triangle.");
    return NULL;
  }
  if ( n11 != n21 )
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    PyErr_SetString(PyExc_TypeError,
                    "TFI: triangles must have the same dimensions.");
    return NULL;
  }
  E_Int n = n11; // nb de sommets sur chaque cote des triangles
  E_Int im3 = nit[0]; E_Int jm3 = njt[0];
  E_Int im4 = nit[1]; E_Int jm4 = njt[1];
  E_Int im5 = nit[2]; E_Int jm5 = njt[2];

  // verification du nb de pts pour les quad
  if ( im3 != n || im4 != n || im5 != n || 
       jm4 != jm3 || jm4 != jm5 || jm3 != jm5 )
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

    PyErr_SetString(PyExc_TypeError,
                    "TFI: quads must have the same j dimensions and the same i dimensions as triangles.");
    return NULL; 
  }  
  E_Int p = jm3;

  E_Float* xtmin = F1->begin(posx);
  E_Float* ytmin = F1->begin(posy);
  E_Float* ztmin = F1->begin(posz);

  E_Float* xtmax = F2->begin(posx);
  E_Float* ytmax = F2->begin(posy);
  E_Float* ztmax = F2->begin(posz);

  E_Float* ximin = F3->begin(posx);
  E_Float* yimin = F3->begin(posy);
  E_Float* zimin = F3->begin(posz);
  
  E_Float* xjmin = F4->begin(posx);
  E_Float* yjmin = F4->begin(posy);
  E_Float* zjmin = F4->begin(posz);
  
  E_Float* xdiag = F5->begin(posx);
  E_Float* ydiag = F5->begin(posy);
  E_Float* zdiag = F5->begin(posz);

  short isok = checkContinuousBndsPENTA(
    n, p, xtmin, ytmin, ztmin,
    xtmax, ytmax, ztmax,
    ximin, yimin, zimin,
    xjmin, yjmin, zjmin,
    xdiag, ydiag, zdiag);
  if ( isok == 0 )
  {
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objs[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

    PyErr_SetString(PyExc_TypeError,
                    "TFI: boundaries don't match.");
    return NULL;
  }

  E_Int npts = (n*(n+1)/2)*p;
  E_Int ind = 0;
  E_Int n1 = n-1;
  E_Int invn1 = 1/n1;
  E_Int ip, jp, kp, ip1, kp1;
  E_Int p1 = p-1;
  E_Int invp1 = 1/p1;
  E_Int indt, inddiagj, inddiagi, indik, indjk;
  E_Int nelts = n1*n1*p1; 
  PyObject* tpl = K_ARRAY::buildArray(3, varString, npts, nelts, -1, "PENTA");
  E_Float* coordp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF coord(npts, 3, coordp, true);
  E_Int* cnp = K_ARRAY::getConnectPtr(tpl);
  FldArrayI cn(nelts,6, cnp, true);
  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);

  for (E_Int k = 0; k < p; k++)
  {
    for (E_Int j = 0; j < n; j++)
      for (E_Int i = 0; i < n-j; i++)
      {
        ip = i*invn1;
        jp = j*invn1;
        kp = k*invp1;
        ip1 = 1-ip-jp;
        kp1 = 1-kp;
        indt = (2*n*j-j*j+j)/2+i;
        indik = i + n*k;
        indjk = j + n*k;
        inddiagj = (2*n*j-j*j+j)/2+n-j-1;
        inddiagi = (2*n*(n1-i)-(n1-i)*(n1-i)+(n1-i))/2+n-(n1-i)-1;

        xt[ind] = 
          kp1 * xtmin[indt] + kp * xtmax[indt] + 
          ip1* xjmin[indik] + ip * xjmin[i+j+n*k] + 
          ip * xdiag[indjk] + jp * xdiag[n1-i+n*k] + jp * ximin[i+j + n*k] +
          ip1 * ximin[indjk] 
          - ip1* ( kp1*xtmin[i] + kp * xtmax[i] ) 
          - ip * ( kp1*xtmin[i+j] + kp * xtmax[i+j]) 
          - ip * ( kp1*xtmin[inddiagj] + kp * xtmax[inddiagj])
          - jp * ( kp1*xtmin[inddiagi] + kp * xtmax[inddiagi] )
          - jp * ( kp1*xtmin[(2*n*(i+j)-(i+j)*(i+j)+(i+j))/2] + kp * xtmax[(2*n*(i+j)-(i+j)*(i+j)+(i+j))/2] ) 
          - ip1*( kp1* xtmin[(2*n*j-j*j+j)/2] + kp * xtmax[(2*n*j-j*j+j)/2] ) 
          - ip1* xjmin[k*n] - ip * xdiag[k*n] - jp * ximin[n1+k*n]
          + ip1*( kp1 * xjmin[0] + kp * xjmin[p1*n]) 
          +       ip *( kp1 * xdiag[0] + kp * xdiag[p1*n]) 
          +       jp *( kp1 * ximin[n1] + kp * ximin[n1+p1*n]); 

        yt[ind] = 
          kp1 * ytmin[indt] + kp * ytmax[indt] + 
          ip1* yjmin[indik] + ip * yjmin[i+j+n*k] + 
          ip * ydiag[indjk] + jp * ydiag[n1-i+n*k] + jp * yimin[i+j + n*k] +
          ip1 * yimin[indjk] 
          - ip1* ( kp1*ytmin[i] + kp * ytmax[i] ) 
          - ip * ( kp1*ytmin[i+j] + kp * ytmax[i+j]) 
          - ip * ( kp1*ytmin[inddiagj] + kp * ytmax[inddiagj])
          - jp * ( kp1*ytmin[inddiagi] + kp * ytmax[inddiagi] )
          - jp * ( kp1*ytmin[(2*n*(i+j)-(i+j)*(i+j)+(i+j))/2] + kp * ytmax[(2*n*(i+j)-(i+j)*(i+j)+(i+j))/2] ) 
          - ip1*( kp1* ytmin[(2*n*j-j*j+j)/2] + kp * ytmax[(2*n*j-j*j+j)/2] ) 
          - ip1* yjmin[k*n] - ip * ydiag[k*n] - jp * yimin[n1+k*n]
          + ip1*( kp1 * yjmin[0] + kp * yjmin[p1*n]) 
          +       ip *( kp1 * ydiag[0] + kp * ydiag[p1*n]) 
          +       jp *( kp1 * yimin[n1] + kp * yimin[n1+p1*n]);

        zt[ind] = 
          kp1 * ztmin[indt] + kp * ztmax[indt] + 
          ip1* zjmin[indik] + ip * zjmin[i+j+n*k] + 
          ip * zdiag[indjk] + jp * zdiag[n1-i+n*k] + jp * zimin[i+j + n*k] +
          ip1 * zimin[indjk] 
          - ip1* ( kp1*ztmin[i] + kp * ztmax[i] ) 
          - ip * ( kp1*ztmin[i+j] + kp * ztmax[i+j]) 
          - ip * ( kp1*ztmin[inddiagj] + kp * ztmax[inddiagj])
          - jp * ( kp1*ztmin[inddiagi] + kp * ztmax[inddiagi] )
          - jp * ( kp1*ztmin[(2*n*(i+j)-(i+j)*(i+j)+(i+j))/2] + kp * ztmax[(2*n*(i+j)-(i+j)*(i+j)+(i+j))/2] ) 
          - ip1*( kp1* ztmin[(2*n*j-j*j+j)/2] + kp * ztmax[(2*n*j-j*j+j)/2] ) 
          - ip1* zjmin[k*n] - ip * zdiag[k*n] - jp * zimin[n1+k*n]
          + ip1*( kp1 * zjmin[0] + kp * zjmin[p1*n]) 
          +       ip *( kp1 * zdiag[0] + kp * zdiag[p1*n]) 
          +       jp *( kp1 * zimin[n1] + kp * zimin[n1+p1*n]); 
        ind++;
      }
  } 
  E_Int et = 0;
  E_Int off1 = 0;
  E_Int off2;
  E_Int off3 = n*(n+1)/2;
  E_Int* cnp1 = cn.begin(1);
  E_Int* cnp2 = cn.begin(2);
  E_Int* cnp3 = cn.begin(3);
  E_Int* cnp4 = cn.begin(4);
  E_Int* cnp5 = cn.begin(5);
  E_Int* cnp6 = cn.begin(6);
  for (E_Int k = 0; k < p1; k++)
  {
    for (E_Int j = n1-1; j >= 0; j--)
    {      
      off2 = off1 + j + 2; 
      for (E_Int i = 0; i < j; i++)
      {
        cnp1[et] = off1 + i + 1;
        cnp2[et] = off1 + i + 2;
        cnp3[et] = off2 + i + 1;
        cnp4[et] = off1 + i + 1 + off3;
        cnp5[et] = off1 + i + 2 + off3;
        cnp6[et] = off2 + i + 1 + off3;
        et++;
        cnp1[et] = off1 + i + 2;
        cnp2[et] = off2 + i + 1;
        cnp3[et] = off2 + i + 2;
        cnp4[et] = off1 + i + 2 + off3;
        cnp5[et] = off2 + i + 1 + off3;
        cnp6[et] = off2 + i + 2 + off3;
        et++;
      }
      cnp1[et] = off1 + j + 1;
      cnp2[et] = off1 + j + 2;
      cnp3[et] = off2 + j + 1;
      cnp4[et] = off1 + j + 1 + off3;
      cnp5[et] = off1 + j + 2 + off3;
      cnp6[et] = off2 + j + 1 + off3;
      et++;
      off1 = off2;
    }
    off1++;
  }
  
  for (E_Int nos = 0; nos < ns0; nos++)
    RELEASESHAREDS(objs[nos], structF[nos]);
  for (E_Int nos = 0; nos < nu0; nos++)
    RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

  return tpl;
}

//=============================================================================
short K_GENERATOR::checkContinuousBndsPENTA(E_Int n, E_Int p,
                                            E_Float* xtmin, E_Float* ytmin, E_Float* ztmin,
                                            E_Float* xtmax, E_Float* ytmax, E_Float* ztmax,
                                            E_Float* ximin, E_Float* yimin, E_Float* zimin,
                                            E_Float* xjmin, E_Float* yjmin, E_Float* zjmin,
                                            E_Float* xdiag, E_Float* ydiag, E_Float* zdiag)
{
  E_Float eps = 1.e-6;
  E_Int i, ind;
  E_Float dx, dy, dz;

  for (i = 0; i < n; i++)
  {
    ind = (2*n*i - i*i + i)/2;
    dx = xtmin[ind] - ximin[i];
    dy = ytmin[ind] - yimin[i];
    dz = ztmin[ind] - zimin[i];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: tmin and imin contours must be continuous.\n");
      return 0;
    }
  }

  for (i = 0; i < n; i++)
  {
    ind = (2*n*i - i*i + i)/2;
    dx = xtmax[ind] - ximin[i + (p-1)*n];
    dy = ytmax[ind] - yimin[i + (p-1)*n];
    dz = ztmax[ind] - zimin[i + (p-1)*n];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: tmax and imin contours must be continuous.\n");
      return 0;
    }
  }

  for (i = 0; i < n; i++)
  {
    dx = xtmin[i] - xjmin[i];
    dy = ytmin[i] - yjmin[i];
    dz = ztmin[i] - zjmin[i];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: tmin and jmin contours must be continuous.\n");
      return 0;
    }
  }

  for (i = 0; i < n; i++)
  {
    dx = xtmax[i] - xjmin[i + (p-1)*n];
    dy = ytmax[i] - yjmin[i + (p-1)*n];
    dz = ztmax[i] - zjmin[i + (p-1)*n];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: tmax and jmin contours must be continuous.\n");
      return 0;
    }
  }
  
  for (i = 0; i < n; i++)
  {
    ind = (2*n*i-i*i+i)/2+n-i-1;
    dx = xtmin[ind] - xdiag[i];
    dy = ytmin[ind] - ydiag[i];
    dz = ztmin[ind] - zdiag[i];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: tmin and diag contours must be continuous.\n");
      return 0;
    }
  }

  for (i = 0; i < n; i++)
  {
    ind = (2*n*i-i*i+i)/2+n-i-1;
    dx = xtmax[ind] - xdiag[i+ (p-1)*n];
    dy = ytmax[ind] - ydiag[i+ (p-1)*n];
    dz = ztmax[ind] - zdiag[i+ (p-1)*n];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: tmax and diag contours must be continuous.\n");
      return 0;
    }
  }

  for (i = 0; i < p; i++)
  {
    dx = ximin[i*n] - xjmin[i*n];
    dy = yimin[i*n] - yjmin[i*n];
    dz = zimin[i*n] - zjmin[i*n];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: imin and jmin contours must be continuous.\n");
      return 0;
    }
  }

  for (i = 0; i < p; i++)
  {
    dx = xjmin[i*n+n-1] - xdiag[i*n];
    dy = yjmin[i*n+n-1] - ydiag[i*n];
    dz = zjmin[i*n+n-1] - zdiag[i*n];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: jmin and diag contours must be continuous.\n");
      return 0;
    }
  }

  for (i = 0; i < p; i++)
  {
    dx = ximin[i*n+n-1] - xdiag[i*n+n-1];
    dy = yimin[i*n+n-1] - ydiag[i*n+n-1];
    dz = zimin[i*n+n-1] - zdiag[i*n+n-1];
    if (K_FUNC::E_abs(dx) > eps || K_FUNC::E_abs(dy) > eps || K_FUNC::E_abs(dz) > eps)
    {
      printf("Error: TFI: imin and diag contours must be continuous.\n");
      return 0;
    }
  }

  return 1;
}
