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

# include "generator.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
PyObject* K_GENERATOR::computeEta(PyObject* self, PyObject* args)
{
  PyObject *arrayc, *arrayn;
  E_Int loop, niter;
  if (!PYPARSETUPLE_(args, OO_ II_,
                    &arrayc, &arrayn, &niter, &loop)) return NULL;

  // Check array of contour
  E_Int ni, nj, nk;
  FldArrayF* coords; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(arrayc, varString, coords, ni, nj, nk, cn, eltType);
  E_Int err = 0;
  if (res != 1){ PyErr_SetString(PyExc_TypeError, "computeEta: 1st array must be structured."); err = 1;}
  if (ni == 1 || nj > 1 || nk > 1) { PyErr_SetString(PyExc_TypeError, "computeEta: 1st arg must be a 1D-array."); err = 1;}

  // check coords
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {err = 1; PyErr_SetString(PyExc_TypeError, "computeEta: 1st arg must contain coordinates.");}
  if (err == 1) {RELEASESHAREDB(res, arrayc, coords, cn); return NULL;}
  posx++; posy++; posz++;
  E_Int api = coords->getApi();
  
  // recup des normales
  E_Int nin, njn, nkn;
  FldArrayF* surf; FldArrayI* cnn;
  char* varStringn; char* eltTypen;
  E_Int resn = K_ARRAY::getFromArray3(arrayn, varStringn, surf, nin, njn, nkn, cnn, eltTypen);
  if (resn != 1){ PyErr_SetString(PyExc_TypeError, "computeEta: 2nd array must be structured."); err = 1;}
  if (nin == 1 || njn > 1 || nkn > 1) { PyErr_SetString(PyExc_TypeError, "computeEta: 2nd arg must be a 1D-array."); err = 1;}
  if (nin != ni) { PyErr_SetString(PyExc_TypeError, "computeEta: 1st and 2nd args must be of same size."); err = 1;}
  if (surf->getNfld() != 3) { PyErr_SetString(PyExc_TypeError, "computeEta: 2nd arg must be a vector."); err = 1;}
  if (err == 1) {RELEASESHAREDB(res, arrayc, coords, cn); RELEASESHAREDB(resn, arrayn, surf, cnn); return NULL;}
  
  // compute Eta
  PyObject* tpl2 = K_ARRAY::buildArray3(3, "etax,etay,etaz", ni, 1, 1, api);
  FldArrayF* eta;
  K_ARRAY::getFromArray3(tpl2, eta);
  computeEta(ni, coords->begin(posx), coords->begin(posy), coords->begin(posz),
             surf->begin(1), surf->begin(2), surf->begin(3),
             loop, niter, eta->begin(1), eta->begin(2), eta->begin(3));

  RELEASESHAREDS(tpl2, eta);
  RELEASESHAREDS(arrayc, coords);  RELEASESHAREDS(arrayn, surf); 
  return tpl2;
}
/*=============================================================================
IN: nic, xc, yc, zc : contour
IN: sxt, syt, szt : surfaces
IN: nxt, nyt, nzt : normales aux surfaces
IN: niter = 0 par defaut nb d ites de lissage
IN: loop = 1 boucle, 0 sinon
OUT: eta: doit etre alloue a nic,1,1
=============================================================================*/
void K_GENERATOR::computeEta(E_Int nic, E_Float* xc, E_Float* yc, E_Float* zc, 
                             E_Float* nxc, E_Float* nyc, E_Float* nzc, 
                             E_Int loop, E_Int niter,
                             E_Float* etax, E_Float* etay, E_Float* etaz)
{
  FldArrayF ksi(nic,3);// vecteur ksi dans la direction du contour
  E_Int nksic = nic-1;  FldArrayF ksic(nksic,3);
  E_Int niksi, njksi, nkksi;
  E_Float* ksixc = ksic.begin(1);
  E_Float* ksiyc = ksic.begin(2);
  E_Float* ksizc = ksic.begin(3);

#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < nksic-1; i++)
  {
    ksixc[i] = xc[i+1]-xc[i];
    ksiyc[i] = yc[i+1]-yc[i];
    ksizc[i] = zc[i+1]-zc[i];
  }
  ksixc[nksic-1] = xc[nksic-1]-xc[nksic-2];
  ksiyc[nksic-1] = yc[nksic-1]-yc[nksic-2];
  ksizc[nksic-1] = zc[nksic-1]-zc[nksic-2];
  if (loop == 1)
  {
    ksixc[nksic-1] = ksixc[0];
    ksiyc[nksic-1] = ksiyc[0];
    ksizc[nksic-1] = ksizc[0];
  }
  K_LOC::center2nodeStruct(ksic, nksic, 1, 1, -1, 0, -1, -1, -1, 
                           ksi, niksi, njksi, nkksi);
  normalize(niksi, ksi.begin(1), ksi.begin(2), ksi.begin(3));
  
  // calcul de eta ; produit vectoriel de n avec ksi
  E_Float* ksix = ksi.begin(1);
  E_Float* ksiy = ksi.begin(2);
  E_Float* ksiz = ksi.begin(3);

#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < nic; i++)
  {
    etax[i] = nyc[i]*ksiz[i]-nzc[i]*ksiy[i];
    etay[i] = nzc[i]*ksix[i]-nxc[i]*ksiz[i];
    etaz[i] = nxc[i]*ksiy[i]-nyc[i]*ksix[i];
  }
  normalize(nic, etax, etay, etaz);

  // lissage
  //E_Float eps = 0.1; 
  E_Int iter = 1;
  E_Float lambda0 = 0.15;
  E_Int im, ip;
  E_Float lx, ly, lz, vx0, vx1, vx2, vy0, vy1, vy2, vz0, vz1, vz2;
  vector<E_Float> etax2(nic);  vector<E_Float> etay2(nic);  vector<E_Float> etaz2(nic);  
  while (iter < niter)
  {
    iter++;
    for (E_Int i = 0; i < nic; i++)
    {etax2[i] = etax[i]; etay2[i] = etay[i]; etaz2[i] = etaz[i];}

    for (E_Int i = 0; i < nic-1; i++)
    {
      im = i-1; ip = i+1;
      if (i == 0)
      {
        if (loop == 1) im = nic-2;
        else im = i;
      }
      vx0 = etax[im]; vy0 = etay[im]; vz0 = etaz[im];
      vx1 = etax[i]; vy1 = etay[i]; vz1 = etaz[i];
      vx2 = etax[ip]; vy2 = etay[ip]; vz2 = etaz[ip];
      lx  = (vx0-2.*vx1+vx2); etax2[i] = vx1 + lambda0*lx;
      ly  = (vy0-2.*vy1+vy2); etay2[i] = vy1 + lambda0*ly;
      lz  = (vz0-2.*vz1+vz2); etaz2[i] = vz1 + lambda0*lz;
    }
    if (loop == 1) {etax2[nic-1] = etax2[0]; etay2[nic-1] = etay2[0]; etaz2[nic-1] = etaz2[0];}                                                                                     
    for (E_Int i = 0; i < nic; i++)
    {etax[i] = etax2[i]; etay[i] = etay2[i]; etaz[i] = etaz2[i];}
    normalize(nic, &etax[0], &etay[0], &etaz[0]);
  }
  return;
}
// ============================================================================
PyObject* K_GENERATOR::straightenVector(PyObject* self, PyObject* args)
{
  PyObject *arrayc, * arrayv, *constrainedPtsa, *constraintsa;
  E_Int loop, niter;
  E_Float toldist;
  if (!PYPARSETUPLE_(args, OOOO_ II_ R_,
                    &arrayc, &arrayv, &constrainedPtsa, &constraintsa, 
                    &loop, &niter, &toldist)) return NULL;
  if (PyList_Check(constrainedPtsa) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "straightenVector: 3rd argument must be a list.");
    return NULL;
  }

  if (PyList_Size(constraintsa) == 0 || PyList_Size(constrainedPtsa) == 0)
    return arrayv;

  // Check array of contour
  E_Int ni, nj, nk;
  FldArrayF* coords;
  char* varString;
  char* eltType;
  FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(arrayc, varString, coords, ni, nj, nk, cn, eltType);
  E_Int err = 0;
  if ( res != 1 ){ PyErr_SetString(PyExc_TypeError, "straightenVector: 1st array must be structured."); err = 1;}
  if ( ni == 1 || nj > 1 || nk > 1 ) { PyErr_SetString(PyExc_TypeError, "straightenVector: 1st arg must be a 1D-array."); err = 1;}
  // check coords
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {err = 1; PyErr_SetString(PyExc_TypeError, "straightenVector: 1st arg must contain coordinates.");}
  if (err == 1) {RELEASESHAREDB(res, arrayc, coords, cn); return NULL;}
  posx++; posy++; posz++;
  E_Int api = coords->getApi();

  // Check array of vector
  E_Int niv, njv, nkv;
  FldArrayF* vectp;
  FldArrayI* cnv;
  char* varStringv;
  E_Int resv = K_ARRAY::getFromArray3(arrayv, varStringv, vectp, niv, njv, nkv, cnv, eltType);
  if ( resv != 1 ) { PyErr_SetString(PyExc_TypeError, "straightenVector: 2nd array must be structured."); err = 1;}
  if ( niv == 1 || njv > 1 || nkv > 1 ) { PyErr_SetString(PyExc_TypeError, "straightenVector: 2nd arg must be a 1D-array."); err = 1;}
  if ( niv != ni ) { PyErr_SetString(PyExc_TypeError, "straightenVector: 1st and 2nd args must be of same size."); err = 1;}
  if ( vectp->getNfld() != 3 ) { PyErr_SetString(PyExc_TypeError, "straightenVector: 2nd arg must be a vector."); err = 1;}
  if ( err == 1 ) {RELEASESHAREDB(res, arrayc, coords, cn); RELEASESHAREDB(resv, arrayv, vectp, cnv); return NULL;}
  
  // Recup des indices 
  E_Int nconsi = PyList_Size(constrainedPtsa);
  vector<E_Int> constrainedPts(nconsi);
  for (int i = 0; i < nconsi; i++)
  {
    PyObject* tpl = PyList_GetItem(constrainedPtsa, i);
    if ( !PyLong_Check(tpl) && !PyInt_Check(tpl))    
    {
      RELEASESHAREDB(resv, arrayv, vectp, cnv); RELEASESHAREDB(res, arrayc, coords, cn); 
      PyErr_SetString(PyExc_TypeError, "straightenVector: 3rd arg must be a list of integers.");
      return NULL;
    }
    else constrainedPts[i] = PyLong_AsLong(tpl);
  }    
  // recuperations des contraintes
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> constraints;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = true;
  E_Bool skipDiffVars = true;
  K_ARRAY::getFromArrays(
    constraintsa, resl, structVarString, unstrVarString,
    constraints, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  
  // Redressage des normales
  PyObject* tpl2 = K_ARRAY::buildArray3(3, varStringv, ni, 1, 1, api);
  FldArrayF* vect2;
  K_ARRAY::getFromArray3(tpl2, vect2);
  straightenVector(ni, coords->begin(posx), coords->begin(posy), coords->begin(posz), 
                   loop, niter, constrainedPts, constraints, *vect2, toldist);
  RELEASESHAREDS(tpl2, vect2);
  RELEASESHAREDB(resv, arrayv, vectp, cnv); RELEASESHAREDB(res, arrayc, coords, cn); 
  for (size_t v = 0; v < objst.size(); v++) RELEASESHAREDS(objst[v], constraints[v]);

  return tpl2;
}
//=============================================================================
// A mettre en commun avec Converter.normalize dans CompGeom ensuite ?
//=============================================================================
void K_GENERATOR::normalize(E_Int npts, E_Float* vx, E_Float* vy, E_Float* vz)
{
#pragma omp parallel default(shared)
  {
    E_Float norm;
#pragma omp for
    for (E_Int i = 0; i < npts; i++)
    {
      norm = sqrt(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
      norm = 1./K_FUNC::E_max(norm, 1.e-12);
      vx[i] *= norm; vy[i] *= norm; vz[i] *= norm; 
    }
  }
}
//=============================================================================
E_Int K_GENERATOR::getClosestIndex(E_Int ni, E_Float x, E_Float y, E_Float z, 
                                   E_Float* xc, E_Float* yc, E_Float* zc)
{
  E_Int isav = -1;
  E_Float dx, dy, dz, dist;
  E_Float dmin = K_CONST::E_MAX_FLOAT;
  for (E_Int i = 0; i < ni; i++)
  {
    dx = K_FUNC::E_abs(x-xc[i]);
    dy = K_FUNC::E_abs(y-yc[i]);
    dz = K_FUNC::E_abs(z-zc[i]);
    dist = dx*dx+dy*dy+dz*dz;
    if (dist < dmin) 
    {dmin = dist; isav = i;}
  }
  return isav;
}
//=============================================================================
void K_GENERATOR::getConstrainedEta(
  E_Float etax0, E_Float etay0, E_Float etaz0, E_Int indA, 
  E_Int nic, E_Float* xc, E_Float* yc, E_Float* zc, 
  E_Float* etam, E_Float toldist)
{
  E_Int ni2 = nic-1;// nic : nb de pts dans la contrainte de coords xc, yc, zc
  E_Int indB = indA-1;
  if (indA == ni2) 
  {
    // determination si c est une boucle
    E_Float dx = xc[ni2]-xc[0];
    E_Float dy = yc[ni2]-yc[0];
    E_Float dz = zc[ni2]-zc[0];
    if (K_FUNC::E_abs(dx) < toldist && K_FUNC::E_abs(dy) < toldist && K_FUNC::E_abs(dz) < toldist) indB = 1;
    else indB = indA-1; // reverse
  }
  else indB = indA+1;
  E_Float etax1 = xc[indB]-xc[indA];
  E_Float etay1 = yc[indB]-yc[indA];
  E_Float etaz1 = zc[indB]-zc[indA];
  E_Float neta = sqrt(etax1*etax1+etay1*etay1+etaz1*etaz1);
  neta = 1./K_FUNC::E_max(neta,1e-10);
  etax1 = etax1*neta; etay1 = etay1*neta; etaz1 = etaz1*neta;
  // verif que le eta est dans le bon sens
  E_Float ps = etax0*etax1+etay1*etay1+etaz0*etaz1;
  if (ps > -toldist) {etam[0] = etax1; etam[1]=etay1;etam[2]=etaz1;}
  else {etam[0]=-etax1; etam[1]=-etay1;etam[2]=-etaz1;}
}

//=============================================================================
/* modifies the vector vx,vy,vz at points connected to constraints
   warning vx,vy,vz must be normalized */
//=============================================================================
void K_GENERATOR::applyConstraintsOnVector(
  E_Int ni, E_Float* xt, E_Float* yt, E_Float* zt,
  E_Int loop, vector<E_Int>& constrainedPts, 
  vector<FldArrayF*>& constraints,
  E_Float* vx, E_Float* vy, E_Float* vz, E_Float toldist)
{
  E_Int ncons = constraints.size();
  if (ncons == 0) return;
  E_Float etam[3];
  for (E_Int noc = 0; noc < ncons; noc++)
  {
    FldArrayF& cons = *constraints[noc]; 
    E_Int ind = constrainedPts[noc];
    if (ind != -1)
    {
      E_Int indc = getClosestIndex(ni, xt[ind], yt[ind], zt[ind], cons.begin(1), cons.begin(2), cons.begin(3));//indice du pt correspondant sur la contrainte
      getConstrainedEta(vx[ind],vy[ind],vz[ind],indc,cons.getSize(), cons.begin(1),cons.begin(2),cons.begin(3), etam, toldist);
      vx[ind] = etam[0]; vy[ind] = etam[1]; vz[ind] = etam[2];
      if (loop == 1)
      {
        if (ind==0){vx[ni-1] = vx[0]; vy[ni-1] = vy[0]; vz[ni-1] = vz[0];}
        else if (ind==ni-1){vx[0] = vx[ni-1]; vy[0] = vy[ni-1]; vz[0] = vz[ni-1];}
      }
    }  
  }
}
//=============================================================================
/*Pres de la contrainte, pour eviter les focalisations, redresse les normales
 des voisins. 
 Par contre dans le cas ou le pt contraint est un point de rebroussement
 ne pas modifier
 IN: ni: nb de pts ds vx, vy, vz */
//=============================================================================
void K_GENERATOR::relaxNearConstraints(
  E_Int ni, E_Float* xt, E_Float* yt, E_Float* zt,
  E_Int loop, vector<E_Int>& constrainedPts, 
  E_Float* alp,  E_Float* vx, E_Float* vy, E_Float* vz, E_Float dalphamax)
{
  E_Float psmax = 0.9;
  for (E_Int i = 0; i < ni-1; i++)
  {
    E_Int ip = i+1;
    E_Int found = 0;
    for (size_t noj = 0; noj < constrainedPts.size(); noj++)
    {
      E_Int j = constrainedPts[noj];
      if ( i == j ) {found = 1; break;}
    }
    if (found == 1)
    {
      E_Float alpha = alp[i];
      if (K_FUNC::E_abs(alpha-180.) < dalphamax)
      {
        E_Float xA1 = xt[i];  E_Float yA1 = yt[i]; E_Float zA1 = zt[i];
        E_Float xA2 = xt[ip]; E_Float yA2 = yt[ip]; E_Float zA2 = zt[ip];
        E_Float xB1 = xA1 + vx[i]; E_Float xB2 = xA2 + vx[ip]; 
        E_Float yB1 = yA1 + vy[i]; E_Float yB2 = yA2 + vy[ip]; 
        E_Float zB1 = zA1 + vz[i]; E_Float zB2 = zA2 + vz[ip]; 
        E_Float xA1A2 = xA2-xA1; E_Float yA1A2 = yA2-yA1; E_Float zA1A2 = zA2-zA1;
        E_Float xB1B2 = xB2-xB1; E_Float yB1B2 = yB2-yB1; E_Float zB1B2 = zB2-zB1;
        E_Float normB = sqrt(xB1B2*xB1B2+yB1B2*yB1B2+zB1B2*zB1B2);
        E_Float normA = sqrt(xA1A2*xA1A2+yA1A2*yA1A2+zA1A2*zA1A2);
        E_Float ps = xA1A2*xB1B2+yA1A2*yB1B2+zA1A2*zB1B2;
        ps = ps/(K_FUNC::E_max(normA*normB,1.e-12));
        if (K_FUNC::E_abs(ps)<psmax)
        {
          vx[ip] = xB1+xA1A2-xA2;
          vy[ip] = yB1+yA1A2-yA2;
          vz[ip] = zB1+zA1A2-zA2;
        }
      }
    }
  }
}
//=============================================================================
/* ni: taille de vx, vy, vz */
//=============================================================================
void K_GENERATOR::straightenVector(E_Int ni, E_Float* xt, E_Float* yt, E_Float* zt,
                                   E_Int loop, E_Int niter, vector<E_Int>& constrainedPts, 
                                   vector<FldArrayF*>& constraints,
                                   FldArrayF& vect, E_Float toldist)
{
  E_Float* vx = vect.begin(1); E_Float* vy = vect.begin(2); E_Float* vz = vect.begin(3);
  if (constraints.size() == 0 || constrainedPts.size() == 0) return;

  E_Float lambda0 = 0.15;
  E_Int iter = 0;
  //1. normalize vect
  normalize(ni, vx, vy, vz);
  
  //2. calcul des contraintes si existence
  applyConstraintsOnVector(ni, xt, yt, zt, loop, constrainedPts, constraints, vx, vy, vz, toldist);

  //3. calcul de l angle max
  E_Float dirVect[3];
  vector<E_Float> alpt(ni);
  E_Float ptA1[3]; E_Float ptB1[3]; E_Float ptA2[3]; E_Float ptB2[3]; 
  dirVect[0] = 0;  dirVect[1] = 0;  dirVect[2] = 1; 
  E_Int i, ip, im;
  E_Float alpha0;
  for (i = 1; i < ni-1; i++)
  {
    im = i-1; ip = i+1;
    ptA1[0] = xt[i]; ptA1[1] = yt[i]; ptA1[2] = zt[i];
    ptB1[0] = xt[im]; ptB1[1] = yt[im]; ptB1[2] = zt[im];
    ptA2[0] = xt[i]; ptA2[1] = yt[i]; ptA2[2] = zt[i];
    ptB2[0] = xt[ip]; ptB2[1] = yt[ip]; ptB2[2] = zt[ip];
    alpha0 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA1, ptB1, ptA2, ptB2, dirVect);
    if (alpha0 != -1000.) alpt[i] = alpha0;
    else alpt[i] = 0.;
  }
  if (loop == 1) // loop
  {
    i = 0; ip = 1; im = ni-2;
    ptA1[0] = xt[i]; ptA1[1] = yt[i]; ptA1[2] = zt[i];
    ptB1[0] = xt[im]; ptB1[1] = yt[im]; ptB1[2] = zt[im];
    ptA2[0] = xt[i]; ptA2[1] = yt[i]; ptA2[2] = zt[i];
    ptB2[0] = xt[ip]; ptB2[1] = yt[ip]; ptB2[2] = zt[ip];
    alpha0 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA1, ptB1, ptA2, ptB2, dirVect);
    if (alpha0 != -1000.) {alpt[0] = alpha0; alpt[ni-1] = alpha0;}
  }

  E_Float dalphamax = 170.;
  vector<E_Float> vxt(ni);  vector<E_Float> vyt(ni); vector<E_Float> vzt(ni);
  FldArrayF vectc(ni-1,3);
  E_Int nin, njn, nkn;
  while (iter < niter)
  {
    iter++;
    //1. lissage de vect
    K_LOC::node2centerStruct(vect, ni, 1, 1, -1, 0, vectc);
    K_LOC::center2nodeStruct(vectc, ni-1, 1, 1, -1, 0, -1, -1, -1, vect, nin, njn, nkn);
    //2.normalize vect
    normalize(ni, vx, vy, vz);

    for (i = 0; i < ni-1; i++)
    {
      im = i-1; ip = i+1;
      if (i==0) 
      {
        if (loop == 1) im = ni-1;
        else im = i;
      }

      E_Float vx0 = vx[im]; E_Float vy0 = vy[im]; E_Float vz0 = vz[im];
      E_Float vx1 = vx[i];  E_Float vy1 = vy[i];  E_Float vz1 = vz[i];
      E_Float vx2 = vx[ip]; E_Float vy2 = vy[ip]; E_Float vz2 = vz[ip];
      vxt[i] = vx1 + lambda0*(vx0-2.*vx1+vx2);
      vyt[i] = vy1 + lambda0*(vy0-2.*vy1+vy2);
      vzt[i] = vz1 + lambda0*(vz0-2.*vz1+vz2);
    }// fin boucle i
  
    if (loop == 1)
    {
      im = ni-2; i = ni-1; ip = 0;
      E_Float vx0 = vx[im]; E_Float vy0 = vy[im]; E_Float vz0 = vz[im];
      E_Float vx1 = vx[i];  E_Float vy1 = vy[i];  E_Float vz1 = vz[i];
      E_Float vx2 = vx[ip]; E_Float vy2 = vy[ip]; E_Float vz2 = vz[ip];
      vxt[i] = vx1 + lambda0*(vx0-2.*vx1+vx2);
      vyt[i] = vy1 + lambda0*(vy0-2.*vy1+vy2);
      vzt[i] = vz1 + lambda0*(vz0-2.*vz1+vz2);
    }
    else{vxt[ni-1] = vx[ni-1]; vyt[ni-1] = vy[ni-1]; vzt[ni-1] = vz[ni-1];}    
    vx = &vxt[0]; vy = &vyt[0]; vz = &vzt[0];
    // normalize
    normalize(ni, vx, vy, vz);
    applyConstraintsOnVector(ni, xt, yt, zt, loop, constrainedPts, constraints, vx, vy, vz, toldist);
    relaxNearConstraints(ni, xt, yt, zt, loop, constrainedPts, &alpt[0], vx, vy, vz, dalphamax);
  }//fin boucle while
                
  //5. normalize
  normalize(ni, vx, vy, vz);  
  for (E_Int i = 0; i < ni; i++)
  {
    vect(i,1) = vx[i];
    vect(i,2) = vy[i];
    vect(i,3) = vz[i];
  }
}
