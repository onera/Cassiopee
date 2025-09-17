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
# include "CompGeom/compGeom.h"
using namespace K_FLD;


//=============================================================================
/* Generates an hexa mesh starting from a1 the surface quad, a2 the external 
   quad, h the height of the mesh, hext the height of the mesh at external
   borders, hf the height near the wall, density the number of points in 
   the normal direction to the wall */
//=============================================================================
PyObject* K_GENERATOR::front2Hexa(PyObject* self, PyObject* args)
{
  PyObject *a1, *a2, *distrib;
  E_Float h0;
  if (!PYPARSETUPLE_(args, OOO_ R_, &a1, &a2, &distrib, &h0))
  {
    return NULL;
  }
  // Check quad1: front projete sur le corps
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray3(a1, varString1, f1, ni1, nj1, nk1, 
                                      cn1, eltType1);
  if (res1 != 2 || strcmp(eltType1, "QUAD") != 0) 
  {
    RELEASESHAREDB(res1, a1, f1, cn1);
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: a1 must be a QUAD.");
    return NULL;
  }
  // Check quad2: exterieur du front
  E_Int ni2, nj2, nk2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray3(a2, varString2, f2, ni2, nj2, nk2, cn2, eltType2);
  if (res2 != 2 || strcmp(eltType2, "QUAD") != 0) 
  {
    RELEASESHAREDB(res1, a1, f1, cn1);
    RELEASESHAREDB(res2, a2, f2, cn2);    
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: a2 must be a QUAD.");
    return NULL;
  }
  // Check distrib
  E_Int nid, njd, nkd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringd; char* eltTyped;
  E_Int resd = K_ARRAY::getFromArray3(distrib, varStringd, fd, nid, njd, nkd, 
                                      cnd, eltTyped);
  if (resd != 1) 
  {
    delete f1; delete cn1; delete f2; delete cn2;
    RELEASESHAREDU(a1, f1, cn1);
    RELEASESHAREDU(a2, f2, cn2);
    RELEASESHAREDB(res2, distrib, fd, cnd);
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: distrib must be structured.");
    return NULL;
  }
  if (nid < 2 || njd != 1 || nkd != 1) 
  {
    RELEASESHAREDU(a1, f1, cn1);
    RELEASESHAREDU(a2, f2, cn2);
    RELEASESHAREDS(distrib, fd);
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: distrib must be an i-array.");
    return NULL;
  }
  
  E_Int npts1 = f1->getSize(); E_Int nelts1 = cn1->getSize();
  E_Int api = f1->getApi();
  if (f2->getSize() != npts1 || cn2->getSize() != nelts1) 
  {
    RELEASESHAREDU(a1, f1, cn1);
    RELEASESHAREDU(a2, f2, cn2);
    RELEASESHAREDS(distrib, fd);
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: a1 and a2 must be of same size.");
    return NULL;    
  }
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    RELEASESHAREDU(a1, f1, cn1);
    RELEASESHAREDU(a2, f2, cn2);
    RELEASESHAREDS(distrib, fd);
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: coords must be present in a1.");
    return NULL;       
  }
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    RELEASESHAREDU(a1, f1, cn1);
    RELEASESHAREDU(a2, f2, cn2);
    RELEASESHAREDS(distrib, fd);
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: coords must be present in a2.");
    return NULL;
  }
  E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
  E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
  E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  if (posxd == -1 || posyd == -1 || poszd == -1)
  {
    RELEASESHAREDU(a1, f1, cn1);
    RELEASESHAREDU(a2, f2, cn2);
    RELEASESHAREDS(distrib, fd);
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: coords must be present in distrib.");
    return NULL;       
  }
  posx1++; posy1++; posz1++; posx2++; posy2++; posz2++; 
  posxd++; posyd++; poszd++;

  E_Float* xt1 = f1->begin(posx1);
  E_Float* yt1 = f1->begin(posy1);
  E_Float* zt1 = f1->begin(posz1);
  E_Float* xt2 = f2->begin(posx2);
  E_Float* yt2 = f2->begin(posy2);
  E_Float* zt2 = f2->begin(posz2);
  E_Int nk = nid;
  FldArrayF dtb(nk, 3);
  dtb.setOneField(*fd, posxd, 1);
  dtb.setOneField(*fd, posyd, 2);
  dtb.setOneField(*fd, poszd, 3);

  E_Float dx, dy, dz;
  E_Float xs, ys, zs, xe, ye, ze;
  E_Int nelts = nelts1*(nk-1); 
  FldArrayF* f = new FldArrayF(npts1*nk, 3);// coords hexa mesh
  FldArrayI* cn = new FldArrayI(nelts, 8);// connectivite hexa mesh
  FldArrayF line(2,3);  
  FldArrayF lineMap(nk,3);
  FldArrayF s(nk); FldArrayF dtbx(nk); FldArrayF dtby(nk); FldArrayF dtbz(nk);

  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);
  E_Float xe2, ye2, ze2;
  E_Float norm;

  for (E_Int vert = 0; vert < npts1; vert++)
  {
    xs = xt1[vert]; ys = yt1[vert]; zs = zt1[vert];
    xe = xt2[vert]; ye = yt2[vert]; ze = zt2[vert];    
    dx = xs-xe; dy = ys-ye; dz = zs-ze;
    norm = 1./sqrt(dx*dx + dy*dy + dz*dz);
    xe2 = h0 * dx * norm + xs;
    ye2 = h0 * dy * norm + ys;
    ze2 = h0 * dz * norm + zs;
    // dbg
    xe2 = xe; ye2 = ye; ze2 = ze;

    // ligne P1P'1
    line(0,1) = xs; line(0,2) = ys; line(0,3) = zs;
    line(1,1) = xe2; line(1,2) = ye2; line(1,3) = ze2;          
    // map dtb on line
    K_COMPGEOM::onedmap(line.getSize(),
			line.begin(1), line.begin(2), line.begin(3),
			nk, dtb.begin(),
			lineMap.begin(1), lineMap.begin(2), lineMap.begin(3),
			s.begin(), dtbx.begin(), dtby.begin(), dtbz.begin());

    for (E_Int k = 0; k < nk; k++)
    {
      xt[vert + k*npts1] = lineMap(k,1);
      yt[vert + k*npts1] = lineMap(k,2);
      zt[vert + k*npts1] = lineMap(k,3);
    }
  }
  //construction de la connectivite HEXA
  E_Int noet = 0;
  E_Int* c1 = cn->begin(1); E_Int* c2 = cn->begin(2);
  E_Int* c3 = cn->begin(3); E_Int* c4 = cn->begin(4);
  E_Int* c5 = cn->begin(5); E_Int* c6 = cn->begin(6);
  E_Int* c7 = cn->begin(7); E_Int* c8 = cn->begin(8);
  E_Int* cp1 = cn1->begin(1); E_Int* cp2 = cn1->begin(2);
  E_Int* cp3 = cn1->begin(3); E_Int* cp4 = cn1->begin(4);
  for (E_Int et = 0; et < nelts1; et++)
  {
    c1[et] = cp1[et];
    c2[et] = cp2[et];
    c3[et] = cp3[et];
    c4[et] = cp4[et];
    c5[et] = cp1[et] + npts1;
    c6[et] = cp2[et] + npts1;
    c7[et] = cp3[et] + npts1;
    c8[et] = cp4[et] + npts1;
    noet++;
  }
  for (E_Int et = noet; et < nelts; et++)
  {
    E_Int etp = et-nelts1;
    c1[et] = c5[etp];
    c2[et] = c6[etp];
    c3[et] = c7[etp];
    c4[et] = c8[etp];
    c5[et] = c5[etp] + npts1;
    c6[et] = c6[etp] + npts1;
    c7[et] = c7[etp] + npts1;
    c8[et] = c8[etp] + npts1;
  }

  RELEASESHAREDU(a1, f1, cn1);
  RELEASESHAREDU(a2, f2, cn2);
  RELEASESHAREDS(distrib, fd);
  PyObject* tpl = K_ARRAY::buildArray3(*f, "x,y,z", *cn, "HEXA", api);
  delete f; delete cn;
  return tpl;
}
