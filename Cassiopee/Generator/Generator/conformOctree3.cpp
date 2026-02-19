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

// Conformize a 27-tree

#include "generator.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Conformize the octree3
   Prend un octree de rapport 3.
   Coupe certaines cellules pour obtenir un maillage conforme */
//=============================================================================
PyObject* K_GENERATOR::conformOctree3(PyObject* self, PyObject* args)
{
  PyObject *octree;
  if (!PYPARSETUPLE_(args, O_, &octree)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(octree, varString, f, ni, nj, nk, cn, 
                                     eltType);
  if (res != 2) 
  {
    if (res == 1) RELEASESHAREDS(octree, f); 
    PyErr_SetString(PyExc_TypeError,
                    "conformOctree3: array must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType, "HEXA") != 0 && strcmp(eltType, "QUAD") != 0 )
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "conformOctree3: the octree must be HEXA or QUAD.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "conformOctree3: the octree must contain coordinates.");
    return NULL;  
  }
  posx++; posy++; posz++;
 
  /*-----------------------------------*/
  // Calcul du niveau de chaque cellule
  /*-----------------------------------*/
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  E_Int nelts = cn->getSize(); E_Int nvert = cn->getNfld();
  E_Int npts = f->getSize();
  E_Int api = f->getApi();

  FldArrayF dht(nelts);
  FldArrayI levels(nelts);
  // calcul de la bbox de l'octree
  E_Float xmin = K_CONST::E_MAX_FLOAT; E_Float xmax =-K_CONST::E_MAX_FLOAT;
  E_Float ymin = K_CONST::E_MAX_FLOAT; E_Float ymax =-K_CONST::E_MAX_FLOAT;
  E_Float zmin = K_CONST::E_MAX_FLOAT; E_Float zmax =-K_CONST::E_MAX_FLOAT;
  FldArrayF bbox(nelts, 6);
  K_COMPGEOM::boundingBoxOfUnstrCells(*cn, xt, yt, zt, bbox);
  E_Float* xmint = bbox.begin(1); E_Float* xmaxt = bbox.begin(4);
  E_Float* ymint = bbox.begin(2); E_Float* ymaxt = bbox.begin(5);
  E_Float* zmint = bbox.begin(3); E_Float* zmaxt = bbox.begin(6);
  for (E_Int et = 0; et < nelts; et++)
  {
    xmin = K_FUNC::E_min(xmin,xmint[et]); xmax = K_FUNC::E_max(xmax,xmaxt[et]);
    ymin = K_FUNC::E_min(ymin,ymint[et]); ymax = K_FUNC::E_max(ymax,ymaxt[et]);
    zmin = K_FUNC::E_min(zmin,zmint[et]); zmax = K_FUNC::E_max(zmax,zmaxt[et]);
  }
  bbox.malloc(0,1);

  // calcul du pas sur chaque grille 
  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
  E_Int ind1, ind2;
  E_Float maxdh = 0.;
  E_Float mindh = K_CONST::E_MAX_FLOAT;
  for (E_Int et = 0; et < nelts; et++)
  {
    ind1 = cn1[et]-1; ind2 = cn2[et]-1;
    dht[et] = xt[ind2]-xt[ind1];
    mindh = K_FUNC::E_min(mindh, dht[et]);
    maxdh = K_FUNC::E_max(maxdh, dht[et]);
  }

  // Initialisations 
  // Recopie des vertex
  E_Int nfld = f->getNfld();
  FldArrayF* fo = new FldArrayF(70*npts, nfld); 
  FldArrayI* cno = new FldArrayI(30*nelts,nvert); 
  FldArrayI levelso(30*nelts);
  for (E_Int v = 1; v <= nfld; v++)
  {
    E_Float* fop = fo->begin(v);
    E_Float* fp = f->begin(v);
    for (E_Int i = 0; i < npts; i++) fop[i] = fp[i]; 
  }
  E_Int l0 = 1;//niveau courant
  E_Float tol = mindh*0.1;
  while (maxdh > mindh-tol )
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      if (K_FUNC::E_abs(dht[et]-maxdh) < tol) levels[et] = l0;
    }
    maxdh = maxdh/3.; l0++;
  }
  FldArrayF floc = *f; FldArrayI cnloc = *cn;

  // Fin initialisations 
  E_Int lmax = l0;
  E_Int no = npts; 
  E_Int eto = 0;
  E_Int count = 0;
  restart:;
  E_Int restart = conformizeElements(mindh, maxdh, lmax, npts, posx, posy, posz, floc, cnloc, levels, *fo, *cno, levelso, no, eto);
  if ( restart == -1 ) 
  {
    printf("Warning: conformOctree3: cannot compute connectivity. Check connectivity."); 
    RELEASESHAREDU(octree, f, cn); return octree;
  }
  count++;
  // Realloc
  fo->reAllocMat(no,fo->getNfld()); cno->reAllocMat(eto,nvert); levelso.resize(eto);
  if (count == 100) restart = 0;//securite
//   //dbx
//   if ( count == 1) restart = 0;
//   // fin debug
  if (restart == 1) 
  {
    floc = *fo; cnloc = *cno; levels = levelso;
    npts = no; nelts = eto;
    fo->reAllocMat(70*npts, nfld);  cno->reAllocMat(30*nelts,nvert); levelso.resize(30*nelts);
    eto = 0;
    goto restart;
  }
  // Sortie
  RELEASESHAREDU(octree, f, cn);
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, *fo, *cno);

  PyObject* tpl = K_ARRAY::buildArray3(*fo, varString, *cno, eltType, api);
  delete fo; delete cno;
  return tpl;
}
//=============================================================================
E_Int K_GENERATOR::conformizeElements(
  E_Float mindh, E_Float maxdh, E_Int lmax, E_Int npts,
  E_Int posx, E_Int posy, E_Int posz,
  FldArrayF& f, FldArrayI& cn, FldArrayI& levels,
  FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, E_Int& no, E_Int& eto)
{
  E_Int restart = 0;
  E_Int nelts = cn.getSize();  E_Int nvert = cn.getNfld();
  vector<E_Int> tag(nvert);
  E_Int ntag, nvoisins, indloc, indv, etv;
  // Connectivite elts/voisins
  vector< vector<E_Int> > cEEN(nelts);
  E_Int ok = getNeighbourElts(npts, f.begin(posx), f.begin(posy), f.begin(posz), cn, cEEN, 1, mindh); 
  if (ok == 0) return -1;

  // Conformisation de chaque element
  eto = 0;

  E_Int l0 = 1;
  while (l0 < lmax)
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      ntag = 0; 
      if (K_FUNC::E_abs(levels[et]) == l0 ) 
      {
        vector<E_Int>& voisins = cEEN[et];
        nvoisins = voisins.size();
        for (E_Int novert = 1; novert <= nvert; novert++)
        {
          indloc = cn(et,novert); tag[novert-1] = 0;
          for (E_Int nv = 0; nv < nvoisins; nv++)
          {
            etv = voisins[nv];
            if (levels[etv] > levels[et] && levels[et] > 0) 
            {
              for (E_Int novertv = 1; novertv <= nvert; novertv++)
              {
                indv = cn(etv,novertv);
                if ( indv == indloc ) 
                {
                  tag[novert-1] = 1; ntag++; 
                  goto nextvert;
                }
              }
            }
          }
          nextvert:;
        }
        if (nvert == 4) // cas 2D
        {
          switch (ntag)
          {
            case 0:
              for (E_Int novert = 1; novert <= nvert; novert++)
                cno(eto,novert) = cn(et,novert);
              levelso[eto] = levels[et];
              eto++;
              break;
              
            case 1:
              insert2DPattern1(et, levels[et], &tag[0], f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;
              
            case 2:
              insert2DPattern2(et, levels[et], &tag[0], f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;
              
            case 3:
              insert2DPattern3(et, levels[et], &tag[0], f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;
              
            case 4:
              insert2DPattern4(et, levels[et], &tag[0], f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;
              
            default:
              for (E_Int novert = 1; novert <= nvert; novert++)
                cno(eto,novert) = cn(et,novert);
              levelso[eto] = levels[et];
              eto++;
              break;
          }
        } 
        else // 3D
        {
          switch (ntag)
          {
            case 0:
              for (E_Int novert = 1; novert <= nvert; novert++)
                cno(eto,novert) = cn(et,novert);
              levelso[eto] = levels[et];
              eto++;
              break;
            
            case 1:
              insert3DPattern1(et, levels[et], &tag[0], f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;

              /*
            case 2:
              ret = find2Template(tag);
              if (ret == 0)
                insert3DPattern2(et, levels[et], tag, f, cn, fo, cno, levelso, eto, no);
              else if (ret == 1)
                insert3DPattern3(et, levels[et], tag, f, cn, fo, cno, levelso, eto, no);
              else insert3DPattern4(et, levels[et], tag, f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;
              
            case 3:
              ret = find3Template(tag);
              if (ret == 0)
                insert3DPattern3(et, levels[et], tag, f, cn, fo, cno, levelso, eto, no);
              else insert3DPattern4(et, levels[et], tag, f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;

            case 4:
              ret = find4Template(tag);
              if (ret == 0) insert3DPattern3(et, levels[et], tag, f, cn, fo, cno, levelso, eto, no);
              else insert3DPattern4(et, levels[et], tag, f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;
            */
              
            case 5:
            case 6:
            case 7:
            case 8:
              insert3DPattern4(et, levels[et], &tag[0], f, cn, fo, cno, levelso, eto, no);
              restart = 1;
              break;
           
            default:
              for (E_Int novert = 1; novert <= nvert; novert++)
                cno(eto,novert) = cn(et,novert);
              levelso[eto] = levels[et];
              eto++;
              break;
          }
        }
      }
    }
    l0++;
  }
  // Nettoyage
  for (E_Int v = 0; v < nelts; v++) cEEN[v].clear();
  cEEN.clear();
  return restart;
}
//=============================================================================
void K_GENERATOR::insert2DPattern1(E_Int et, E_Int let, E_Int* tag,
                                   FldArrayF& f, FldArrayI& cn,
                                   FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, 
                                   E_Int& eto, E_Int& no)
{
  E_Int c1, c2, c3, c4;
  c1 = cn(et,1)-1; c2 = cn(et,2)-1; c3 = cn(et,3)-1; c4 = cn(et,4)-1;
  if (tag[0] == 1)
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
    no++;
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
    no++;
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;
    cno(eto,1) = c1+1; cno(eto,2) = no-2; cno(eto,3) = no; cno(eto,4) = no-1; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-2; cno(eto,2) = c2+1; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] = -let; eto++;
    cno(eto,1) = no-1; cno(eto,2) = no; cno(eto,3) = c3+1; cno(eto,4) = c4+1; levelso[eto] = -let; eto++;
  }
  else if (tag[1] == 1)
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c1,nf);
    no++;
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c3,nf);
    no++;
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;
    cno(eto,1) = c1+1; cno(eto,2) = no-2; cno(eto,3) = no; cno(eto,4) = c4+1; levelso[eto] = -let; eto++;
    cno(eto,1) = no-2; cno(eto,2) = c2+1; cno(eto,3) = no-1; cno(eto,4) = no; levelso[eto] = let+1; eto++;
    cno(eto,1) = no; cno(eto,2) = no-1; cno(eto,3) = c3+1; cno(eto,4) = c4+1; levelso[eto] = -let; eto++;
  }
  else if (tag[2] == 1)
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c2,nf);
    no++;
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c4,nf);
    no++;
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c4,nf);
    no++;
    cno(eto,1) = c1+1; cno(eto,2) = c2+1; cno(eto,3) = no-2; cno(eto,4) = no; levelso[eto] = -let; eto++;
    cno(eto,1) = no; cno(eto,2) = no-2; cno(eto,3) = c3+1; cno(eto,4) = no-1; levelso[eto] = let+1; eto++;
    cno(eto,1) = c1+1; cno(eto,2) = no; cno(eto,3) = no-1; cno(eto,4) = c4+1; levelso[eto] = -let; eto++;
  }
  else if (tag[3] == 1)
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c1,nf);
    no++;
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c3,nf);
    no++;
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf);
    no++;
    cno(eto,1) = c1+1; cno(eto,2) = c2+1; cno(eto,3) = no; cno(eto,4) = no-2; levelso[eto] = -let; eto++;
    cno(eto,1) = no; cno(eto,2) = c2+1; cno(eto,3) = c3+1; cno(eto,4) = no-1; levelso[eto] = -let; eto++;
    cno(eto,1) = no-2; cno(eto,2) = no; cno(eto,3) = no-1; cno(eto,4) = c4+1; levelso[eto] = let+1; eto++;
  }
}
//=============================================================================
void K_GENERATOR::insert2DPattern2(E_Int et, E_Int let, E_Int* tag,
                                   FldArrayF& f, FldArrayI& cn,
                                   FldArrayF& fo, FldArrayI& cno,  FldArrayI& levelso,
                                   E_Int& eto, E_Int& no)
{
  E_Int c1, c2, c3, c4;
  c1 = cn(et,1)-1; c2 = cn(et,2)-1; c3 = cn(et,3)-1; c4 = cn(et,4)-1;
  if (tag[0] == 1 && tag[1] == 1) // cas 2a
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
    no++;// no-7
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c2,nf);
    no++;// no-6
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
    no++;// no-5
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c3,nf);
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c4,nf);
    no++;// no
    cno(eto,1) = c1+1; cno(eto,2) = no-7; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-7; cno(eto,2) = no-6; cno(eto,3) = no-3; cno(eto,4) = no-4; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-6; cno(eto,2) = c2+1; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = no-1; cno(eto,4) = c4+1; levelso[eto] =-let; eto++;
    cno(eto,1) = no-4; cno(eto,2) = no-3; cno(eto,3) = no; cno(eto,4) = no-1; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] =-let; eto++;
    cno(eto,1) = no-1; cno(eto,2) = no; cno(eto,3) = c3+1; cno(eto,4) = c4+1; levelso[eto] =-let; eto++;
  }
  else if (tag[1] == 1 && tag[2] == 1)//cas 2a
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c2,nf);
    no++;// no-7
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-6
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-5
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c3,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c4,nf);
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c2,nf)+2./3.*f(c3,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c4,nf)+2./3.*f(c3,nf);
    no++;// no
    cno(eto,1) = c1+1; cno(eto,2) = no-7; cno(eto,3) = no-5; cno(eto,4) = no-6; levelso[eto] =-let; eto++;
    cno(eto,1) = no-7; cno(eto,2) = c2+1; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] = let+1; eto++;
    cno(eto,1) = c1+1; cno(eto,2) = no-6; cno(eto,3) = no-3; cno(eto,4) = c4+1; levelso[eto] =-let; eto++;
    cno(eto,1) = no-6; cno(eto,2) = no-5; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = no-1; cno(eto,4) = no-2; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = no; cno(eto,4) = c4+1; levelso[eto] =-let; eto++;
    cno(eto,1) = no-2; cno(eto,2) = no-1; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] = let+1; eto++;
  }
  else if (tag[2] == 1 && tag[3] == 1)//cas 2a
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-7
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-6
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c4,nf);
    no++;// no-5
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c4,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c2,nf)+2./3.*f(c3,nf);
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c3,nf)+2./3.*f(c4,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c4,nf)+2./3.*f(c3,nf);
    no++;// no
    cno(eto,1) = c1+1; cno(eto,2) = c2+1; cno(eto,3) = no-6; cno(eto,4) = no-7; levelso[eto] =-let; eto++;
    cno(eto,1) = c1+1; cno(eto,2) = no-7; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] =-let; eto++;
    cno(eto,1) = no-7; cno(eto,2) = no-6; cno(eto,3) = no-3; cno(eto,4) = no-4; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-6; cno(eto,2) = c2+1; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] =-let; eto++;
    cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = no-1; cno(eto,4) = c4+1; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-4; cno(eto,2) = no-3; cno(eto,3) = no; cno(eto,4) = no-1; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] = let+1; eto++;
  }
  else if (tag[3] == 1 && tag[0] == 1)//cas 2a
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
    no++;// no-7
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
    no++;// no-6
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-5
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c4,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf);
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c4,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c3,nf)+2./3.*f(c4,nf);
    no++;// no
    cno(eto,1) = c1+1; cno(eto,2) = no-7; cno(eto,3) = no-5; cno(eto,4) = no-6; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-7; cno(eto,2) = c2+1; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] =-let; eto++;
    cno(eto,1) = no-6; cno(eto,2) = no-5; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = no-1; cno(eto,4) = no-2; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-4; cno(eto,2) = c2+1; cno(eto,3) = c3+1; cno(eto,4) = no-1; levelso[eto] =-let; eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = no; cno(eto,4) = c4+1; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-2; cno(eto,2) = no-1; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] =-let; eto++;
  }
  else if ( tag[0] == 1 && tag[2] == 1 ) // cas 2b
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./4.*(f(c1,nf)+f(c2,nf)+f(c3,nf)+f(c4,nf));
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c2,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c4,nf);
    no++;// no
    cno(eto,1) = c1+1; cno(eto,2) = no-4; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] =-let; eto++;
    cno(eto,1) = no-4; cno(eto,2) = c2+1; cno(eto,3) = no-1; cno(eto,4) = no-2; levelso[eto] =-let; eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = no; cno(eto,4) = c4+1; levelso[eto] =-let; eto++;
    cno(eto,1) = no-2; cno(eto,2) = no-1; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] =-let; eto++;
  }
  else if ( tag[1] == 1 && tag[3] == 1 ) // cas 2b
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c2,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c3,nf)+2./3.*f(c2,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./4.*(f(c1,nf)+f(c2,nf)+f(c3,nf)+f(c4,nf));
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c4,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c3,nf);
    no++;// no
    cno(eto,1) = c1+1; cno(eto,2) = no-4; cno(eto,3) = no-2; cno(eto,4) = no-1; levelso[eto] =-let;eto++;
    cno(eto,1) = no-4; cno(eto,2) = c2+1; cno(eto,3) = no-3; cno(eto,4) = no-2; levelso[eto] =-let;eto++;
    cno(eto,1) = no-1; cno(eto,2) = no-2; cno(eto,3) = no; cno(eto,4) = c4+1; levelso[eto] =-let;eto++;
    cno(eto,1) = no-2; cno(eto,2) = no-3; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] =-let;eto++;
  }
}
//=============================================================================
void K_GENERATOR::insert2DPattern3(E_Int et, E_Int let, E_Int* tag,
                                   FldArrayF& f, FldArrayI& cn,
                                   FldArrayF& fo, FldArrayI& cno,  FldArrayI& levelso,
                                   E_Int& eto, E_Int& no)
{
  E_Int c1, c2, c3, c4;
  c1 = cn(et,1)-1; c2 = cn(et,2)-1; c3 = cn(et,3)-1; c4 = cn(et,4)-1;
  if (tag[3] == 0)
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
    no++;// no-9
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c2,nf);
    no++;// no-8
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
    no++;// no-7
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-6
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-5
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c3,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c1,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c4,nf);
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c2,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c4,nf);
    no++;// no
    cno(eto,1) = c1+1; cno(eto,2) = no-9; cno(eto,3) = no-6; cno(eto,4) = no-7; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-9; cno(eto,2) = no-8; cno(eto,3) = no-5; cno(eto,4) = no-6; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-8; cno(eto,2) = c2+1; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-7; cno(eto,2) = no-6; cno(eto,3) = no-3; cno(eto,4) = c4+1; levelso[eto] =-let; eto++;
    cno(eto,1) = no-6; cno(eto,2) = no-5; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = no-1; cno(eto,4) = no-2; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = no; cno(eto,4) = c4+1; levelso[eto] =-let; eto++;
    cno(eto,1) = no-2; cno(eto,2) = no-1; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] = let+1; eto++;
  }
  else if (tag[2] == 0)
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
    no++;// no-9
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c2,nf);
    no++;// no-8
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
    no++;// no-7
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-6
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-5
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c3,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c1,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c1,nf);
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c4,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c3,nf);
    no++;// no    
    cno(eto,1) = c1+1; cno(eto,2) = no-9; cno(eto,3) = no-6; cno(eto,4) = no-7; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-9; cno(eto,2) = no-8; cno(eto,3) = no-5; cno(eto,4) = no-6; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-8; cno(eto,2) = c2+1; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-7; cno(eto,2) = no-6; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-6; cno(eto,2) = no-5; cno(eto,3) = no-1; cno(eto,4) = no-2; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = c3+1; cno(eto,4) = no-1; levelso[eto] =-let; eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = no; cno(eto,4) = c4+1; levelso[eto] = let+1; eto++;
    cno(eto,1) = no-2; cno(eto,2) = no-1; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] =-let; eto++;
  }
  else if (tag[1] == 0)
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
    no++;// no-9
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
    no++;// no-8
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-7
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-6
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c1,nf);
    no++;// no-5
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c1,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c4,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c2,nf);
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c3,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c4,nf);
    no++;// no    

    cno(eto,1) = c1+1; cno(eto,2) = no-9; cno(eto,3) = no-7; cno(eto,4) = no-8; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-9; cno(eto,2) = c2+1; cno(eto,3) = no-6; cno(eto,4) = no-7; levelso[eto] =-let;eto++;
    cno(eto,1) = no-8; cno(eto,2) = no-7; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-7; cno(eto,2) = no-6; cno(eto,3) = no-3; cno(eto,4) = no-4; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-6; cno(eto,2) = c2+1; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] =-let;eto++;
    cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = no-1; cno(eto,4) = c4+1; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-4; cno(eto,2) = no-3; cno(eto,3) = no; cno(eto,4) = no-1; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] = let+1;eto++;
  }
  else if (tag[0] == 0)
  {
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c1,nf);
    no++;// no-9
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-8
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
    no++;// no-7
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c3,nf);
    no++;// no-6
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c1,nf);
    no++;// no-5
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c1,nf);
    no++;// no-4
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c4,nf);
    no++;// no-3
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c2,nf);
    no++;// no-2
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c3,nf);
    no++;// no-1
    for (E_Int nf = 1; nf <= f.getNfld(); nf++)
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c4,nf);
    no++;// no  
    cno(eto,1) = c1+1; cno(eto,2) = no-9; cno(eto,3) = no-7; cno(eto,4) = no-8; levelso[eto] =-let;eto++;
    cno(eto,1) = no-9; cno(eto,2) = c2+1; cno(eto,3) = no-6; cno(eto,4) = no-7; levelso[eto] = let+1;eto++;
    cno(eto,1) = c1+1; cno(eto,2) = no-8; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] =-let;eto++;
    cno(eto,1) = no-8; cno(eto,2) = no-7; cno(eto,3) = no-3; cno(eto,4) = no-4; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-7; cno(eto,2) = no-6; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = no-1; cno(eto,4) = c4+1; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-4; cno(eto,2) = no-3; cno(eto,3) = no; cno(eto,4) = no-1; levelso[eto] = let+1;eto++;
    cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] = let+1;eto++;
  }
}
//=============================================================================
void K_GENERATOR::insert2DPattern4(E_Int et, E_Int let, E_Int* tag,
                                   FldArrayF& f, FldArrayI& cn,
                                   FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso,
                                   E_Int& eto, E_Int& no)
{
  E_Int c1, c2, c3, c4;
  E_Int nfld = f.getNfld();
  c1 = cn(et,1)-1; c2 = cn(et,2)-1; c3 = cn(et,3)-1; c4 = cn(et,4)-1;
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
  no++;// no-11
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 1./3.*f(c1,nf)+2./3.*f(c2,nf);
  no++;// no-10
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
  no++;// no-9
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
  no++;// no-8
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
  no++;// no-7
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c3,nf);
  no++;// no-6
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c1,nf);
  no++;// no-5
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 1./2.*f(c4,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c1,nf);
  no++;// no-4
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 1./2.*f(c3,nf)+1./6.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c4,nf);
  no++;// no-3
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c2,nf);
  no++;// no-2
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c3,nf);
  no++;// no-1
  for (E_Int nf = 1; nf <= nfld; nf++)
    fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c4,nf);
  no++;// no  
  cno(eto,1) = c1+1; cno(eto,2) = no-11; cno(eto,3) = no-8; cno(eto,4) = no-9; levelso[eto] = let+1;eto++;
  cno(eto,1) = no-11; cno(eto,2) = no-10; cno(eto,3) = no-7; cno(eto,4) = no-8; levelso[eto] = let+1;eto++;
  cno(eto,1) = no-10; cno(eto,2) = c2+1; cno(eto,3) = no-6; cno(eto,4) = no-7; levelso[eto] = let+1;eto++;
  cno(eto,1) = no-9; cno(eto,2) = no-8; cno(eto,3) = no-4; cno(eto,4) = no-5; levelso[eto] = let+1;eto++;
  cno(eto,1) = no-8; cno(eto,2) = no-7; cno(eto,3) = no-3; cno(eto,4) = no-4; levelso[eto] = let+1;eto++;
  cno(eto,1) = no-7; cno(eto,2) = no-6; cno(eto,3) = no-2; cno(eto,4) = no-3; levelso[eto] = let+1;eto++;
  cno(eto,1) = no-5; cno(eto,2) = no-4; cno(eto,3) = no-1; cno(eto,4) = c4+1; levelso[eto] = let+1;eto++;
  cno(eto,1) = no-4; cno(eto,2) = no-3; cno(eto,3) = no; cno(eto,4) = no-1; levelso[eto] = let+1;eto++;
  cno(eto,1) = no-3; cno(eto,2) = no-2; cno(eto,3) = c3+1; cno(eto,4) = no; levelso[eto] = let+1;eto++;
}

//=============================================================================
// Determine le template a utiliser pour le cas ntags=2 et 3D
// Retourne 0 (A), 1(B), 2(C)
//=============================================================================
E_Int K_GENERATOR::find2Template(E_Int* tags)
{
  // Cas 2C
  if (tags[0] == 1 && tags[6] == 1) return 2;
  if (tags[1] == 1 && tags[7] == 1) return 2;
  if (tags[2] == 1 && tags[4] == 1) return 2;
  if (tags[3] == 1 && tags[5] == 1) return 2;

  // Cas 2A
  if (tags[0] == 1 && tags[1] == 1) return 0;
  if (tags[1] == 1 && tags[2] == 1) return 0;
  if (tags[2] == 1 && tags[3] == 1) return 0;
  if (tags[3] == 1 && tags[0] == 1) return 0;
  if (tags[4] == 1 && tags[5] == 1) return 0;
  if (tags[5] == 1 && tags[6] == 1) return 0;
  if (tags[6] == 1 && tags[7] == 1) return 0;
  if (tags[7] == 1 && tags[4] == 1) return 0;
  if (tags[0] == 1 && tags[4] == 1) return 0;
  if (tags[3] == 1 && tags[7] == 1) return 0;
  if (tags[1] == 1 && tags[5] == 1) return 0;
  if (tags[2] == 1 && tags[6] == 1) return 0;

  // Cas 2B
  return 1;
}

//=============================================================================
// Determine le template a utiliser pour le cas ntags=3 et 3D
// Retourne 0 (A), 1(B,C)
//=============================================================================
E_Int K_GENERATOR::find3Template(E_Int* tags)
{
  if (tags[0] == 1 && tags[1] == 1 && tags[2] == 1) return 0;
  if (tags[1] == 1 && tags[2] == 1 && tags[3] == 1) return 0;
  if (tags[2] == 1 && tags[3] == 1 && tags[0] == 1) return 0;
  if (tags[3] == 1 && tags[0] == 1 && tags[1] == 1) return 0;

  if (tags[4] == 1 && tags[5] == 1 && tags[6] == 1) return 0;
  if (tags[5] == 1 && tags[6] == 1 && tags[7] == 1) return 0;
  if (tags[6] == 1 && tags[7] == 1 && tags[4] == 1) return 0;
  if (tags[7] == 1 && tags[4] == 1 && tags[5] == 1) return 0;

  if (tags[0] == 1 && tags[3] == 1 && tags[7] == 1) return 0;
  if (tags[3] == 1 && tags[7] == 1 && tags[4] == 1) return 0;
  if (tags[7] == 1 && tags[4] == 1 && tags[0] == 1) return 0;
  if (tags[4] == 1 && tags[0] == 1 && tags[3] == 1) return 0;

  if (tags[1] == 1 && tags[2] == 1 && tags[6] == 1) return 0;
  if (tags[2] == 1 && tags[6] == 1 && tags[5] == 1) return 0;
  if (tags[6] == 1 && tags[5] == 1 && tags[1] == 1) return 0;
  if (tags[5] == 1 && tags[1] == 1 && tags[2] == 1) return 0;

  if (tags[0] == 1 && tags[1] == 1 && tags[5] == 1) return 0;
  if (tags[1] == 1 && tags[5] == 1 && tags[4] == 1) return 0;
  if (tags[5] == 1 && tags[4] == 1 && tags[0] == 1) return 0;
  if (tags[4] == 1 && tags[0] == 1 && tags[1] == 1) return 0;

  if (tags[3] == 1 && tags[2] == 1 && tags[6] == 1) return 0;
  if (tags[2] == 1 && tags[6] == 1 && tags[7] == 1) return 0;
  if (tags[6] == 1 && tags[7] == 1 && tags[3] == 1) return 0;
  if (tags[7] == 1 && tags[3] == 1 && tags[2] == 1) return 0;

  return 1;
}

//=============================================================================
// Determine le template a utiliser pour le cas ntags=4 et 3D
// Retourne 0 (A), 1(B,C,D,E,F)
//=============================================================================
E_Int K_GENERATOR::find4Template(E_Int* tags)
{
  if (tags[0] == 1 && tags[1] == 1 && tags[2] == 1 && tags[3] == 1) return 0;
  if (tags[4] == 1 && tags[5] == 1 && tags[6] == 1 && tags[7] == 1) return 0;
  if (tags[0] == 1 && tags[3] == 1 && tags[7] == 1 && tags[4] == 1) return 0;
  if (tags[1] == 1 && tags[2] == 1 && tags[6] == 1 && tags[5] == 1) return 0;
  if (tags[0] == 1 && tags[1] == 1 && tags[5] == 1 && tags[4] == 1) return 0;
  if (tags[3] == 1 && tags[2] == 1 && tags[6] == 1 && tags[7] == 1) return 0;
  return 1;
}
//=============================================================================
void K_GENERATOR::insert3DPattern1(E_Int et, E_Int let, E_Int* tag,
                                   FldArrayF& f, FldArrayI& cn,
                                   FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, 
                                   E_Int& eto, E_Int& no)
{
  E_Int c1, c2, c3, c4, c5, c6, c7, c8;
  c1 = cn(et,1)-1; c2 = cn(et,2)-1; c3 = cn(et,3)-1; c4 = cn(et,4)-1;
  c5 = cn(et,5)-1; c6 = cn(et,6)-1; c7 = cn(et,7)-1; c8 = cn(et,8)-1;
  E_Int nfld = f.getNfld();

  if (tag[0] == 1)
  {
    for (E_Int nf = 1; nf <= nfld; nf++)
    {
      fo(no, nf) = 2./3.*f(c1,nf)+1./3.*f(c2,nf);
      fo(no+1, nf) = 1./2.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf);
      fo(no+2, nf) = 2./3.*f(c1,nf)+1./3.*f(c4,nf);
      fo(no+3, nf) = 2./3.*f(c1,nf)+1./3.*f(c5,nf);
      fo(no+4, nf) = 1./2.*f(c1,nf)+1./6.*f(c5,nf)+1./6.*f(c6,nf)+1./6.*f(c2,nf);
      fo(no+5, nf) = 5./12.*f(c1,nf)+1./12.*f(c2,nf)+1./12.*f(c3,nf)+1./12.*f(c4,nf)+1./12.*f(c5,nf)+1./12.*f(c6,nf)+1./12.*f(c7,nf)+1./12.*f(c8,nf);
      fo(no+6, nf) = 1./2.*f(c1,nf)+1./6.*f(c4,nf)+1./6.*f(c5,nf)+1./6.*f(c8,nf);
    }
    cn(eto,1) = c1+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c2+1; cn(eto,6) = c3+1; cn(eto,7) = c7+1; cn(eto,8) = c6+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c3+1; cn(eto,6) = c4+1; cn(eto,7) = c8+1; cn(eto,8) = c7+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c5+1; cn(eto,6) = c6+1; cn(eto,7) = c7+1; cn(eto,8) = c8+1; levelso[eto] = -let; eto++;
    no += 7;
  }
  else if (tag[1] == 1)
  {
    for (E_Int nf = 1; nf <= nfld; nf++)
    {
      fo(no, nf) = 2./3.*f(c2,nf)+1./3.*f(c3,nf);
      fo(no+1, nf) = 1./2.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf)+1./6.*f(c1,nf);
      fo(no+2, nf) = 2./3.*f(c2,nf)+1./3.*f(c1,nf);
      fo(no+3, nf) = 2./3.*f(c2,nf)+1./3.*f(c6,nf);
      fo(no+4, nf) = 1./2.*f(c2,nf)+1./6.*f(c6,nf)+1./6.*f(c7,nf)+1./6.*f(c3,nf);
      fo(no+5, nf) = 5./12.*f(c2,nf)+1./12.*f(c1,nf)+1./12.*f(c3,nf)+1./12.*f(c4,nf)+1./12.*f(c5,nf)+1./12.*f(c6,nf)+1./12.*f(c7,nf)+1./12.*f(c8,nf);
      fo(no+6, nf) = 1./2.*f(c2,nf)+1./6.*f(c1,nf)+1./6.*f(c5,nf)+1./6.*f(c6,nf);
    }
    cn(eto,1) = c2+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c3+1; cn(eto,6) = c4+1; cn(eto,7) = c8+1; cn(eto,8) = c7+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c4+1; cn(eto,6) = c1+1; cn(eto,7) = c5+1; cn(eto,8) = c8+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c6+1; cn(eto,6) = c7+1; cn(eto,7) = c8+1; cn(eto,8) = c5+1; levelso[eto] = -let; eto++;
    no += 7;
  }
  else if (tag[2] == 1)
  {
    for (E_Int nf = 1; nf <= nfld; nf++)
    {
      fo(no, nf) = 2./3.*f(c3,nf)+1./3.*f(c4,nf);
      fo(no+1, nf) = 1./2.*f(c3,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c4,nf);
      fo(no+2, nf) = 2./3.*f(c3,nf)+1./3.*f(c2,nf);
      fo(no+3, nf) = 2./3.*f(c3,nf)+1./3.*f(c7,nf);
      fo(no+4, nf) = 1./2.*f(c3,nf)+1./6.*f(c4,nf)+1./6.*f(c7,nf)+1./6.*f(c8,nf);
      fo(no+5, nf) = 5./12.*f(c3,nf)+1./12.*f(c1,nf)+1./12.*f(c2,nf)+1./12.*f(c4,nf)+1./12.*f(c5,nf)+1./12.*f(c6,nf)+1./12.*f(c7,nf)+1./12.*f(c8,nf);
      fo(no+6, nf) = 1./2.*f(c3,nf)+1./6.*f(c2,nf)+1./6.*f(c6,nf)+1./6.*f(c7,nf);
    }
    cn(eto,1) = c3+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c4+1; cn(eto,6) = c1+1; cn(eto,7) = c5+1; cn(eto,8) = c8+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c1+1; cn(eto,6) = c2+1; cn(eto,7) = c6+1; cn(eto,8) = c5+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c7+1; cn(eto,6) = c8+1; cn(eto,7) = c5+1; cn(eto,8) = c6+1; levelso[eto] = -let; eto++;
    no += 7;
  }
  else if (tag[3] == 1)
  {
    for (E_Int nf = 1; nf <= nfld; nf++)
    {
      fo(no, nf) = 2./3.*f(c4,nf)+1./3.*f(c1,nf);
      fo(no+1, nf) = 1./2.*f(c4,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf);
      fo(no+2, nf) = 2./3.*f(c4,nf)+1./3.*f(c3,nf);
      fo(no+3, nf) = 2./3.*f(c4,nf)+1./3.*f(c8,nf);
      fo(no+4, nf) = 1./2.*f(c4,nf)+1./6.*f(c1,nf)+1./6.*f(c5,nf)+1./6.*f(c8,nf);
      fo(no+5, nf) = 5./12.*f(c4,nf)+1./12.*f(c1,nf)+1./12.*f(c2,nf)+1./12.*f(c3,nf)+1./12.*f(c5,nf)+1./12.*f(c6,nf)+1./12.*f(c7,nf)+1./12.*f(c8,nf);
      fo(no+6, nf) = 1./2.*f(c4,nf)+1./6.*f(c3,nf)+1./6.*f(c7,nf)+1./6.*f(c8,nf);
    }
    cn(eto,1) = c4+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c1+1; cn(eto,6) = c2+1; cn(eto,7) = c6+1; cn(eto,8) = c5+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c2+1; cn(eto,6) = c3+1; cn(eto,7) = c7+1; cn(eto,8) = c6+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c8+1; cn(eto,6) = c5+1; cn(eto,7) = c6+1; cn(eto,8) = c7+1; levelso[eto] = -let; eto++;
    no += 7;
  }
  else if (tag[4] == 1)
  {
    for (E_Int nf = 1; nf <= nfld; nf++)
    {
      fo(no, nf) = 2./3.*f(c5,nf)+1./3.*f(c6,nf);
      fo(no+1, nf) = 1./2.*f(c5,nf)+1./6.*f(c6,nf)+1./6.*f(c7,nf)+1./6.*f(c8,nf);
      fo(no+2, nf) = 2./3.*f(c5,nf)+1./3.*f(c8,nf);
      fo(no+3, nf) = 2./3.*f(c5,nf)+1./3.*f(c1,nf);
      fo(no+4, nf) = 1./2.*f(c5,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c6,nf);
      fo(no+5, nf) = 5./12.*f(c5,nf)+1./12.*f(c1,nf)+1./12.*f(c2,nf)+1./12.*f(c3,nf)+1./12.*f(c4,nf)+1./12.*f(c6,nf)+1./12.*f(c7,nf)+1./12.*f(c8,nf);
      fo(no+6, nf) = 1./2.*f(c5,nf)+1./6.*f(c8,nf)+1./6.*f(c4,nf)+1./6.*f(c1,nf);
    }
    cn(eto,1) = c5+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c6+1; cn(eto,6) = c7+1; cn(eto,7) = c3+1; cn(eto,8) = c2+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c7+1; cn(eto,6) = c8+1; cn(eto,7) = c4+1; cn(eto,8) = c3+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c1+1; cn(eto,6) = c2+1; cn(eto,7) = c3+1; cn(eto,8) = c4+1; levelso[eto] = -let; eto++;
    no += 7;
  }
  else if (tag[5] == 1)
  {
    for (E_Int nf = 1; nf <= nfld; nf++)
    {
      fo(no, nf) = 2./3.*f(c6,nf)+1./3.*f(c7,nf);
      fo(no+1, nf) = 1./2.*f(c6,nf)+1./6.*f(c5,nf)+1./6.*f(c7,nf)+1./6.*f(c8,nf);
      fo(no+2, nf) = 2./3.*f(c6,nf)+1./3.*f(c5,nf);
      fo(no+3, nf) = 2./3.*f(c6,nf)+1./3.*f(c2,nf);
      fo(no+4, nf) = 1./2.*f(c6,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c7,nf);
      fo(no+5, nf) = 5./12.*f(c6,nf)+1./12.*f(c1,nf)+1./12.*f(c2,nf)+1./12.*f(c3,nf)+1./12.*f(c4,nf)+1./12.*f(c5,nf)+1./12.*f(c7,nf)+1./12.*f(c8,nf);
      fo(no+6, nf) = 1./2.*f(c6,nf)+1./6.*f(c1,nf)+1./6.*f(c2,nf)+1./6.*f(c5,nf);
    }
    cn(eto,1) = c6+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c7+1; cn(eto,6) = c8+1; cn(eto,7) = c4+1; cn(eto,8) = c3+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c8+1; cn(eto,6) = c5+1; cn(eto,7) = c1+1; cn(eto,8) = c4+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c2+1; cn(eto,6) = c3+1; cn(eto,7) = c4+1; cn(eto,8) = c1+1; levelso[eto] = -let; eto++;
    no += 7;
  }
  else if (tag[6] == 1)
  {
    for (E_Int nf = 1; nf <= nfld; nf++)
    {
      fo(no, nf) = 2./3.*f(c7,nf)+1./3.*f(c8,nf);
      fo(no+1, nf) = 1./2.*f(c7,nf)+1./6.*f(c5,nf)+1./6.*f(c6,nf)+1./6.*f(c8,nf);
      fo(no+2, nf) = 2./3.*f(c7,nf)+1./3.*f(c6,nf);
      fo(no+3, nf) = 2./3.*f(c7,nf)+1./3.*f(c3,nf);
      fo(no+4, nf) = 1./2.*f(c7,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf)+1./6.*f(c8,nf);
      fo(no+5, nf) = 5./12.*f(c7,nf)+1./12.*f(c1,nf)+1./12.*f(c2,nf)+1./12.*f(c3,nf)+1./12.*f(c4,nf)+1./12.*f(c5,nf)+1./12.*f(c6,nf)+1./12.*f(c8,nf);
      fo(no+6, nf) = 1./2.*f(c7,nf)+1./6.*f(c2,nf)+1./6.*f(c3,nf)+1./6.*f(c6,nf);
    }
    cn(eto,1) = c7+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c8+1; cn(eto,6) = c5+1; cn(eto,7) = c1+1; cn(eto,8) = c4+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c5+1; cn(eto,6) = c6+1; cn(eto,7) = c2+1; cn(eto,8) = c1+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c3+1; cn(eto,6) = c4+1; cn(eto,7) = c1+1; cn(eto,8) = c2+1; levelso[eto] = -let; eto++;
    no += 7;
  }
  else if (tag[7] == 1)
  {
    for (E_Int nf = 1; nf <= nfld; nf++)
    {
      fo(no, nf) = 2./3.*f(c8,nf)+1./3.*f(c5,nf);
      fo(no+1, nf) = 1./2.*f(c8,nf)+1./6.*f(c5,nf)+1./6.*f(c6,nf)+1./6.*f(c7,nf);
      fo(no+2, nf) = 2./3.*f(c8,nf)+1./3.*f(c7,nf);
      fo(no+3, nf) = 2./3.*f(c8,nf)+1./3.*f(c4,nf);
      fo(no+4, nf) = 1./2.*f(c8,nf)+1./6.*f(c4,nf)+1./6.*f(c1,nf)+1./6.*f(c5,nf);
      fo(no+5, nf) = 5./12.*f(c8,nf)+1./12.*f(c1,nf)+1./12.*f(c2,nf)+1./12.*f(c3,nf)+1./12.*f(c4,nf)+1./12.*f(c5,nf)+1./12.*f(c6,nf)+1./12.*f(c7,nf);
      fo(no+6, nf) = 1./2.*f(c8,nf)+1./6.*f(c3,nf)+1./6.*f(c4,nf)+1./6.*f(c7,nf);
    }
    cn(eto,1) = c8+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c5+1; cn(eto,6) = c6+1; cn(eto,7) = c2+1; cn(eto,8) = c1+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c6+1; cn(eto,6) = c7+1; cn(eto,7) = c3+1; cn(eto,8) = c2+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c4+1; cn(eto,6) = c1+1; cn(eto,7) = c2+1; cn(eto,8) = c3+1; levelso[eto] = -let; eto++;
    no += 7;
  }

}
//=============================================================================
void K_GENERATOR::insert3DPattern3(E_Int et, E_Int let, E_Int* tag,
                                   FldArrayF& f, FldArrayI& cn,
                                   FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, 
                                   E_Int& eto, E_Int& no)
{
  E_Int c1, c2, c3, c4, c5, c6, c7, c8;
  c1 = cn(et,1)-1; c2 = cn(et,2)-1; c3 = cn(et,3)-1; c4 = cn(et,4)-1;
  c5 = cn(et,5)-1; c6 = cn(et,6)-1; c7 = cn(et,7)-1; c8 = cn(et,8)-1;
  E_Int nfld = f.getNfld();
  E_Float a1, a2, a3, b1, b2, b3;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;

  if (tag[0] == 1 && tag[1] == 1 && tag[2] == 1 && tag[3] == 1)
  {
    E_Int inc = no;
    for (E_Int k = 0; k < 2; k++)
      for (E_Int j = 0; j < 4; j++)
        for (E_Int i = 0; i < 4; i++)
        {
          a1 = i*1./3.; a2 = j*1./3.; a3 = k*1./3.;
          b1 = 1.-a1; b2 = 1.-a2; b3 = 1.-a3;
          for (E_Int nf = 1; nf <= nfld; nf++)
            fo(no, nf) = b1*b2*b3*f(c1,nf)+a1*b2*b3*f(c2,nf)+a1*a2*b3*f(c3,nf)+b1*a2*b3*f(c4,nf)+
              b1*b2*a3*f(c5,nf)+a1*b2*a3*f(c6,nf)+a1*a2*a3*f(c7,nf)+b1*a2*a3*f(c8,nf);
          no++;
        }
    for (E_Int k = 2; k < 3; k++)
      for (E_Int j = 0; j < 4; j++)
        for (E_Int i = 0; i < 4; i++)
        {
          if ((i != 0 || j != 0) && (i != 0 || j != 3) && (i != 3 || j != 0) && (i != 3 || j != 3))
          {
            a1 = i*1./3.; a2 = j*1./3.; a3 = k*1./3.;
            b1 = 1.-a1; b2 = 1.-a2; b3 = 1.-a3;
            for (E_Int nf = 1; nf <= nfld; nf++)
              fo(no, nf) = b1*b2*b3*f(c1,nf)+a1*b2*b3*f(c2,nf)+a1*a2*b3*f(c3,nf)+b1*a2*b3*f(c4,nf)+
                b1*b2*a3*f(c5,nf)+a1*b2*a3*f(c6,nf)+a1*a2*a3*f(c7,nf)+b1*a2*a3*f(c8,nf);
            no++;
          }
        }
    for (E_Int nf = 1; nf <= nfld; nf++) fo(no, nf) = f(c5, nf); 
    no++;
    for (E_Int nf = 1; nf <= nfld; nf++) fo(no, nf) = f(c6, nf); 
    no++;
    for (E_Int nf = 1; nf <= nfld; nf++) fo(no, nf) = f(c7, nf); 
    no++;
    for (E_Int nf = 1; nf <= nfld; nf++) fo(no, nf) = f(c8, nf); 
    no++;
    for (E_Int kc = 0; kc < 1; kc++)
      for (E_Int jc = 0; jc < 3; jc++)
        for (E_Int ic = 0; ic < 3; ic++)
        {
          ind1 = ic+jc*4+kc*16+inc; ind2 = ind1+1; ind3 = ind2 + 4; ind4 = ind1+4; 
          ind5 = ind1+16; ind6 = ind5+1; ind7 = ind6+4; ind8 = ind5+4; 
          cno(eto,1) = ind1+1; cno(eto,2) = ind2+1; cno(eto,3) = ind3+1; cno(eto,4) = ind4+1;
          cno(eto,5) = ind5+1; cno(eto,6) = ind6+1; cno(eto,7) = ind7+1; cno(eto,8) = ind8+1;
          levelso[eto] = let+1; eto++;
        }
    

    cn(eto,1) = c1+1; cn(eto,2) = no+1; cn(eto,3) = no+2; cn(eto,4) = no+3; cn(eto,5) = no+4; cn(eto,6) = no+5; cn(eto,7) = no+6; cn(eto,8) = no+7; levelso[eto] =let+1; eto++;
    cn(eto,1) = no+1; cn(eto,2) = no+2; cn(eto,3) = no+6; cn(eto,4) = no+5; cn(eto,5) = c2+1; cn(eto,6) = c3+1; cn(eto,7) = c7+1; cn(eto,8) = c6+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+2; cn(eto,2) = no+3; cn(eto,3) = no+7; cn(eto,4) = no+6; cn(eto,5) = c3+1; cn(eto,6) = c4+1; cn(eto,7) = c8+1; cn(eto,8) = c7+1; levelso[eto] = -let; eto++;
    cn(eto,1) = no+4; cn(eto,2) = no+5; cn(eto,3) = no+6; cn(eto,4) = no+7; cn(eto,5) = c5+1; cn(eto,6) = c6+1; cn(eto,7) = c7+1; cn(eto,8) = c8+1; levelso[eto] = -let; eto++;
    no += 7;
  }
 
}
//=============================================================================
void K_GENERATOR::insert3DPattern4(E_Int et, E_Int let, E_Int* tag,
                                   FldArrayF& f, FldArrayI& cn,
                                   FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, 
                                   E_Int& eto, E_Int& no)
{
  // les coins sont dupliques
  E_Int c1, c2, c3, c4, c5, c6, c7, c8;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
  c1 = cn(et,1)-1; c2 = cn(et,2)-1; c3 = cn(et,3)-1; c4 = cn(et,4)-1;
  c5 = cn(et,5)-1; c6 = cn(et,6)-1; c7 = cn(et,7)-1; c8 = cn(et,8)-1;
  E_Float a1, a2, a3, b1, b2, b3;
  E_Int nfld = f.getNfld();
  E_Int inc = no;
  for (E_Int k = 0; k < 4; k++)
    for (E_Int j = 0; j < 4; j++)
      for (E_Int i = 0; i < 4; i++)
      {
        a1 = i*1./3.; a2 = j*1./3.; a3 = k*1./3.;
        b1 = 1.-a1; b2 = 1.-a2; b3 = 1.-a3;
         for (E_Int nf = 1; nf <= nfld; nf++)
           fo(no, nf) = b1*b2*b3*f(c1,nf)+a1*b2*b3*f(c2,nf)+a1*a2*b3*f(c3,nf)+b1*a2*b3*f(c4,nf)+
             b1*b2*a3*f(c5,nf)+a1*b2*a3*f(c6,nf)+a1*a2*a3*f(c7,nf)+b1*a2*a3*f(c8,nf);
         no++;
      }
  E_Int n = 4; E_Int n2 = 16;
  for (E_Int kc = 0; kc < n-1; kc++)
    for (E_Int jc = 0; jc < n-1; jc++)
      for (E_Int ic = 0; ic < n-1; ic++)
      {
        ind1 = ic+jc*n+kc*n2+inc; ind2 = ind1+1; ind3 = ind2 + n; ind4 = ind1+n; 
        ind5 = ind1+n2; ind6 = ind5+1; ind7 = ind6+n; ind8 = ind5+n; 
        cno(eto,1) = ind1+1; cno(eto,2) = ind2+1; cno(eto,3) = ind3+1; cno(eto,4) = ind4+1;
        cno(eto,5) = ind5+1; cno(eto,6) = ind6+1; cno(eto,7) = ind7+1; cno(eto,8) = ind8+1;
        levelso[eto] = let+1; eto++;
      }
}
