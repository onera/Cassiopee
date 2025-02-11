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

# include "transform.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

#define NORM2(vx,vy,vz) vx*vx+vy*vy+vz*vz
#define SCAL(v1x,v1y,v1z,v2x,v2y,v2z) v1x*v2x+v1y*v2y+v1z*v2z
#define VECTX(v1x,v1y,v1z,v2x,v2y,v2z) v1x*v2z-v1z*v2x
#define VECTY(v1x,v1y,v1z,v2x,v2y,v2z) v1z*v2y-v1y*v2z
#define VECTZ(v1x,v1y,v1z,v2x,v2y,v2z) v1x*v2y-v1y*v2x

// ============================================================================
/* Deform mesh by moving surface of a given vector, with collapse en expand 
 */
// ============================================================================
PyObject* K_TRANSFORM::deform2(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* normal;

  if (!PYPARSETUPLE_(args, OO_, &array, &normal))
    return NULL;
  
  // Check array
  E_Int im1, jm1, km1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int im2, jm2, km2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  E_Int res1 = K_ARRAY::getFromArray(array, varString1, f1, 
                                     im1, jm1, km1, cn1, eltType1, true);
  
  E_Int res2 = K_ARRAY::getFromArray(normal, varString2, f2, 
                                     im2, jm2, km2, cn2, eltType2, true);
  
  // Vecteur et array valides (structure ou non structure) ?
  if (res1 == -1)
  {
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_TypeError,
                    "deform: 1st argument is invalid.");
    return NULL;
  }
  
  // Only for NGONS
  if (res1 != 2 || strcmp(eltType1, "NGON") != 0)
  {
    RELEASESHAREDB(res1, normal,f1,cn1);
    PyErr_SetString(PyExc_TypeError,
                    "deform: 1st argument must be a NGON.");
    return NULL;
  }
  
  if (res2 == -1)
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    PyErr_SetString(PyExc_TypeError,
                    "deform: 2nd argument is invalid.");
    return NULL;
  }

  if (res1 != res2) 
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_TypeError,
                    "deform: 1st and 2nd argument must be both structured or unstructured.");
    return NULL;
  }
  // Presence des coordonnees ?
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_ValueError,
                    "deform: coordinates not found in 1st argument.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  E_Int npts = f1->getSize();
  E_Int nFld = f1->getNfld(); 
  
  // Vecteur normal de dimension 3 ?
  E_Int dim = f2->getNfld();
  if (dim != 3)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    RELEASESHAREDB(res2, normal, f2, cn2);
    PyErr_SetString(PyExc_ValueError,
                    "deform: 2nd argument must have 3 variables defining the vector.");
    return NULL;
  }
  
  // Vecteur et array de la meme taille ?
  if (npts != f2->getSize())
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_ValueError, 
                    "deform: sizes of 1st and 2nd arguments are not equal.");
    return NULL;
  }

  // Connectivite des faces incidentes a un noeud
  E_Int* cNG = cn1->begin();
  E_Int nfaces = cNG[0];
  E_Int nelts = cNG[2+cNG[1]];
  E_Int* cNGe = cNG + 2 + cNG[1];
  //printf("input surface has " SF_D_ " faces and " SF_D_ " elements\n", nfaces, nelts);
  vector< vector<E_Int> > cVF(npts); K_CONNECT::connectNG2VF(*cn1, cVF);
  FldArrayI cFE; K_CONNECT::connectNG2FE(*cn1, cFE);
  vector< vector<E_Int> > cEV(nelts); K_CONNECT::connectNG2EV(*cn1,cEV);

  FldArrayI posFaces(nfaces); K_CONNECT::getPosFaces(*cn1, posFaces);
  E_Int* posFacesp = posFaces.begin();
  FldArrayI posElts(nelts); K_CONNECT::getPosElts(*cn1, posElts);
  //E_Int* posEltsp = posElts.begin();

  // tableau pour stocker la face triangle par noeuds du maillage initial
  vector<E_Int>** ft = new vector<E_Int>* [npts];
  for (E_Int i = 0; i < npts; i++) ft[i] = NULL;

  E_Int indf, e1, e2, current, posindf, posindn, n1, n, found, eprev;
  E_Int ind1, ind2, nbn, ind, nf, jloc, node, indft, if1, if2;
  E_Int e3, e4, e, elt;
  E_Int* pt; E_Int* pt1; E_Int* pt2; E_Int* ptf1; E_Int* ptf2;
  E_Float x1, x2, y1, y2, z1, z2, l1, l2, d, xp1, xp2, yp1, yp2, zp1, zp2;
  E_Float v1x,v1y,v1z,v2x,v2y,v2z,d1x,d1y,d1z,r1x,r1y,r1z,r2x,r2y,r2z,s;
  bool boundary;
  E_Float* f1x = f1->begin(posx);
  E_Float* f1y = f1->begin(posy);
  E_Float* f1z = f1->begin(posz);
  E_Float* f2x = f2->begin(1);
  E_Float* f2y = f2->begin(2);
  E_Float* f2z = f2->begin(3);

  // 1. Tag input surface (-1: collapse, 0: normal, +1: expand)
  FldArrayI tag(npts); tag.setAllValuesAtNull();
  E_Int* tagp = tag.begin();
  //tagp[10] = -1; tagp[11] = -1; tagp[2]= 0;
  // Essai pour trouver tag a partir des longeurs des faces
  E_Int collapsed = 0; E_Int expanded = 0;
  for (E_Int i = 0; i < nfaces; i++)
  {
   pt = cNG+posFaces[i];
   ind1 = pt[1]-1; ind2 = pt[2]-1; 
   x1 = f1x[ind1]; y1 = f1y[ind1]; z1 = f1z[ind1];
   x2 = f1x[ind2]; y2 = f1y[ind2]; z2 = f1z[ind2];
   xp1 = x1 + f2x[ind1]; yp1 = y1 + f2y[ind1]; zp1 = z1 + f2z[ind1];
   xp2 = x2 + f2x[ind2]; yp2 = y2 + f2y[ind2]; zp2 = z2 + f2z[ind2];
   l1 = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
   l2 = (xp2-xp1)*(xp2-xp1)+(yp2-yp1)*(yp2-yp1)+(zp2-zp1)*(zp2-zp1);
   l2 = std::max(l2, 1.e-12);
   // Critere de rapport de longeurs d'aretes
   d = l1/l2;
   d = sqrt(d);
   if (d > 1.5) { tagp[ind1] = -1; tagp[ind2] = -1; collapsed++; } // collapse
   else if (d < -0.5) {tagp[ind1] = +1; tagp[ind2] = +1; expanded++; } // expand
   // Critere de retournement d'aretes
   v1x = x2-x1; v1y = y2-y1; v1z = z2-z1;
   v2x = xp2-xp1; v2y = yp2-yp1; v2z = zp2-zp1;
   d1x =  xp1-x1; d1y = yp1-y1; d1z = zp1-z1;
   r1x = VECTX(v1x,v1y,v1z,d1x,d1y,d1z);
   r1y = VECTY(v1x,v1y,v1z,d1x,d1y,d1z);
   r1z = VECTZ(v1x,v1y,v1z,d1x,d1y,d1z);
   r2x = VECTX(v2x,v2y,v2z,d1x,d1y,d1z);
   r2y = VECTY(v2x,v2y,v2z,d1x,d1y,d1z);
   r2z = VECTZ(v2x,v2y,v2z,d1x,d1y,d1z);
   s = SCAL(r1x,r1y,r1z,r2x,r2y,r2z);
   if (s < 0)  { tagp[ind1] = -1; tagp[ind2] = -1; collapsed++; } // collapse
  }
  printf("Expanded=" SF_D_ ", collapsed=" SF_D_ "\n", expanded, collapsed);

  // Reperage des elements touches par l'expand (qui possede un noeud 
  // tagge pour l'expand)
  // em contient ne nbre de noeuds a expander dans l'element
  FldArrayI em(nelts); em.setAllValuesAtNull();
  for (E_Int i = 0; i < nelts; i++)
  {
    vector<E_Int>& n = cEV[i];
    nbn = n.size();
    for (E_Int j = 0; j < nbn; j++)
    {
      ind = n[j];
      if (tagp[ind-1] == 1)
      {
        em[i] += 1;
      }
    }
  }
  //printf("======================\n");
  //for (E_Int i = 0; i < nelts; i++) printf("em[" SF_D_ "] = " SF_D_ "\n", i, em[i]);
  //printf("======================\n");

  // Compte le nbre de pts, nbre de faces, d'elements en plus
  // par les triangles ajoutes par expand
  E_Int nptsPlus = 0; E_Int nfacesPlus = 0; E_Int neltsPlus = 0;
  E_Int sizeNGonPlus = 0; E_Int sizeNFacePlus = 0;
  for (E_Int i = 0; i < npts; i++)
  {
    if (tagp[i] == 1) // expand
    {
      vector<E_Int>& f = cVF[i]; // faces incidente a i
      E_Int nbf = f.size(); 
      boundary = false;
      for (E_Int k = 0; k < nbf; k++) // frontiere?
      {
        indf = f[k];
        e1 = cFE(indf-1, 1); // premier element de la face
        e2 = cFE(indf-1, 2);
        if (e1 == 0 || e2 == 0) { boundary = true; break; }
      }

      nptsPlus += nbf;
      if (boundary == false) 
      { 
        nfacesPlus += nbf+nbf; sizeNGonPlus += 2*nbf*3; 
        neltsPlus += nbf; sizeNFacePlus += nbf*4;
        //sizeNFacePlus += neltsPlus*4 + neltsPlus;
      }
      else 
      { 
        nfacesPlus += 2*nbf-1; sizeNGonPlus += (2*nbf-1)*3;
        neltsPlus += nbf-1; sizeNFacePlus += (nbf-1)*4;
        //sizeNFacePlus += neltsPlus*4 + neltsPlus;
      }
    }
  }
  
  // Compte le nbre d'elements en moins par supression em puis rajout a la fin
  E_Int sizeNFaceMoins = 0; // raccourci la connectivite precedente
  for (E_Int i = 0; i < nelts; i++)
  {
    if (em[i] > 0) 
    { 
      // suppression dans la connectivite precedente
      neltsPlus += -1; pt = cNG+posElts[i];
      sizeNFacePlus += -pt[0]-1;
      sizeNFaceMoins += -pt[0]-1;
      // ajout a la fin de la nouvelle connectivite
      neltsPlus += 1;
      sizeNFacePlus += pt[0]+1+em[i];
    }
  }

  //printf("on cree " SF_D_ " pts, " SF_D_ " faces, " SF_D_ " elts\n", nptsPlus, nfacesPlus, neltsPlus);
  //printf("sizeNGonPlus=" SF_D_ " sizeNFacePluis=" SF_D_ " sizeNFaceMoins=" SF_D_ "\n", sizeNGonPlus, sizeNFacePlus, sizeNFaceMoins);
  // Construit l'array resultat et l'initialise par copie
  PyObject* tpl;
  E_Int sizeConnect = 4+cNG[1]+cNGe[1]+sizeNGonPlus+sizeNFacePlus;
  E_Int nptsNew = npts+nptsPlus;
  tpl = K_ARRAY::buildArray(nFld, varString1, nptsNew, sizeConnect,-1,eltType1,false,sizeConnect);      
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);

  // copie des faces (identiques a l'input)
  cnnp[0] = nfaces+nfacesPlus; // nbre de faces
  cnnp[1] = cNG[1]+sizeNGonPlus;
  E_Int* cne = cnnp+cnnp[1]+2;
  cne[0] = nelts+neltsPlus;
  cne[1] = cNGe[1]+sizeNFacePlus;
  for (E_Int i = 0; i < cNG[1]; i++) cnnp[i+2] = cNG[i+2];

  // Recopie avec suppression des elements em
  pt1 = cNGe+2; pt2 = cne+2;
  for (E_Int i = 0; i < nelts; i++)
  {
    n = pt1[0];
    if (em[i] == 0)
    {
      pt2[0] = n;
      for (E_Int j = 0; j < n; j++) pt2[j+1] = pt1[j+1];
      pt2 += n+1;
    }
    pt1 += n+1;
  }

  // Init des coordonnees
  E_Float* newFp = K_ARRAY::getFieldPtr(tpl);
  E_Float* fx = newFp;
  E_Float* fy = newFp+nptsNew;
  E_Float* fz = newFp+2*nptsNew;
  //#pragma omp parallel for default(shared)
  for (E_Int i = 0;  i < npts; i++) 
  {   
    fx[i] = f1x[i] + f2x[i];
    fy[i] = f1y[i] + f2y[i];
    fz[i] = f1z[i] + f2z[i];
  }

  printf("input surface has " SF_D_ " points\n", npts);

  // 4. Collapse some nodes
  E_Float fxm = 0.; E_Float fym = 0.; E_Float fzm = 0.;
  E_Int collapse;
  for (E_Int i = 0; i < nfaces; i++)
  {
    pt = cNG+posFaces[i];
    ind1 = pt[1]-1; ind2 = pt[2]-1;
    if (tagp[ind1] == -1 && tagp[ind2] ==  -1)
    {
      collapse = 0; // collapse -1:gauche, 0: milieu, 1:droite, 2: rien
      // Face frontiere ?
      e1 = cFE(i, 1); e2 = cFE(i, 2);
      if (e1 == 0 || e2 == 0) collapse = 2; // no collapse

      // node frontiere ?
      vector<E_Int>& f = cVF[ind1];
      E_Int nbf = f.size();
      for (E_Int k = 0; k < nbf; k++)
      {
        indf = f[k];
        e1 = cFE(indf-1, 1); // premier element de la face
        e2 = cFE(indf-1, 2);
        if (e1 == 0 || e2 == 0) { collapse = -1; break; }
      }
      vector<E_Int>& f2 = cVF[ind2];
      nbf = f2.size();
      for (E_Int k = 0; k < nbf; k++)
      {
        indf = f2[k];
        e1 = cFE(indf-1, 1); // premier element de la face
        e2 = cFE(indf-1, 2);
        if (e1 == 0 || e2 == 0) { collapse = +1; break; }
      }
      if (collapse == 0)
      {
        fxm = 0.5*(fx[ind1]+fx[ind2]);
        fym = 0.5*(fy[ind1]+fy[ind2]);
        fzm = 0.5*(fz[ind1]+fz[ind2]);
        fx[ind1] = fxm; fy[ind1] = fym; fz[ind1] = fzm;
        fx[ind2] = fxm; fy[ind2] = fym; fz[ind2] = fzm;
      }
      else if (collapse == -1)
      {
        fx[ind2] = fx[ind1]; fy[ind2] = fy[ind1]; fz[ind2] = fz[ind1];
      }
      else if (collapse == +1)
      {
        fx[ind1] = fx[ind2]; fy[ind1] = fy[ind2]; fz[ind1] = fz[ind2];
      }
      tagp[ind1] = 0; tagp[ind2] = 0; // forced
      //printf("" SF_D_ " collapse " SF_D_ "\n", i, collapse);
    }
  }

  // 2. Expand input surface (coord+faces)
  E_Int currNpts = npts;
  E_Int* ptFace = cnnp+2+cNG[1];
  E_Int currFace = nfaces;
  E_Int* ptElem = cnnp+4+cnnp[1]+cNGe[1]+sizeNFaceMoins;
  for (E_Int i = 0; i < npts; i++)
  {
    if (tagp[i] == 1)
    {
      vector<E_Int>& f = cVF[i]; // faces incidente a i
      E_Int nbf = f.size(); // nbre de faces incidentes

      // Construction de la liste des faces dans le bon ordre
      vector<E_Int> ordered(nbf);
      // Premiere face incidente
      current = 0; // face courante
      boundary = false;
      indf = f[0];
      for (E_Int k = 0; k < nbf; k++)
      {
        indf = f[k];
        e1 = cFE(indf-1, 1); // premier element de la face
        e2 = cFE(indf-1, 2);
        if (e1 == 0 || e2 == 0) { boundary = true; break; }
      }
      ordered[0] = indf;
      eprev = -1;

      while (current < nbf)
      {
        //printf("current face " SF_D_ "\n", indf);
        e1 = cFE(indf-1, 1); // premier element de la face
        e2 = cFE(indf-1, 2);
        if (e1 == eprev && e2 != 0) e1 = e2;
        if (e1 == 0) e1 = e2;
        
        //printf("current element: " SF_D_ "\n", e1);
        pt = cNG+posElts[e1-1]; // faces de e1
        n1 = pt[0]; // nbre de faces de e1
        //printf("nbre de face de l element " SF_D_ "\n", n1);
        
        posindf = 0; posindn = 0; found = 0;
        for (E_Int k = 0; k < n1; k++) // pour toutes faces de e1
        {
          if (pt[k+1] == indf) // position de indf dans e1
            posindf = k;
          else
          {
            for (E_Int j = 0; j < nbf; j++) // pour toutes faces incidentes
            {
              if (pt[k+1] == f[j])
              {
                // j est la suivante
                posindn = k;
                found = pt[k+1];
              }
            }
          }
        }
        if (posindn > posindf) ordered[current+1] = found;
        else ordered[current+1] = -found;

        current += 1;
        indf = found;
        eprev = e1;
      }

      // Construction des pts sur les facettes incidentes
      for (E_Int k = 0; k < nbf; k++)
      {
        indf = abs(ordered[k]);
        pt = cNG+posFacesp[indf-1];
        ind1 = pt[1]-1; ind2 = pt[2]-1;
        //printf("face " SF_D_ " : " SF_D2_ "\n",indf,ind1,ind2);
        fx[currNpts+k] = 0.5*(fx[ind1]+fx[ind2]);
        fy[currNpts+k] = 0.5*(fy[ind1]+fy[ind2]);
        fz[currNpts+k] = 0.5*(fz[ind1]+fz[ind2]);
        //fx[currNpts+k] = f1x[i];
        //fy[currNpts+k] = f1y[i];
        //fz[currNpts+k] = f1z[i];

        //f2nx[currNpts+k] = 0.5*(f2x[ind1]+f2x[ind2]);
        //f2ny[currNpts+k] = 0.5*(f2y[ind1]+f2y[ind2]);
        //f2nz[currNpts+k] = 0.5*(f2z[ind1]+f2z[ind2]);

        ptFace[0] = 2;
        ptFace[1] = i+1; ptFace[2] = currNpts+k+1;
        ptFace += 3;
      }

      // Construction des faces opposees sur les facettes incidentes
      for (E_Int k = 0; k < nbf-1; k++)
      {
        ptFace[0] = 2;
        ptFace[1] = currNpts+k+1; ptFace[2] = currNpts+k+2;
        ptFace += 3;
      }

      if (boundary == false)
      {
        ptFace[0] = 2;
        ptFace[1] = currNpts+nbf; ptFace[2] = currNpts+1;
        ptFace += 3;
      }

      // Modification des facette initiales (racourcies)
      for (E_Int k = 0; k < nbf; k++)
      {
        indf = abs(ordered[k]);
        pt = cnnp+posFacesp[indf-1];
        if (pt[1] == i+1) pt[1] = currNpts+k+1;
        else if (pt[2] == i+1) pt[2] = currNpts+k+1;
        else printf("error\n");       
      }

      currNpts += nbf;

      // Construction des elements triangulaires
      for (E_Int k = 0; k < nbf-1; k++)
      {
        ptElem[0] = 3;
        if (ordered[k] < 0)
        {
          ptElem[1] = currFace+k+1;
          ptElem[2] = currFace+nbf+k+1;
          ptElem[3] = currFace+k+2;
        }
        else
        {
          ptElem[1] = currFace+k+2;
          ptElem[2] = currFace+nbf+k+1;
          ptElem[3] = currFace+k+1;
        }
        ptElem += 4;

        if (ft[i] == NULL) ft[i] = new vector<E_Int>;
        if1 = abs(ordered[k]);
        if2 = abs(ordered[k+1]);
        e1 = cFE(if1-1,1);
        e2 = cFE(if1-1,2);
        e3 = cFE(if2-1,1);
        e4 = cFE(if2-1,2);
        e = 0;
        if (e1 == e3 && e1 != 0) e = e1;
        else if (e1 == e4 && e1 != 0) e = e1;
        else if (e2 == e3 && e2 != 0) e = e2;
        else if (e2 == e4 && e2 != 0) e = e2;
        ft[i]->push_back(currFace+nbf+k+1); // stockage de la face opposee pour le noeud i
        ft[i]->push_back(e); 
      }

      if (boundary == false)
      {
        ptElem[0] = 3;
        if (ordered[nbf-1] < 0)
        {
          ptElem[1] = currFace+nbf;
          ptElem[2] = currFace+2*nbf;
          ptElem[3] = currFace+1;
        }
        else
        {
          ptElem[1] = currFace+1;
          ptElem[2] = currFace+2*nbf;
          ptElem[3] = currFace+nbf;
        }
        ptElem += 4;

        if (ft[i] == NULL) ft[i] = new vector<E_Int>;
        if1 = abs(ordered[nbf-1]);
        if2 = abs(ordered[0]);
        e1 = cFE(if1-1,1);
        e2 = cFE(if1-1,2);
        e3 = cFE(if2-1,1);
        e4 = cFE(if2-1,2);
        e = 0;
        if (e1 == e3 && e1 != 0) e = e1;
        else if (e1 == e4 && e1 != 0) e = e1;
        else if (e2 == e3 && e2 != 0) e = e2;
        else if (e2 == e4 && e2 != 0) e = e2;
        ft[i]->push_back(currFace+2*nbf); // stockage de la face opposee pour le noeud i 
        ft[i]->push_back(e); 
      }

      if (boundary == false) currFace += 2*nbf;
      else currFace += 2*nbf-1;
    }
  }

  // Re-Creation des elements disparus (a la fin de la connectivite)
  for (E_Int i = 0; i < nelts; i++)
  {
    if (em[i] > 0)
    {
      pt = cNG+posElts[i]; // ancien element

      ptElem[0] = pt[0]+em[i]; // des faces en plus

      //for (E_Int j = 0; j < pt[0]; j++) ptElem[j+1] = pt[j+1]; // anciennes faces
      //for (E_Int j = 0; j < em[i]; j++) // nouvelles faces
      //  ptElem[j+pt[0]+1] = 1; // fake
      nf = pt[0]; jloc = 0; // nbre de faces initiales
      for (E_Int j = 0; j < nf; j++)
      {
        ptElem[jloc+1] = pt[j+1]; // recopie de la face precedente

        // Recherche si il faut inserer une face triangulaire
        if1 = pt[j+1];
        if (j == nf-1) if2 = pt[1];
        else if2 = pt[j+2];
        ptf1 = cNG+posFaces[if1-1];
        ptf2 = cNG+posFaces[if2-1];
        node = 0;
        if (ptf1[1] == ptf2[1]) node = ptf1[1];
        else if (ptf1[1] == ptf2[2]) node = ptf1[1];
        else if (ptf1[2] == ptf2[1]) node = ptf1[2];
        else if (ptf1[2] == ptf2[2]) node = ptf1[2];
        vector<E_Int>* tft = ft[node-1];
        if (tft != NULL)
        {
          for (std::size_t k = 0; k < tft->size()/2; k++)
          {
            indft = (*tft)[2*k]; // no de la face triangulaire opposee
            elt = (*tft)[2*k+1]; // element ou est le triangle
            if (elt == i+1)
            {
              //printf("for element " SF_D_ " : found tft " SF_D2_ "\n", i+1, indft, elt);
              //printf("found tri face " SF_D_ "\n", indft);
              ptElem[jloc+2] = indft; jloc++;    
              break;
            }
          } 
        }
        jloc++;
      } 
      ptElem += pt[0]+em[i]+1;
    }
  }

  for (E_Int i = 0; i < npts; i++) 
  {
    if (ft[i] != NULL) delete ft[i];
  } 
  delete [] ft;
  RELEASESHAREDB(res1, array, f1, cn1);
  RELEASESHAREDB(res2, normal, f2, cn2);
  return tpl;
}
