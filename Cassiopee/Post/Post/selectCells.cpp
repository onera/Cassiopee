/*    
    Copyright 2013-2024 Onera.

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

// selectCells
# include "stdio.h"
# include "post.h"
# include "Nuga/include/ngon_t.hxx"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Selectionne les cellules taggees d'un array */
// ============================================================================
PyObject* K_POST::selectCellsBoth(PyObject* self, PyObject* args)
{
  PyObject* arrayNodes;  PyObject* arrayCenters; PyObject* tag;
  PyObject* PE;
  E_Int strict; E_Int cleanConnectivity;
  
  if (!PYPARSETUPLE_(args, OOO_ I_ O_ I_,
		                 &arrayNodes, &arrayCenters, &tag, &strict, &PE,
                     &cleanConnectivity))
  {
      return NULL;
  }
  // Extract array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cnp;
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray3(arrayNodes, varString, f,
                                     ni, nj, nk, cnp, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: arrayNodes is invalid.");
    return NULL;
  }

  // Extract arrayCenters 
  char* varStringC; char* eltTypeC;
  FldArrayF* fC; FldArrayI* cnpC;
  E_Int resC, niC, njC, nkC;
  resC = K_ARRAY::getFromArray3(arrayCenters, varStringC, fC,
  				                      niC, njC, nkC, cnpC, eltTypeC);
  
  if (resC != 1 && resC != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: arrayCenters is invalid.");
  	return NULL;
  }

  // Extract tag
  char* varString2; char* eltType2;
  FldArrayF* f2; FldArrayI* cnp2;
  E_Int ni2, nj2, nk2;
  E_Int res2 = K_ARRAY::getFromArray3(tag, varString2, f2,
                                      ni2, nj2, nk2, cnp2, eltType2);

  if (res2 != 1 && res2 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag array is invalid.");
    RELEASESHAREDB(res, arrayNodes, f, cnp); 
    return NULL;
  }
  if (res != res2)
  {  
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag and array must represent the same grid.");
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    RELEASESHAREDB(res2, tag, f2, cnp2);
    return NULL;
  }
  if (f->getSize() != f2->getSize())
  {
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    RELEASESHAREDB(res2, tag, f2, cnp2);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag and array must represent the same grid.");
    return NULL;
  }

  if (res == 2 && (cnp->getSize() != cnp2->getSize() ||
                   cnp->getNfld() != cnp2->getNfld())) 
  {
    RELEASESHAREDU(arrayNodes, f, cnp); 
    RELEASESHAREDB(res2, tag, f2, cnp2);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag and array must represent the same grid.");
    return NULL;  
  }

  if (f2->getNfld() != 1) 
  {
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    RELEASESHAREDB(res2, tag, f2, cnp2);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag must have one variable only.");
    return NULL;
  }

  E_Float oneEps = 1.-1.e-10;
  E_Int elt = -1;
  E_Int vertex;
  // no check of coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (res == 1) // Create connectivity cnp (HEXA or QUAD)
  {
    E_Int dim0 = 3;
    if (ni == 1 || nj == 1 || nk == 1) dim0 = 2;
    if (nj == 1 && nk == 1) dim0 = 1;
    else if (ni == 1 && nk == 1) dim0 = 1;
    else if (ni == 1 && nj == 1) dim0 = 1;
    eltType = new char [128];
    if (dim0 == 3) strcpy(eltType, "HEXA");
    else if (dim0 == 2) strcpy(eltType, "QUAD");
    else strcpy(eltType, "BAR");
    E_Int ni1 = E_max(1, E_Int(ni)-1);
    E_Int nj1 = E_max(1, E_Int(nj)-1);
    E_Int nk1 = E_max(1, E_Int(nk)-1);
    E_Int ninj = ni*nj;
    E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
    cnp = new FldArrayI();
    FldArrayI& cn = *cnp;
    E_Int nelts; // nb d'elements non structures
    
    switch (dim0)
    {
      case 1:
      {
        nelts = ncells;
        cn.malloc(nelts, 2);
        E_Int* cn1 = cn.begin(1);
        E_Int* cn2 = cn.begin(2);
        E_Int ind1, ind2;
        elt = 1; //BAR
        if (nk1 == 1 && nj1 == 1)
        {
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = i + 1;
            ind2 = ind1 + 1;
            cn1[i] = ind1; cn2[i] = ind2;
          }
        }
        else if (ni1 == 1 && nj1 == 1)
        {
          for (E_Int k = 0; k < nk1; k++)
          {
            ind1 = k*ni*nj + 1;
            ind2 = ind1 + ni*nj;
            cn1[k] = ind1; cn2[k] = ind2;
          }
        }
        else if (ni1 == 1 && nk1 == 1)
        {
          for (E_Int j = 0; j < nj1; j++)
          {
            ind1 = j*ni + 1;
            ind2 = ind1 + ni;
            cn1[j] = ind1; cn2[j] = ind2;
          }
        }
      }
      break;
    
      case 2:
      {
        nelts = ncells;
        cn.malloc(nelts, 4);
        elt = 3; // QUAD
        E_Int* cn1 = cn.begin(1);
        E_Int* cn2 = cn.begin(2);
        E_Int* cn3 = cn.begin(3);
        E_Int* cn4 = cn.begin(4);

        if (nk1 == 1)
        {
#pragma omp parallel default(shared) if (nj1 > __MIN_SIZE_MEAN__)
          {
            E_Int ind, ind1, ind2, ind3, ind4;
#pragma omp for
            for (E_Int j = 0; j < nj1; j++)
              for (E_Int i = 0; i < ni1; i++)
              {
                //starts from 1
                ind1 = i + j*ni + 1; //(i,j,1)
                ind2 = ind1 + 1;  //(i+1,j,1)
                ind3 = ind2 + ni; //(i+1,j+1,1)
                ind4 = ind3 - 1;  //(i,j+1,1)
                ind = i + j*ni1;
                cn1[ind] = ind1; cn2[ind] = ind2;
                cn3[ind] = ind3; cn4[ind] = ind4;
              }
          }
        }
        else if (nj1 == 1)
        {
#pragma omp parallel default(shared) if (nk1 > __MIN_SIZE_MEAN__)
          {
            E_Int ind, ind1, ind2, ind3, ind4;
#pragma omp for
            for (E_Int k = 0; k < nk1; k++)
              for (E_Int i = 0; i < ni1; i++)
              {
                ind1 = i + k*ninj + 1;  //(i,1,k)
                ind2 = ind1 + ninj; //(i,1,k+1)
                ind3 = ind2 + 1;    //(i+1,1,k+1)
                ind4 = ind3 - 1;    //(i,1,k+1)
                ind = i + k*ni1;
                cn1[ind] = ind1; cn2[ind] = ind2;
                cn3[ind] = ind3; cn4[ind] = ind4;      
              }
          }
        }
        else // i1 = 1 
        {
#pragma omp parallel default(shared) if (nk1 > __MIN_SIZE_MEAN__)
          {
            E_Int ind, ind1, ind2, ind3, ind4;
#pragma omp for   
            for (E_Int k = 0; k < nk1; k++)
              for (E_Int j = 0; j < nj1; j++)
              {
                ind1 = 1 + j*ni + k*ninj; //(1,j,k)
                ind2 = ind1 + ni;   //(1,j+1,k)
                ind3 = ind2 + ninj; //(1,j+1,k+1)
                ind4 = ind3 - ni;   //(1,j,k+1)
                ind = j+k*nj1;
                cn1[ind] = ind1; cn2[ind] = ind2;
                cn3[ind] = ind3; cn4[ind] = ind4;
              }
          }
        }// i1 = 1
      }
      break;

      case 3:
      { 
        nelts = ncells;
        cn.malloc(nelts,8);
        elt = 7; //HEXA
        E_Int* cn1 = cn.begin(1);
        E_Int* cn2 = cn.begin(2);
        E_Int* cn3 = cn.begin(3);
        E_Int* cn4 = cn.begin(4);
        E_Int* cn5 = cn.begin(5);
        E_Int* cn6 = cn.begin(6);
        E_Int* cn7 = cn.begin(7);
        E_Int* cn8 = cn.begin(8);
#pragma omp parallel default(shared) if (nk1 > __MIN_SIZE_MEAN__)
        {
          E_Int ind, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
#pragma omp for
            for (E_Int k = 0; k < nk1; k++)
              for (E_Int j = 0; j < nj1; j++)
                for (E_Int i = 0; i < ni1; i++)
                {
                  ind1 = 1 + i + j*ni + k*ninj; //A(  i,  j,k)
                  ind2 = ind1 + 1;              //B(i+1,  j,k)
                  ind3 = ind2 + ni;            //C(i+1,j+1,k)
                  ind4 = ind3 - 1;              //D(  i,j+1,k)
                  ind5 = ind1 + ninj;           //E(  i,  j,k+1)
                  ind6 = ind2 + ninj;           //F(i+1,  j,k+1)
                  ind7 = ind3 + ninj;           //G(i+1,j+1,k+1)
                  ind8 = ind4 + ninj;           //H(  i,j+1,k+1) 
                  ind = i+j*ni1+k*ni1*nj1;
                  cn1[ind] = ind1; cn2[ind] = ind2;
                  cn3[ind] = ind3; cn4[ind] = ind4;
                  cn5[ind] = ind5; cn6[ind] = ind6;
                  cn7[ind] = ind7; cn8[ind] = ind8;
                }
          }
      }
      break;
    }
  }

  PyObject* l = PyList_New(0); 
  
  // Infos sur le type d'element
  E_Int isNGon = 1; E_Int isNode = 1;
  isNGon = strcmp(eltType, "NGON"); // vaut 0 si l'elmt est un NGON
  isNode = strcmp(eltType, "NODE"); // vaut 0 si l'elmt est un NODE
    
  // Selection
  PyObject* tpl;
  PyObject* tplc;
  if (isNGon != 0 && isNode != 0) // tous les elements sauf NGON et NODE
  {
    E_Int nfld = f->getNfld();
    E_Int nt = cnp->getNfld();
    E_Int csize = cnp->getSize();
    FldArrayF* an = new FldArrayF();
    FldArrayF& coord = *an;
    FldArrayI* acn = new FldArrayI();
    FldArrayI& cn = *acn;

    // Selection des vertex
    FldArrayI selected(f->getSize(), 2);
    selected.setAllValuesAtNull();
    E_Int* selected1 = selected.begin(1);
    E_Int* selected2 = selected.begin(2);
    E_Float* tagp = f2->begin();
    E_Int isSel = 0;

    // Selection des champs aux centres
    E_Int nfldC         = fC->getNfld(); // nombre de champs en centres 
    FldArrayF* foutC    = new FldArrayF(*fC);
    FldArrayF& fcenter  = *foutC; 
    FldArrayF& fcenter0 = *fC; 
    E_Int ii            = 0;
   
    if (strict == 0) // cell selectionnee des qu'un sommet est tag=1
    {
      for (E_Int i = 0; i < csize; i++)
      {
        for (E_Int nv = 1; nv <= nt; nv++)
        {
          vertex = (*cnp)(i, nv);
          if (tagp[vertex-1] >= oneEps) 
          {  
            for (E_Int nv = 1; nv <= nt; nv++)
            {
              vertex = (*cnp)(i, nv);
              selected1[vertex-1] = 1;
            }
	          // champs en centres
            for (E_Int k = 1; k <= nfldC; k++)
            {
              fcenter(ii,k) = fcenter0(i,k); 
            }
            ii++;
            break; 
          }
        }
      }
    }
    else // cell selectionnee si tous les sommets sont tag=1
    {
      for (E_Int i = 0; i < csize; i++)
      {
        isSel = 0;
        for (E_Int nv = 1; nv <= nt; nv++)
        {
          vertex = (*cnp)(i, nv);
          if (tagp[vertex-1] >= oneEps) isSel++;
        }
        if (isSel == nt)
        {
          for (E_Int nv = 1; nv <= nt; nv++)
          {
            vertex = (*cnp)(i, nv);
            selected1[vertex-1] = 1;
          }
	       // champs en centres
	       for (E_Int k = 1; k <= nfldC; k++)
	       {
	         fcenter(ii,k) = fcenter0(i,k); 
	       }
	       ii++;
        }
      }    
    }

    // Tableau des champs en centre
    foutC->reAllocMat(ii,nfldC);
    
    E_Int cprev = 0;
    selected2[0] = cprev;
    for (E_Int i = 1; i < selected.getSize(); i++)
    {
      if (selected1[i-1] == 1) selected2[i] = cprev+1; // real position
      else selected2[i] = cprev;
      cprev = selected2[i];
    }

    // Tableau des vertex
    if (selected1[selected.getSize()-1] == 1) cprev++;
    coord.malloc(cprev, nfld);

    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* coordp = coord.begin(n);
      E_Float* fp = f->begin(n);
      cprev = 0;
      for (E_Int i = 0; i < f->getSize(); i++)
      {
        if (selected1[i] == 1)
        {
          coordp[cprev] = fp[i]; cprev++;
        }
      }
    }
  
    // Tableau des connectivites
    cn.malloc(csize, nt);
  
    cprev = 0;
    E_Int p;
    for (E_Int i = 0; i < csize; i++)
    {
      p = 1;
      for (E_Int n = 1; n <= nt; n++)
      {
        p = p*selected1[(*cnp)(i,n)-1];
      }
      if (p == 1)
      {
        for (E_Int n = 1; n <= nt; n++)
        {
          vertex = (*cnp)(i, n);
          cn(cprev,n) = selected2[vertex-1]+1;
        }
        cprev++;
      }
    }

    RELEASESHAREDB(res2, tag, f2, cnp2);
    if (res == 1) delete cnp;
   
    cn.reAllocMat(cprev, nt);
    if (cprev == 0) an->reAllocMat(0,nfld);
    else 
    {
      if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
        K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, *an, *acn);
    }
    tpl  = K_ARRAY::buildArray(*an,    varString,  *acn, elt, eltType);
    tplc = K_ARRAY::buildArray(*foutC, varStringC, *acn, elt, eltTypeC);
    
    delete foutC ; 
    delete an; delete acn;
    if (res == 1) delete[] eltType;
  }
  else if (isNode == 0) // NODE
  {
    E_Int nfld = f->getNfld();
    FldArrayF* an = new FldArrayF();
    FldArrayF& coord = *an;

    // Selection des vertex
    FldArrayI selected(f->getSize(), 2);
    selected.setAllValuesAtNull();
    E_Int* selectedp = selected.begin();
    E_Float* tagp = f2->begin();
    E_Int count = 0;
    for (E_Int i = 0; i < f->getSize(); i++)
    {
      if (tagp[i] >= oneEps) { selectedp[i] = 1; count++; }
    }

    // Tableau des vertex
    coord.malloc(count, nfld);

    for (E_Int n = 1; n <= nfld; n++)
    {
      count = 0;
      E_Float* coordn = coord.begin(n);
      E_Float* fn = f->begin(n);
      for (E_Int i = 0; i < f->getSize(); i++)
      {
        if (selectedp[i] == 1)
        {
          coordn[count] = fn[i];
          count++;
        }
      }
    }
  
    // Tableau des connectivites
    FldArrayI* acn = new FldArrayI();
    FldArrayI& cn = *acn; cn.malloc(0, 1);
  
    RELEASESHAREDB(res2, tag, f2, cnp2);
    if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, *an, *acn);
    tpl = K_ARRAY::buildArray(*an, varString, *acn, elt, eltType);
    tplc = tpl ;
    delete an; delete acn;
  }
  else // elements NGON
  {
    FldArrayF* fout = new FldArrayF(*f);
    // Algo: 
    // Parcours des faces pour selectionner les valides (selectionnees selon le critere "strict")
    // Parcours les elmts. 
    // Parcours des faces de l'elmt. 
    // Si strict=0 et 1 face valide ou si strict=1 et toutes les faces valides, 
    // On ajoute l'elmt dans une nouvelle connectivite Elmt/Face
    // Au final, on reconstruit une connectivite complete en ajoutant l'ancienne connectivite Face/Noeuds 
    // Et on apelle cleanConnectivity
    // ----------------------------
    E_Float* tagp = f2->begin();   // champ indiquant le critere de selection
    E_Int* cnpp = cnp->begin();    // pointeur sur l ancienne connectivite
    E_Int nbFaces = cnpp[0];       // nombre de faces de l'ancienne connectivite
    E_Int sizeFN = cnpp[1];        // taille de l'ancienne connectivite Face/Noeuds
    E_Int nbElements = cnpp[sizeFN+2]; // nombre d'elts de l'ancienne connectivite
    E_Int sizeEF = cnpp[sizeFN+3];     // taille de l'ancienne connectivite Elmt/Faces
    E_Int* cnEFp = cnpp + 4+sizeFN;  // pointeur sur l'ancienne connectivite Elmt/Faces
    FldArrayI cn2 = FldArrayI(sizeEF); // nouvelle connectivite Elmt/Faces
    E_Int* cn2p = cn2.begin();       // pointeur sur la nouvelle connectivite
    E_Int size2 = 0;                 // compteur pour la nouvelle connectivite
    E_Int next = 0;                  // nbre d'elts selectionnes
    FldArrayI selectedFaces(nbFaces);  // tableau indiquant les faces valides 
    E_Int* selectedFacesp = selectedFaces.begin();
    E_Int nbfaces, nbnodes;
    E_Int isSel;

    // Champs en centre  
    FldArrayF* foutC    = new FldArrayF(*fC);
    FldArrayF& fcenter  = *foutC;
    FldArrayF& fcenter0 = *fC;
    E_Int nfldC         = fC->getNfld();
    E_Int ii            = 0;
    
    // Selection des faces valides
    E_Int fa = 0; E_Int numFace = 0; cnpp += 2;

    if (strict == 0)  // cell selectionnee des qu'un sommet est tag=1
    {
      while (fa < sizeFN) // parcours de la connectivite face/noeuds
      {
        nbnodes = cnpp[0];
        selectedFacesp[numFace] = 0;
        for (E_Int n = 1; n <= nbnodes; n++)
        {
          if (tagp[cnpp[n]-1] >= oneEps) { selectedFacesp[numFace] = 1; break; }   
        }
        cnpp += nbnodes+1; fa += nbnodes+1; numFace++;
      }
    }
    else //strict=1, cell selectionnee si tous les sommets sont tag
    {
      while (fa < sizeFN) // parcours de la connectivite face/noeuds
      {
        nbnodes = cnpp[0];
        selectedFacesp[numFace] = 1;
        for (E_Int n = 1; n <= nbnodes; n++)
        {
          if (tagp[cnpp[n]-1] < oneEps) { selectedFacesp[numFace] = 0; break; }   
        }
        cnpp += nbnodes+1; fa += nbnodes+1; numFace++;
      }
    }
	   
    
    // Si mise a jour du ParentElement, tab d'indirection des faces et des elmts
    // ------------------------------------------------------------------------
    E_Int newNumFace = 0;
    FldArrayI new_pg_ids;
    FldArrayI keep_pg;
    FldArrayI new_ph_ids;
    
    if (PE != Py_None)
    {   
      new_pg_ids.malloc(nbFaces);    // Tableau d'indirection des faces (pour maj PE)
      keep_pg.malloc(nbFaces);       // Flag de conservation des faces 
      new_ph_ids.malloc(nbElements); // Tableau d'indirection des elmts (pour maj PE)

      new_pg_ids = -1;
      new_ph_ids = -1;
      keep_pg    = -1; 
      
      // Selection des elements en fonction des faces valides
      if (strict == 0)  // cell selectionnee des qu'un sommet est tag=1
      {
        for (E_Int i = 0; i < nbElements; i++)
        {
	        new_ph_ids[i] = -1;
          nbfaces       = cnEFp[0];
          for (E_Int n = 1; n <= nbfaces; n++)
          {
            if (selectedFacesp[cnEFp[n]-1] == 1) 
            { 
              cn2p[0] = nbfaces; size2 += 1;
              for (E_Int n = 1; n <= nbfaces; n++)
              {
                cn2p[n] = cnEFp[n];
                keep_pg[cnEFp[n]-1] = +1;
              }
              size2 += nbfaces; cn2p += nbfaces+1; next++;
        
              for (E_Int k = 1; k <= nfldC; k++)
                fcenter(ii,k) = fcenter0(i,k);
    
              new_ph_ids[i] = ii;
              ii++;  
              break; 
            }
          }
          cnEFp += nbfaces+1; 
        }	  
      }
      else //strict=1, cell selectionnee si tous les sommets sont tag
      {
        for (E_Int i = 0; i < nbElements; i++)
        {
          isSel         = 0;
          new_ph_ids[i] = -1;
          nbfaces       = cnEFp[0];
          for (E_Int n = 1; n <= nbfaces; n++)
          {
            if (selectedFacesp[cnEFp[n]-1] == 1) isSel++;
          }
          if (isSel == nbfaces) //cell selectionnee si tous les sommets sont tag
          {
            cn2p[0] = nbfaces; size2 +=1;
            for (E_Int n = 1; n <= nbfaces; n++)
            {
              cn2p[n] = cnEFp[n];
              keep_pg[cnEFp[n]-1] = +1;	
            }
            size2 += nbfaces; cn2p += nbfaces+1; next++;
      
            for (E_Int k = 1; k <= nfldC; k++) fcenter(ii,k) = fcenter0(i,k);
  
            new_ph_ids[i] = ii;
  
            ii++;
          }
          cnEFp += nbfaces+1; 
        } 
      }
      cn2.reAlloc(size2);
      foutC->reAllocMat(ii,nfldC);

      E_Int nn = 0; 
      for (E_Int n = 0; n<new_pg_ids.getSize(); n++)
      {
        if (keep_pg[n]>0){ new_pg_ids[n] = nn; nn++; newNumFace++;}
      }	
    }
    else  // PE == Py_None - pas de creation de tab d'indirection 
    {
      // Selection des elements en fonction des faces valides
      if (strict == 0)  // cell selectionnee des qu'un sommet est tag=1
      {
        for (E_Int i = 0; i < nbElements; i++)
        {
          nbfaces = cnEFp[0];
          for (E_Int n = 1; n <= nbfaces; n++)
          {
            if (selectedFacesp[cnEFp[n]-1] == 1) 
            { 
              cn2p[0] = nbfaces; size2 += 1;
              for (E_Int n = 1; n <= nbfaces; n++) cn2p[n] = cnEFp[n];
              size2 += nbfaces; cn2p += nbfaces+1; next++;
        
              for (E_Int k = 1; k <= nfldC; k++)
                fcenter(ii,k) = fcenter0(i,k);
              ii++;
                
              break; 
            }
          }
          cnEFp += nbfaces+1; 
        }
      }
      else //strict=1, cell selectionnee si tous les sommets sont tag
      {
        for (E_Int i = 0; i < nbElements; i++)
        {
          isSel = 0;
          nbfaces = cnEFp[0];
          for (E_Int n = 1; n <= nbfaces; n++)
          {
            if (selectedFacesp[cnEFp[n]-1] == 1) isSel++;
          }
          if (isSel == nbfaces) //cell selectionnee si tous les sommets sont tag
          {
            cn2p[0] = nbfaces; size2 +=1;
            for (E_Int n = 1; n <= nbfaces; n++) cn2p[n] = cnEFp[n];
            size2 += nbfaces; cn2p += nbfaces+1; next++;
      
            for (E_Int k = 1; k <= nfldC; k++) fcenter(ii,k) = fcenter0(i,k);
            ii++;
          }
          cnEFp += nbfaces+1; 
        } 
      }
      cn2.reAlloc(size2);
      foutC->reAllocMat(ii,nfldC);
    }

    // Cree la nouvelle connectivite complete
    E_Int coutsize = sizeFN+4+size2;
    FldArrayI* cout = new FldArrayI(coutsize);
    E_Int* coutp = cout->begin();
    cnpp = cnp->begin(); cn2p = cn2.begin(); 
    for (E_Int i = 0; i < sizeFN+2; i++) coutp[i] = cnpp[i];
    coutp += sizeFN+2;
    coutp[0] = next;
    coutp[1] = size2; coutp += 2;
    for (E_Int i = 0; i < size2; i++) coutp[i] = cn2p[i];

    RELEASESHAREDB(res2, tag, f2, cnp2);


    if (PE != Py_None)
    {
      // Check numpy (parentElement)
      FldArrayI* cFE;
      E_Int res = K_NUMPY::getFromNumpyArray(PE, cFE, true);
      
      if (res == 0)
      {
	      RELEASESHAREDN(PE, cFE);
	      PyErr_SetString(PyExc_TypeError, "selectCellsBoth: PE numpy is invalid.");
	      return NULL;
      }
      
      ngon_t<K_FLD::FldArrayI> ng(*cout); // construction d'un ngon_t à partir d'un FldArrayI

      FldArrayI* cFEp_new = new FldArrayI(newNumFace,2);
      FldArrayI& cFE_new  = *cFEp_new ; 

      E_Int* cFEl = cFE_new.begin(1);
      E_Int* cFEr = cFE_new.begin(2);

      E_Int* cFEl_old = cFE->begin(1);
      E_Int* cFEr_old = cFE->begin(2);
      
      E_Int old_ph_1, old_ph_2;

      for (E_Int pgi = 0; pgi < nbFaces; pgi++)
      {
	     if (new_pg_ids[pgi]>=0)
	     {
	       old_ph_1 = cFEl_old[pgi]-1;
	       old_ph_2 = cFEr_old[pgi]-1;

	       if (old_ph_1 >= 0) // l'elmt gauche existe 
	       {
	         cFEl[new_pg_ids[pgi]] = new_ph_ids[old_ph_1]+1;
 
	         if (old_ph_2 >= 0) // l'elmt droit existe 
	           cFEr[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1; 
	         else
	           cFEr[new_pg_ids[pgi]] = 0;
	       } 
	       else // l'elmt gauche a disparu - switch droite/gauche
	       {
	         cFEl[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1;
	         cFEr[new_pg_ids[pgi]] = 0;
	         // reverse
	         E_Int s = ng.PGs.stride(pgi);
	         E_Int* p = ng.PGs.get_facets_ptr(pgi);
	         std::reverse(p, p + s);
	       }
	     }
      } // boucle pgi      
      
      // export ngon
      ng.export_to_array(*cout);

      // objet Python de sortie
      PyObject* pyPE = K_NUMPY::buildNumpyArray(cFE_new, 1);

      PyList_Append(l,pyPE);
      
      delete cFEp_new;
      RELEASESHAREDN(PE, cFE);
    }

   
    // close
    if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
      K_CONNECT::cleanConnectivityNGon(posx, posy, posz, 1.e-10, *fout, *cout);
    
    tpl  = K_ARRAY::buildArray(*fout,   varString, *cout, 8);
    tplc = K_ARRAY::buildArray(*foutC, varStringC, *cout, 8);
    
    delete fout; delete foutC; delete cout;
  }

  RELEASESHAREDB(res, arrayNodes, f, cnp);
  RELEASESHAREDB(resC, arrayCenters, fC, cnpC);
  
  PyList_Append(l,tpl) ; Py_DECREF(tpl);
  PyList_Append(l,tplc); Py_DECREF(tplc);

  return l;  
}

//==============================================================================
// Retourne les indices des elements correspondants a tag=1
// IN: tout type 
// IN: tag au centres
// OUT: liste d''elements
//==============================================================================
PyObject* K_POST::selectCells3(PyObject* self, PyObject* args)
{
  PyObject* tag; E_Int flag;
  if (!PYPARSETUPLE_(args, O_ I_, &tag, &flag))
  {
      return NULL;
  }

  // Extract tag (centers)
  char* varString2; char* eltType2;
  FldArrayF* f2; FldArrayI* cn2;
  E_Int ni2, nj2, nk2;
  E_Int res2;

  if (flag == 0) // array
  {
    res2 = K_ARRAY::getFromArray3(tag, varString2, f2,
                                  ni2, nj2, nk2, cn2, eltType2);
    if (res2 != 1 && res2 != 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "selectCells3: tag array is invalid.");
      return NULL;
    }
  }
  else // numpy
  {
    res2 = K_NUMPY::getFromNumpyArray(tag, f2, true);
    if (res2 == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "selectCells3: tag numpy is invalid.");
      return NULL;
    }
  }

  E_Int nelts = f2->getSize();
  printf("nelts=" SF_D_ "\n", nelts);
  E_Float* tagp = f2->begin(); 
  E_Float oneEps = 1.-1.e-10;

  // Compte les tag > 1-eps
  E_Int nthreads = __NUMTHREADS__;
  E_Int* nb = new E_Int [nthreads];
  E_Int* pos = new E_Int[nthreads];
  for (E_Int i = 0; i < nthreads; i++) nb[i] = 0;

#pragma omp parallel default(shared)
  {
    E_Int t = __CURRENT_THREAD__;

#pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      if (tagp[i] > oneEps) nb[t]++;
    }
  }

  // allocate
  E_Int size = 0;
  for (E_Int i = 0; i < nthreads; i++) { pos[i] = size; size += nb[i]; }
  printf("size = " SF_D_ "\n", size);

  PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1, 0);
  E_Int* ind = K_NUMPY::getNumpyPtrI(o);

#pragma omp parallel default(shared)
  {
    E_Int t = __CURRENT_THREAD__;
    E_Int c = 0;
    E_Int* ptr = ind+pos[t];

#pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      if (tagp[i] > oneEps) { ptr[c] = i; c++; }
    }
  }

  delete [] nb; delete [] pos;
  
  if (flag == 0) { RELEASESHAREDB(res2, tag, f2, cn2); }
  else { RELEASESHAREDN(tag, f2); }
  return o;
}

//=============================================================================
/* Selectionne les cellules taggees d'un array */
// ============================================================================
PyObject* K_POST::selectCells(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* tag;
  PyObject* PE;
  E_Int strict; E_Int cleanConnectivity;
  
  if (!PYPARSETUPLE_(args, OO_ I_ O_ I_,
                    &array, &tag, &strict, &PE, &cleanConnectivity))
  {
      return NULL;
  }
  
  // Extract array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cnp;
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, 
                                    cnp, eltType, true);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: array is invalid.");
    return NULL;
  }

  // Extract tag
  char* varString2; char* eltType2;
  FldArrayF* f2; FldArrayI* cnp2;
  E_Int ni2, nj2, nk2;
  E_Int res2 = 
    K_ARRAY::getFromArray(tag, varString2, f2, ni2, nj2, nk2, cnp2, eltType2, 
                          true);

  if (res2 != 1 && res2 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag array is invalid.");
    RELEASESHAREDB(res, array, f, cnp); 
    return NULL;
  }
  if (res != res2)
  {  
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag and array must represent the same grid.");
    RELEASESHAREDB(res, array, f, cnp);
    RELEASESHAREDB(res2, tag, f2, cnp2);
    return NULL;
  }
  if (f->getSize() != f2->getSize())
  {
    RELEASESHAREDB(res, array, f, cnp);
    RELEASESHAREDB(res2, tag, f2, cnp2);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag and array must represent the same grid.");
    return NULL;
  }

  if (res == 2 && (cnp->getSize() != cnp2->getSize() ||
                   cnp->getNfld() != cnp2->getNfld())) 
  {
    RELEASESHAREDU(array, f, cnp); 
    RELEASESHAREDB(res2, tag, f2, cnp2);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag and array must represent the same grid.");
    return NULL;  
  }

  if (f2->getNfld() != 1) 
  {
    RELEASESHAREDB(res, array, f, cnp);
    RELEASESHAREDB(res2, tag, f2, cnp2);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag must have one variable only.");
    return NULL;
  }

  E_Float oneEps = 1.-1.e-10;
  E_Int elt = -1;
  E_Int vertex;
  // no check of coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (res == 1) // Create connectivity cnp (HEXA or QUAD)
  {
    E_Int dim0 = 3;
    if (ni == 1 || nj == 1 || nk == 1) dim0 = 2;
    if (nj == 1 && nk == 1) dim0 = 1;
    else if (ni == 1 && nk == 1) dim0 = 1;
    else if (ni == 1 && nj == 1) dim0 = 1;
    eltType = new char [128];
    if (dim0 == 3) strcpy(eltType, "HEXA");
    else if (dim0 == 2) strcpy(eltType, "QUAD");
    else strcpy(eltType, "BAR");
    E_Int ni1 = E_max(1, E_Int(ni)-1);
    E_Int nj1 = E_max(1, E_Int(nj)-1);
    E_Int nk1 = E_max(1, E_Int(nk)-1);
    E_Int ninj = ni*nj;
    E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
    cnp = new FldArrayI();
    FldArrayI& cn = *cnp;
    E_Int nelts; // nb d'elements non structures
    
    switch (dim0)
    {
      case 1:
      {
        nelts = ncells;
        cn.malloc(nelts, 2);
        E_Int* cn1 = cn.begin(1);
        E_Int* cn2 = cn.begin(2);
        E_Int ind1, ind2;
        elt = 1; //BAR
        if (nk1 == 1 && nj1 == 1)
        {
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = i + 1;
            ind2 = ind1 + 1;
            cn1[i] = ind1; cn2[i] = ind2;
          }
        }
        else if (ni1 == 1 && nj1 == 1)
        {
          for (E_Int k = 0; k < nk1; k++)
          {
            ind1 = k*ni*nj + 1;
            ind2 = ind1 + ni*nj;
            cn1[k] = ind1; cn2[k] = ind2;
          }
        }
        else if (ni1 == 1 && nk1 == 1)
        {
          for (E_Int j = 0; j < nj1; j++)
          {
            ind1 = j*ni + 1;
            ind2 = ind1 + ni;
            cn1[j] = ind1; cn2[j] = ind2;
          }
        }
      }
      break;
    
      case 2:
      {
        nelts = ncells;
        cn.malloc(nelts, 4);
        elt = 3; // QUAD
        E_Int* cn1 = cn.begin(1);
        E_Int* cn2 = cn.begin(2);
        E_Int* cn3 = cn.begin(3);
        E_Int* cn4 = cn.begin(4);

        if (nk1 == 1)
        {
#pragma omp parallel default(shared) if (nj1 > __MIN_SIZE_MEAN__)
          {
            E_Int ind, ind1, ind2, ind3, ind4;
#pragma omp for
            for (E_Int j = 0; j < nj1; j++)
              for (E_Int i = 0; i < ni1; i++)
              {
                //starts from 1
                ind1 = i + j*ni + 1; //(i,j,1)
                ind2 = ind1 + 1;  //(i+1,j,1)
                ind3 = ind2 + ni; //(i+1,j+1,1)
                ind4 = ind3 - 1;  //(i,j+1,1)
                ind = i + j*ni1;
                cn1[ind] = ind1; cn2[ind] = ind2;
                cn3[ind] = ind3; cn4[ind] = ind4;
              }
          }
        }
        else if (nj1 == 1)
        {
#pragma omp parallel default(shared) if (nk1 > __MIN_SIZE_MEAN__)
          {
            E_Int ind, ind1, ind2, ind3, ind4;
#pragma omp for
            for (E_Int k = 0; k < nk1; k++)
              for (E_Int i = 0; i < ni1; i++)
              {
                ind1 = i + k*ninj + 1;  //(i,1,k)
                ind2 = ind1 + ninj; //(i,1,k+1)
                ind3 = ind2 + 1;    //(i+1,1,k+1)
                ind4 = ind3 - 1;    //(i,1,k+1)
                ind = i + k*ni1;
                cn1[ind] = ind1; cn2[ind] = ind2;
                cn3[ind] = ind3; cn4[ind] = ind4;      
              }
          }
        }
        else // i1 = 1 
        {
#pragma omp parallel default(shared) if (nk1 > __MIN_SIZE_MEAN__)
          {
            E_Int ind, ind1, ind2, ind3, ind4;
#pragma omp for   
            for (E_Int k = 0; k < nk1; k++)
              for (E_Int j = 0; j < nj1; j++)
              {
                ind1 = 1 + j*ni + k*ninj; //(1,j,k)
                ind2 = ind1 + ni;   //(1,j+1,k)
                ind3 = ind2 + ninj; //(1,j+1,k+1)
                ind4 = ind3 - ni;   //(1,j,k+1)
                ind = j+k*nj1;
                cn1[ind] = ind1; cn2[ind] = ind2;
                cn3[ind] = ind3; cn4[ind] = ind4;
              }
          }
        }// i1 = 1
      }
      break;

      case 3:
      { 
        nelts = ncells;
        cn.malloc(nelts,8);
        elt = 7; //HEXA
        E_Int* cn1 = cn.begin(1);
        E_Int* cn2 = cn.begin(2);
        E_Int* cn3 = cn.begin(3);
        E_Int* cn4 = cn.begin(4);
        E_Int* cn5 = cn.begin(5);
        E_Int* cn6 = cn.begin(6);
        E_Int* cn7 = cn.begin(7);
        E_Int* cn8 = cn.begin(8);
#pragma omp parallel default(shared) if (nk1 > __MIN_SIZE_MEAN__)
        {
          E_Int ind, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
#pragma omp for
            for (E_Int k = 0; k < nk1; k++)
              for (E_Int j = 0; j < nj1; j++)
                for (E_Int i = 0; i < ni1; i++)
                {
                  ind1 = 1 + i + j*ni + k*ninj; //A(  i,  j,k)
                  ind2 = ind1 + 1;              //B(i+1,  j,k)
                  ind3 = ind2 + ni;            //C(i+1,j+1,k)
                  ind4 = ind3 - 1;              //D(  i,j+1,k)
                  ind5 = ind1 + ninj;           //E(  i,  j,k+1)
                  ind6 = ind2 + ninj;           //F(i+1,  j,k+1)
                  ind7 = ind3 + ninj;           //G(i+1,j+1,k+1)
                  ind8 = ind4 + ninj;           //H(  i,j+1,k+1) 
                  ind = i+j*ni1+k*ni1*nj1;
                  cn1[ind] = ind1; cn2[ind] = ind2;
                  cn3[ind] = ind3; cn4[ind] = ind4;
                  cn5[ind] = ind5; cn6[ind] = ind6;
                  cn7[ind] = ind7; cn8[ind] = ind8;
                }
          }
      }
      break;
    }
  }

  PyObject* l = PyList_New(0);
  
  // Infos sur le type d'element
  E_Int isNGon = 1; E_Int isNode = 1;
  isNGon = strcmp(eltType, "NGON"); // vaut 0 si l'elmt est un NGON
  isNode = strcmp(eltType, "NODE"); // vaut 0 si l'elmt est un NODE
    
  // Selection
  PyObject* tpl;
  if (isNGon != 0 && isNode != 0) // tous les elements sauf NGON et NODE
  {
    E_Int nfld = f->getNfld();
    E_Int nt = cnp->getNfld();
    E_Int csize = cnp->getSize();
    FldArrayF* an = new FldArrayF();
    FldArrayF& coord = *an;
    FldArrayI* acn = new FldArrayI();
    FldArrayI& cn = *acn;

    // Selection des vertex
    FldArrayI selected(f->getSize(), 2);
    selected.setAllValuesAtNull();
    E_Int* selected1 = selected.begin(1);
    E_Int* selected2 = selected.begin(2);
    E_Float* tagp = f2->begin();
    E_Int isSel = 0;
    
    if (strict == 0) // cell selectionnee des qu'un sommet est tag=1
    {
      for (E_Int i = 0; i < csize; i++)
      {
        for (E_Int nv = 1; nv <= nt; nv++)
        {
          vertex = (*cnp)(i, nv);
          if (tagp[vertex-1] >= oneEps) 
          {  
            for (E_Int nv = 1; nv <= nt; nv++)
            {
              vertex = (*cnp)(i, nv);
              selected1[vertex-1] = 1;
            }
            break; 
          }
        }
      }
    }
    else // cell selectionnee si tous les sommets sont tag=1
    {
      for (E_Int i = 0; i < csize; i++)
      {
        isSel = 0;
        for (E_Int nv = 1; nv <= nt; nv++)
        {
          vertex = (*cnp)(i, nv);
          if (tagp[vertex-1] >= oneEps) isSel++;
        }
        if (isSel == nt)
        {
          for (E_Int nv = 1; nv <= nt; nv++)
          {
            vertex = (*cnp)(i, nv);
            selected1[vertex-1] = 1;
          }
        }
      }    
    }
    E_Int cprev = 0;
    selected2[0] = cprev;
    for (E_Int i = 1; i < selected.getSize(); i++)
    {
      if (selected1[i-1] == 1) selected2[i] = cprev+1; // real position
      else selected2[i] = cprev;
      cprev = selected2[i];
    }

    // Tableau des vertex
    if (selected1[selected.getSize()-1] == 1) cprev++;
    coord.malloc(cprev, nfld);

    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* coordp = coord.begin(n);
      E_Float* fp = f->begin(n);
      cprev = 0;
      for (E_Int i = 0; i < f->getSize(); i++)
      {
        if (selected1[i] == 1)
        {
          coordp[cprev] = fp[i]; cprev++;
        }
      }
    }
  
    // Tableau des connectivites
    cn.malloc(csize, nt);
  
    cprev = 0;
    E_Int p;
    for (E_Int i = 0; i < csize; i++)
    {
      p = 1;
      for (E_Int n = 1; n <= nt; n++)
      {
        p = p*selected1[(*cnp)(i,n)-1];
      }
      if (p == 1)
      {
        for (E_Int n = 1; n <= nt; n++)
        {
          vertex = (*cnp)(i, n);
          cn(cprev,n) = selected2[vertex-1]+1;
        }
        cprev++;
      }
    }

    RELEASESHAREDB(res2, tag, f2, cnp2);
    if (res == 1) delete cnp;
   
    cn.reAllocMat(cprev, nt);
    if (cprev == 0) an->reAllocMat(0,nfld);
    else 
    {
      if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
        K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, *an, *acn);
    }
    tpl = K_ARRAY::buildArray(*an, varString, *acn, elt, eltType);
    delete an; delete acn;
    if (res == 1) delete[] eltType;
  }
  else if (isNode == 0) // NODE
  {
    E_Int nfld = f->getNfld();
    FldArrayF* an = new FldArrayF();
    FldArrayF& coord = *an;

    // Selection des vertex
    FldArrayI selected(f->getSize(), 2);
    selected.setAllValuesAtNull();
    E_Int* selectedp = selected.begin();
    E_Float* tagp = f2->begin();
    E_Int count = 0;
    for (E_Int i = 0; i < f->getSize(); i++)
    {
      if (tagp[i] >= oneEps) { selectedp[i] = 1; count++; }
    }

    // Tableau des vertex
    coord.malloc(count, nfld);

    for (E_Int n = 1; n <= nfld; n++)
    {
      count = 0;
      E_Float* coordn = coord.begin(n);
      E_Float* fn = f->begin(n);
      for (E_Int i = 0; i < f->getSize(); i++)
      {
        if (selectedp[i] == 1)
        {
          coordn[count] = fn[i];
          count++;
        }
      }
    }
  
    // Tableau des connectivites
    FldArrayI* acn = new FldArrayI();
    FldArrayI& cn = *acn; cn.malloc(0, 1);
  
    RELEASESHAREDB(res2, tag, f2, cnp2);
    if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, *an, *acn);
    tpl = K_ARRAY::buildArray(*an, varString, *acn, elt, eltType);
    delete an; delete acn;
  }
  else // elements NGON
  {
    FldArrayF* fout = new FldArrayF(*f);
    // Algo: 
    // Parcours des faces pour selectionner les valides (selectionnees selon le critere "strict")
    // Parcours les elmts. 
    // Parcours des faces de l'elmt. 
    // Si strict=0 et 1 face valide ou si strict=1 et toutes les faces valides, 
    // On ajoute l'elmt dans une nouvelle connectivite Elmt/Face
    // Au final, on reconstruit une connectivite complete en ajoutant l'ancienne connectivite Face/Noeuds 
    // Et on apelle cleanConnectivity
    // ----------------------------
    E_Float* tagp = f2->begin();   // champ indiquant le critere de selection
    E_Int* cnpp = cnp->begin();    // pointeur sur l ancienne connectivite
    E_Int nbFaces = cnpp[0];       // nombre de faces de l'ancienne connectivite
    E_Int sizeFN = cnpp[1];        // taille de l'ancienne connectivite Face/Noeuds
    E_Int nbElements = cnpp[sizeFN+2]; // nombre d'elts de l'ancienne connectivite
    E_Int sizeEF = cnpp[sizeFN+3];     // taille de l'ancienne connectivite Elmt/Faces
    E_Int* cnEFp = cnpp + 4+sizeFN;  // pointeur sur l'ancienne connectivite Elmt/Faces
    FldArrayI cn2 = FldArrayI(sizeEF); // nouvelle connectivite Elmt/Faces
    E_Int* cn2p = cn2.begin();       // pointeur sur la nouvelle connectivite
    E_Int size2 = 0;                 // compteur pour la nouvelle connectivite
    E_Int next = 0;                  // nbre d'elts selectionnes
    FldArrayI selectedFaces(nbFaces);  // tableau indiquant les faces valides 
    E_Int* selectedFacesp = selectedFaces.begin();
    E_Int nbfaces, nbnodes;
    E_Int isSel;
    
    E_Int newNumFace = 0;
    FldArrayI new_pg_ids; // Tableau d'indirection des faces (pour maj PE)
    FldArrayI keep_pg;    // Flag de conservation des faces 
    FldArrayI new_ph_ids; // Tableau d'indirection des elmts (pour maj PE)
    E_Int ii = 0 ;
 
    new_pg_ids = -1 ;
    new_ph_ids = -1 ;
    keep_pg    = -1 ;   

    // Selection des faces valides
    E_Int fa = 0; E_Int numFace = 0; cnpp += 2;
    

    if (strict == 0)  // cell selectionnee des qu'un sommet est tag=1
    {
      while (fa < sizeFN) // parcours de la connectivite face/noeuds
      {
        nbnodes = cnpp[0];
        selectedFacesp[numFace] = 0;
        for (E_Int n = 1; n <= nbnodes; n++)
        {
          if (tagp[cnpp[n]-1] >= oneEps) { selectedFacesp[numFace] = 1; break; }   
        }
        cnpp += nbnodes+1; fa += nbnodes+1; numFace++;
      }
    }
    else //strict=1, cell selectionnee si tous les sommets sont tag
    {
      while (fa < sizeFN) // parcours de la connectivite face/noeuds
      {
        nbnodes = cnpp[0];
        selectedFacesp[numFace] = 1;
        for (E_Int n = 1; n <= nbnodes; n++)
        {
          if (tagp[cnpp[n]-1] < oneEps) { selectedFacesp[numFace] = 0; break; }   
        }
        cnpp += nbnodes+1; fa += nbnodes+1; numFace++;
      }
    }

    
    // Si mise a jour du ParentElement, tab d'indirection des faces et des elmts
    // ------------------------------------------------------------------------
    if (PE != Py_None)
    {      
        new_pg_ids.malloc(nbFaces);    // Tableau d'indirection des faces (pour maj PE)
        keep_pg.malloc(nbFaces);       // Flag de conservation des faces 
        new_ph_ids.malloc(nbElements); // Tableau d'indirection des elmts (pour maj PE)

        new_pg_ids = -1;
        new_ph_ids = -1;
        keep_pg    = -1;
      
        // Selection des elements en fonction des faces valides
        if (strict == 0)  // cell selectionnee des qu'un sommet est tag=1
        {
          for (E_Int i = 0; i < nbElements; i++)
          {
            nbfaces       = cnEFp[0];
            for (E_Int n = 1; n <= nbfaces; n++)
            {
              if (selectedFacesp[cnEFp[n]-1] == 1) 
              { 
                cn2p[0] = nbfaces; size2 += 1;
                for (E_Int n = 1; n <= nbfaces; n++)
                {
                  cn2p[n] = cnEFp[n];
                  keep_pg[cnEFp[n]-1] = +1;
                }
		
                size2 += nbfaces; cn2p += nbfaces+1; next++;
                new_ph_ids[i] = ii;
                ii++;
		
                break; 
              }
            }
            cnEFp += nbfaces+1; 
          }
        }
        else //strict=1, cell selectionnee si tous les sommets sont tag
        {
          for (E_Int i = 0; i < nbElements; i++)
          {
            isSel = 0;
            nbfaces = cnEFp[0];
            for (E_Int n = 1; n <= nbfaces; n++)
            {
              if (selectedFacesp[cnEFp[n]-1] == 1) isSel++;
            }
            if (isSel == nbfaces) //cell selectionnee si tous les sommets sont tag
            {
              cn2p[0] = nbfaces; size2 +=1;
              for (E_Int n = 1; n <= nbfaces; n++)
	      {
		cn2p[n] = cnEFp[n];
		keep_pg[cnEFp[n]-1] = +1;
	      }
              size2 += nbfaces; cn2p += nbfaces+1; next++;
	      
	      new_ph_ids[i] = ii;	
	      ii++;
            }
            cnEFp += nbfaces+1; 
          } 
        }
        cn2.reAlloc(size2);	
	
	E_Int nn = 0 ; 
	for (E_Int n = 0; n<new_pg_ids.getSize(); n++)
	{
	  if (keep_pg[n]>0){ new_pg_ids[n] = nn; nn++; newNumFace++;}
	}

    }
    else  // PE == Py_None - pas de creation de tab d'indirection
    {
        // Selection des elements en fonction des faces valides
        if (strict == 0)  // cell selectionnee des qu'un sommet est tag=1
        {
          for (E_Int i = 0; i < nbElements; i++)
          {
            nbfaces = cnEFp[0];
            for (E_Int n = 1; n <= nbfaces; n++)
            {
              if (selectedFacesp[cnEFp[n]-1] == 1) 
              { 
                cn2p[0] = nbfaces; size2 += 1;
                for (E_Int n = 1; n <= nbfaces; n++) cn2p[n] = cnEFp[n];
                size2 += nbfaces; cn2p += nbfaces+1; next++;
                break; 
              }
            }
            cnEFp += nbfaces+1; 
          }
        }
        else //strict=1, cell selectionnee si tous les sommets sont tag
        {
          for (E_Int i = 0; i < nbElements; i++)
          {
            isSel = 0;
            nbfaces = cnEFp[0];
            for (E_Int n = 1; n <= nbfaces; n++)
            {
              if (selectedFacesp[cnEFp[n]-1] == 1) isSel++;
            }
            if (isSel == nbfaces) //cell selectionnee si tous les sommets sont tag
            {
              cn2p[0] = nbfaces; size2 +=1;
              for (E_Int n = 1; n <= nbfaces; n++) cn2p[n] = cnEFp[n];
              size2 += nbfaces; cn2p += nbfaces+1; next++;
            }
            cnEFp += nbfaces+1; 
          } 
        }
        cn2.reAlloc(size2);
    }
        
    // Cree la nouvelle connectivite complete
    E_Int coutsize = sizeFN+4+size2;
    FldArrayI* cout = new FldArrayI(coutsize);
    E_Int* coutp = cout->begin();
    cnpp = cnp->begin(); cn2p = cn2.begin();
    for (E_Int i = 0; i < sizeFN+2; i++) coutp[i] = cnpp[i];
    coutp += sizeFN+2;
    coutp[0] = next;
    coutp[1] = size2; coutp += 2;
    for (E_Int i = 0; i < size2; i++) coutp[i] = cn2p[i];

    RELEASESHAREDB(res2, tag, f2, cnp2);    

    if (PE != Py_None)
    {
      // Check numpy (parentElement)
      FldArrayI* cFE;
      E_Int res = K_NUMPY::getFromNumpyArray(PE, cFE, true);
          
      if (res == 0)
      {
	RELEASESHAREDN(PE, cFE);
	PyErr_SetString(PyExc_TypeError, 
			"selectCells: PE numpy is invalid.");
	return NULL;
      }
            
      ngon_t<K_FLD::FldArrayI> ng(*cout); // construction d'un ngon_t à partir d'un FldArrayI
      
      FldArrayI* cFEp_new = new FldArrayI(newNumFace,2);
      FldArrayI& cFE_new  = *cFEp_new ; 
      
      E_Int* cFEl = cFE_new.begin(1);
      E_Int* cFEr = cFE_new.begin(2);
      
      E_Int* cFEl_old = cFE->begin(1);
      E_Int* cFEr_old = cFE->begin(2);
      
      E_Int old_ph_1, old_ph_2;

      for (E_Int pgi = 0; pgi < nbFaces; pgi++)
      {
	      if (new_pg_ids[pgi]>=0)
	      {
	        old_ph_1 = cFEl_old[pgi]-1;
	        old_ph_2 = cFEr_old[pgi]-1;

	        if (old_ph_1 >= 0) // l'elmt gauche existe 
	        {
	          cFEl[new_pg_ids[pgi]] = new_ph_ids[old_ph_1]+1;
       
	          if (old_ph_2 >= 0) // l'elmt droit existe 
	            cFEr[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1; 
	          else
	            cFEr[new_pg_ids[pgi]] = 0;
	        }
	        else // l'elmt gauche a disparu - switch droite/gauche
	        {
	          cFEl[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1;
	          cFEr[new_pg_ids[pgi]] = 0;
	          // reverse
	          E_Int s = ng.PGs.stride(pgi);
	          E_Int* p = ng.PGs.get_facets_ptr(pgi);
	          std::reverse(p, p + s);
	        }
	        
	      }
      } // boucle pgi     
      
      // export ngon
      ng.export_to_array(*cout);

      // objet Python de sortie
      PyObject* pyPE = K_NUMPY::buildNumpyArray(cFE_new, 1);

      PyList_Append(l,pyPE);
      
      delete cFEp_new;
      RELEASESHAREDN(PE, cFE);
    }

    // close
    if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
      K_CONNECT::cleanConnectivityNGon(posx, posy, posz, 1.e-10, *fout, *cout);
    tpl = K_ARRAY::buildArray(*fout, varString, *cout, 8);
    delete fout; delete cout;
  }
  RELEASESHAREDB(res, array, f, cnp);
  
  PyList_Append(l,tpl) ; Py_DECREF(tpl);
  
  return l;  
}
