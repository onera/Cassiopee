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
PyObject* K_POST::selectCells(PyObject* self, PyObject* args)
{
  PyObject* arrayNodes; PyObject* tag;
  PyObject* arrayCenters = NULL;
  PyObject* PE;
  E_Int strict; E_Int cleanConnectivity;

  // Check different signatures
  Py_ssize_t nargs = PyTuple_GET_SIZE(args);
  if (nargs == 5)
  {
    if (!PYPARSETUPLE_(args, OO_ I_ O_ I_,
                       &arrayNodes, &tag, &strict, &PE, &cleanConnectivity))
      return NULL;
  }
  else
  {
    if (!PYPARSETUPLE_(args, OOO_ I_ O_ I_,
		                   &arrayNodes, &arrayCenters, &tag, &strict, &PE,
                       &cleanConnectivity))
      return NULL;
  }

  // Extract arrayNodes
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cnp;
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray3(arrayNodes, varString, f,
                                     ni, nj, nk, cnp, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: array is invalid.");
    return NULL;
  }

  // Extract arrayCenters 
  char* varStringc; char* eltTypec;
  FldArrayF* fc = NULL; FldArrayI* cnpc;
  E_Int resc = -1, nic = 0, njc = 0, nkc = 0, nfldc = 0;
  if (arrayCenters != NULL)
  {
    resc = K_ARRAY::getFromArray3(arrayCenters, varStringc, fc,
                                  nic, njc, nkc, cnpc, eltTypec);

    if (resc != 1 && resc != 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "selectCells: arrayCenters is invalid.");
      RELEASESHAREDB(res, arrayNodes, f, cnp);
      return NULL;
    }

    nfldc = fc->getNfld();  // nombre de champs en centres
  }

  // Extract tag
  char* varStringa; char* eltTypea;
  FldArrayF* fa; FldArrayI* cnpa;
  E_Int nia, nja, nka;
  E_Int resa = K_ARRAY::getFromArray3(tag, varStringa, fa,
                                      nia, nja, nka, cnpa, eltTypea);

  if (resa != 1 && resa != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag array is invalid.");
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    return NULL;
  }
  if (res != resa)
  {  
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag and array must represent the same grid.");
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    RELEASESHAREDB(resa, tag, fa, cnpa);
    return NULL;
  }

  E_Int api = f->getApi();
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();

  if (npts != fa->getSize())
  {
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    RELEASESHAREDB(resa, tag, fa, cnpa);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag and array must represent the same grid.");
    return NULL;
  }

  if (res == 2 && (cnp->getSize() != cnpa->getSize() ||
                   cnp->getNfld() != cnpa->getNfld())) 
  {
    RELEASESHAREDU(arrayNodes, f, cnp); 
    RELEASESHAREDB(resa, tag, fa, cnpa);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag and array must represent the same grid.");
    return NULL;  
  }

  if (fa->getNfld() != 1) 
  {
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    RELEASESHAREDB(resa, tag, fa, cnpa);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells: tag must have one variable only.");
    return NULL;
  }

  E_Float oneEps = 1.-1.e-10;
  // no check of coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;

  if (res == 1)
  {
    // Create BE connectivity
    E_Int dim0 = 0;
    if (ni > 1) dim0 += 1;
    if (nj > 1) dim0 += 1;
    if (nk > 1) dim0 += 1;
    E_Int nvpe;
    eltType = new char [K_ARRAY::VARSTRINGLENGTH];
    if (dim0 == 3) { strcpy(eltType, "HEXA"); nvpe = 8; }
    else if (dim0 == 2) { strcpy(eltType, "QUAD"); nvpe = 4; }
    else { strcpy(eltType, "BAR"); nvpe = 2; }

    E_Int ni1 = K_FUNC::E_max(1, ni-1);
    E_Int nj1 = K_FUNC::E_max(1, nj-1);
    E_Int nk1 = K_FUNC::E_max(1, nk-1);
    E_Int ninj = ni*nj;
    E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees

    cnp = new FldArrayI(ncells, nvpe);
    FldArrayI& cm = *(cnp->getConnect(0));

    if (dim0 == 1)
    {
      E_Int ind1, ind2;
      if (nj1 == 1 && nk1 == 1)
      {
        for (E_Int i = 0; i < ni1; i++)
        {
          ind1 = i + 1;     // (i,1,1)
          ind2 = ind1 + 1;  // (i+1,1,1)
          cm(i,1) = ind1; cm(i,2) = ind2;
        }
      }
      else if (ni1 == 1 && nj1 == 1)
      {
        for (E_Int k = 0; k < nk1; k++)
        {
          ind1 = 1 + k*ninj;   // (1,1,k)
          ind2 = ind1 + ninj;  // (1,1,k+1)
          cm(k,1) = ind1; cm(k,2) = ind2;
        }
      }
      else  // ni1 == 1 && nk1 == 1
      {
        for (E_Int j = 0; j < nj1; j++)
        {
          ind1 = j*ni + 1;   // (1,j,1)
          ind2 = ind1 + ni;  // (1,j+1,1)
          cm(j,1) = ind1; cm(j,2) = ind2;
        }
      }
    }
    else if (dim0 == 2)
    {
      if (nk1 == 1)
      {
        #pragma omp parallel if (nj1 > __MIN_SIZE_MEAN__)
        {
          E_Int c, ind1, ind2, ind3, ind4;
          #pragma omp for collapse(2)
          for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = i + j*ni + 1; // (i,j,1)
            ind2 = ind1 + 1;     // (i+1,j,1)
            ind3 = ind2 + ni;    // (i+1,j+1,1)
            ind4 = ind3 - 1;     // (i,j+1,1)
            c = i + j*ni1;
            cm(c,1) = ind1; cm(c,2) = ind2;
            cm(c,3) = ind3; cm(c,4) = ind4;
          }
        }
      }
      else if (nj1 == 1)
      {
        #pragma omp parallel if (nk1 > __MIN_SIZE_MEAN__)
        {
          E_Int c, ind1, ind2, ind3, ind4;
          #pragma omp for collapse(2)
          for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = i + k*ninj + 1;  // (i,1,k)
            ind2 = ind1 + ninj;     // (i,1,k+1)
            ind3 = ind2 + 1;        // (i+1,1,k+1)
            ind4 = ind3 - 1;        // (i,1,k+1)
            c = i + k*ni1;
            cm(c,1) = ind1; cm(c,2) = ind2;
            cm(c,3) = ind3; cm(c,4) = ind4;
          }
        }
      }
      else // i1 = 1 
      {
        #pragma omp parallel if (nk1 > __MIN_SIZE_MEAN__)
        {
          E_Int c, ind1, ind2, ind3, ind4;
          #pragma omp for collapse(2)
          for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj1; j++)
          {
            ind1 = 1 + j*ni + k*ninj; // (1,j,k)
            ind2 = ind1 + ni;         // (1,j+1,k)
            ind3 = ind2 + ninj;       // (1,j+1,k+1)
            ind4 = ind3 - ni;         // (1,j,k+1)
            c = j+k*nj1;
            cm(c,1) = ind1; cm(c,2) = ind2;
            cm(c,3) = ind3; cm(c,4) = ind4;
          }
        }
      }
    }
    else  // dim0 = 3
    {
      #pragma omp parallel if (nk1 > __MIN_SIZE_MEAN__)
      {
        E_Int c, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
        #pragma omp for collapse(3)
        for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          ind1 = 1 + i + j*ni + k*ninj; // A(  i,  j,k)
          ind2 = ind1 + 1;              // B(i+1,  j,k)
          ind3 = ind2 + ni;             // C(i+1,j+1,k)
          ind4 = ind3 - 1;              // D(  i,j+1,k)
          ind5 = ind1 + ninj;           // E(  i,  j,k+1)
          ind6 = ind2 + ninj;           // F(i+1,  j,k+1)
          ind7 = ind3 + ninj;           // G(i+1,j+1,k+1)
          ind8 = ind4 + ninj;           // H(  i,j+1,k+1) 
          c = i+j*ni1+k*ni1*nj1;
          cm(c,1) = ind1; cm(c,2) = ind2;
          cm(c,3) = ind3; cm(c,4) = ind4;
          cm(c,5) = ind5; cm(c,6) = ind6;
          cm(c,7) = ind7; cm(c,8) = ind8;
        }
      }
    }
  }

  // Selection
  PyObject* l = PyList_New(0);
  PyObject* tpln; PyObject* tplc = NULL;

  if (K_STRING::cmp(eltType, "NGON") == 0)
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
    E_Float* tagp = fa->begin();   // champ indiquant le critere de selection
    E_Int* cnpp = cnp->begin();    // pointeur sur l ancienne connectivite
    E_Int nbFaces = cnpp[0];       // nombre de faces de l'ancienne connectivite
    E_Int sizeFN = cnpp[1];        // taille de l'ancienne connectivite Face/Noeuds
    E_Int nbElements = cnpp[sizeFN+2]; // nombre d'elts de l'ancienne connectivite
    E_Int sizeEF = cnpp[sizeFN+3];     // taille de l'ancienne connectivite Elmt/Faces
    E_Int* cnEFp = cnpp+4+sizeFN;  // pointeur sur l'ancienne connectivite Elmt/Faces
    FldArrayI cn2 = FldArrayI(sizeEF); // nouvelle connectivite Elmt/Faces
    E_Int* cn2p = cn2.begin();       // pointeur sur la nouvelle connectivite
    E_Int size2 = 0;                 // compteur pour la nouvelle connectivite
    E_Int next = 0;                  // nbre d'elts selectionnes
    FldArrayI selectedFaces(nbFaces);  // tableau indiquant les faces valides 
    E_Int* selectedFacesp = selectedFaces.begin();
    E_Int nbfaces, nbnodes;
    E_Int cnt;

    // Champs en centre
    FldArrayF* foutc = NULL;
    if (arrayCenters != NULL) foutc = new FldArrayF(*fc);
    E_Int ii = 0;

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
    FldArrayI new_pg_ids; // Tableau d'indirection des faces (pour maj PE)
    FldArrayI keep_pg;    // Flag de conservation des faces 
    FldArrayI new_ph_ids; // Tableau d'indirection des elmts (pour maj PE)
    new_pg_ids = -1;
    new_ph_ids = -1;
    keep_pg = -1;

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

              for (E_Int k = 1; k <= nfldc; k++) (*foutc)(ii,k) = (*fc)(i,k);

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
          cnt = 0;
          new_ph_ids[i] = -1;
          nbfaces = cnEFp[0];
          for (E_Int n = 1; n <= nbfaces; n++)
          {
            if (selectedFacesp[cnEFp[n]-1] == 1) cnt++;
          }
          if (cnt == nbfaces) //cell selectionnee si tous les sommets sont tag
          {
            cn2p[0] = nbfaces; size2 +=1;
            for (E_Int n = 1; n <= nbfaces; n++)
            {
              cn2p[n] = cnEFp[n];
              keep_pg[cnEFp[n]-1] = +1;
            }
            size2 += nbfaces; cn2p += nbfaces+1; next++;
      
            for (E_Int k = 1; k <= nfldc; k++) (*foutc)(ii,k) = (*fc)(i,k);
  
            new_ph_ids[i] = ii;
            ii++;
          }
          cnEFp += nbfaces+1;
        }
      }

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

              for (E_Int k = 1; k <= nfldc; k++) (*foutc)(ii,k) = (*fc)(i,k);
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
          cnt = 0;
          nbfaces = cnEFp[0];
          for (E_Int n = 1; n <= nbfaces; n++)
          {
            if (selectedFacesp[cnEFp[n]-1] == 1) cnt++;
          }
          if (cnt == nbfaces) //cell selectionnee si tous les sommets sont tag
          {
            cn2p[0] = nbfaces; size2 +=1;
            for (E_Int n = 1; n <= nbfaces; n++) cn2p[n] = cnEFp[n];
            size2 += nbfaces; cn2p += nbfaces+1; next++;

            for (E_Int k = 1; k <= nfldc; k++) (*foutc)(ii,k) = (*fc)(i,k);
            ii++;
          }
          cnEFp += nbfaces+1; 
        } 
      }
    }

    cn2.reAlloc(size2);
    if (arrayCenters != NULL) foutc->reAllocMat(ii, nfldc);
    
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

    if (PE != Py_None)
    {
      // Check numpy (parentElement)
      FldArrayI* cFE;
      E_Int res = K_NUMPY::getFromNumpyArray(PE, cFE);
      
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
            {
              cFEr[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1;
            }
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

      PyList_Append(l, pyPE);
      
      delete cFEp_new;
      RELEASESHAREDN(PE, cFE);
    }

    // close
    if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
      K_CONNECT::cleanConnectivityNGon(posx, posy, posz, 1.e-10, *fout, *cout);

    cout->setNGonType(1);
    tpln = K_ARRAY::buildArray3(*fout, varString, *cout, eltType, api);
    if (arrayCenters != NULL)
    {
      tplc = K_ARRAY::buildArray3(*foutc, varStringc, *cout, eltType, api);
      delete foutc;
    }
    delete fout; delete cout;
  }
  else if (K_STRING::cmp(eltType, "NODE") == 0)
  {
    FldArrayF* an = new FldArrayF();
    FldArrayF& coord = *an;

    // Selection des vertex
    FldArrayI selected(npts, 2);
    selected.setAllValuesAtNull();
    E_Int* selectedp = selected.begin();
    E_Float* tagp = fa->begin();
    E_Int count = 0;
    for (E_Int i = 0; i < npts; i++)
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
      for (E_Int i = 0; i < npts; i++)
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
  
    if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, *an, *acn);
    tpln = K_ARRAY::buildArray3(*an, varString, *acn, eltType, api);
    if (arrayCenters != NULL) tplc = tpln;
    delete an; delete acn;
  }
  else  // other BE/ME
  {
    E_Float* tagp = fa->begin();

    // Number of elements per connectivity and cumulative offsets
    E_Int nc = cnp->getNConnect();
    std::vector<E_Int> nepc(nc);
    std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cnp->getConnect(ic));
      nepc[ic] = cm.getSize();
      cumnepc[ic+1] = cumnepc[ic] + nepc[ic];
    }
    E_Int ntotElts = cumnepc[nc];

    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);

    // In a first pass, tag selected elements and their vertices
    std::vector<E_Int> vindir(npts, 0);
    std::vector<E_Int> eindir(ntotElts, 0);

    #pragma omp parallel
    {
      E_Int indv, nvpe, eidx;
      E_Bool selected;

      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(cnp->getConnect(ic));
        nvpe = cm.getNfld();

        if (strict == 0)  // selected if at least one vertex is tagged with 1
        {
          #pragma omp for nowait schedule(static)
          for (E_Int i = 0; i < nepc[ic]; i++)
          {
            selected = false;
            for (E_Int j = 1; j <= nvpe; j++)
            {
              indv = cm(i,j) - 1;
              if (tagp[indv] >= oneEps) { selected = true; break; }
            }

            if (selected)
            {
              // cell selected, tag all vertices
              eidx = cumnepc[ic] + i;
              eindir[eidx] = 1;
              for (E_Int j = 1; j <= nvpe; j++)
              {
                indv = cm(i,j) - 1;
                vindir[indv] = 1;
              }
            }
          }
        }
        else  // selected if all vertices are tagged with 1
        {
          #pragma omp for nowait schedule(static)
          for (E_Int i = 0; i < nepc[ic]; i++)
          {
            selected = true;
            for (E_Int j = 1; j <= nvpe; j++)
            {
              indv = cm(i,j) - 1;
              if (tagp[indv] < oneEps) { selected = false; break; }
            }

            if (selected)
            {
              // cell selected, tag all vertices
              eidx = cumnepc[ic] + i;
              eindir[eidx] = 1;
              for (E_Int j = 1; j <= nvpe; j++)
              {
                indv = cm(i,j) - 1;
                vindir[indv] = 1;
              }
            }
          }
        }
      }
    }

    // Transform the masks into maps from old to new numbering.
    // eindir is renumbered per connectivity (bucket), vindir globally
    std::vector<E_Int> tmp_nepc2 = K_CONNECT::mask2Indir(eindir, nepc);
    E_Int npts2 = K_CONNECT::mask2Indir(vindir);

    // Nothing selected: return an empty NODE connectivity
    if (npts2 == 0)
    {
      tpln = K_ARRAY::buildArray3(nfld, varString, 0, 0, "NODE", false, api);
      if (arrayCenters != NULL) tplc = tpln;
    }
    else
    {
      // Build new eltType and nepc2 from conns with at least one element
      // ('tmp_' is uncompressed: same number of connectivities as the input)
      E_Int nc2 = 0;
      char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH];
      eltType2[0] = '\0';
      std::vector<E_Int> nepc2;
      for (E_Int ic = 0; ic < nc; ic++)
      {
        if (tmp_nepc2[ic] <= 0) continue;
        nc2++;
        nepc2.push_back(tmp_nepc2[ic]);
        if (eltType2[0] == '\0') strcpy(eltType2, eltTypes[ic]);
        else { strcat(eltType2, ","); strcat(eltType2, eltTypes[ic]); }
      }

      // Cumulative number of selected elements, uncompressed (offset of each
      // connectivity in the output center fields)
      std::vector<E_Int> cumnepc2(nc+1); cumnepc2[0] = 0;
      for (E_Int ic = 0; ic < nc; ic++)
        cumnepc2[ic+1] = cumnepc2[ic] + tmp_nepc2[ic];

      // Build new ME connectivity
      tpln = K_ARRAY::buildArray3(nfld, varString, npts2, nepc2,
                                  eltType2, false, api);
      FldArrayF* f2; FldArrayI* cn2;
      K_ARRAY::getFromArray3(tpln, f2, cn2);

      #pragma omp parallel
      {
        E_Int indv, nvpe, ind, ic2;

        // Copy fields at nodes
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* fp = f->begin(n);
          E_Float* f2p = f2->begin(n);
          #pragma omp for nowait
          for (E_Int i = 0; i < npts; i++)
          {
            indv = vindir[i];
            if (indv > 0) f2p[indv-1] = fp[i];
          }
        }

        // Connectivity
        ic2 = 0;
        for (E_Int ic = 0; ic < nc; ic++)
        {
          if (tmp_nepc2[ic] == 0) continue;  // no selected elements in this conn
          FldArrayI& cm = *(cnp->getConnect(ic));
          FldArrayI& cm2 = *(cn2->getConnect(ic2));
          nvpe = cm.getNfld();

          #pragma omp for nowait schedule(static)
          for (E_Int i = 0; i < nepc[ic]; i++)
          {
            ind = eindir[cumnepc[ic] + i];
            if (ind > 0)
            {
              for (E_Int j = 1; j <= nvpe; j++)
              {
                indv = cm(i,j) - 1;
                cm2(ind-1, j) = vindir[indv];
              }
            }
          }
          ic2++;
        }
      }

      // Free memory
      std::vector<E_Int>().swap(vindir);

      FldArrayF* fc2 = NULL;
      if (arrayCenters != NULL)
      {
        tplc = K_ARRAY::buildArray3(nfldc, varStringc, npts2,
                                    *cn2, eltType2, true, api, true);
        K_ARRAY::getFromArray3(tplc, fc2);

        // Copy fields at centers
        #pragma omp parallel
        {
          E_Int eidx, ind;
          for (E_Int n = 1; n <= nfldc; n++)
          {
            E_Float* fcp = fc->begin(n);
            E_Float* fc2p = fc2->begin(n);
            for (E_Int ic = 0; ic < nc; ic++)
            {
              #pragma omp for nowait
              for (E_Int i = 0; i < nepc[ic]; i++)
              {
                eidx = cumnepc[ic] + i;
                ind = eindir[eidx];
                if (ind > 0) fc2p[cumnepc2[ic] + ind-1] = fcp[eidx];
              }
            }
          }
        }
      }

      for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
      if (res == 1) { delete cnp; delete[] eltType; }
      if (arrayCenters != NULL) RELEASESHAREDS(tplc, fc2);
    }
  }

  RELEASESHAREDB(resa, tag, fa, cnpa);
  RELEASESHAREDB(res, arrayNodes, f, cnp);
  PyList_Append(l, tpln); Py_DECREF(tpln);
  if (arrayCenters != NULL)
  {
    RELEASESHAREDB(resc, arrayCenters, fc, cnpc);
    PyList_Append(l, tplc); Py_DECREF(tplc);
  }
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
  if (!PYPARSETUPLE_(args, O_ I_, &tag, &flag)) return NULL;

  // Extract tag (centers)
  char* varStringa; char* eltTypea;
  FldArrayF* fa; FldArrayI* cna;
  E_Int resa, nia, nja, nka;

  if (flag == 0) // array
  {
    resa = K_ARRAY::getFromArray3(tag, varStringa, fa,
                                  nia, nja, nka, cna, eltTypea);
    if (resa != 1 && resa != 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "selectCells3: tag array is invalid.");
      return NULL;
    }
  }
  else // numpy
  {
    resa = K_NUMPY::getFromNumpyArray(tag, fa);
    if (resa == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "selectCells3: tag numpy is invalid.");
      return NULL;
    }
  }

  E_Int nelts = fa->getSize();
  E_Float* tagp = fa->begin(); 
  E_Float oneEps = 1.-1.e-10;

  // Compte les tag > 1-eps
  E_Int nthreads = __NUMTHREADS__;
  E_Int* nb = new E_Int [nthreads];
  E_Int* pos = new E_Int [nthreads];
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

  delete[] nb; delete[] pos;

  if (flag == 0) { RELEASESHAREDB(resa, tag, fa, cna); }
  else RELEASESHAREDN(tag, fa);
  return o;
}
