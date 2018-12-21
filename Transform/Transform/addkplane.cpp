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
# include "transform.h"
 
using namespace K_FLD;

// ============================================================================
/* Ajoute un plan k a un maillage */
// ============================================================================
PyObject* K_TRANSFORM::addkplane(PyObject* self, PyObject* args)
{
  PyObject* array; 
  if (!PyArg_ParseTuple(args, "O", &array))
  {
      return NULL;
  }
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, 
                                    cn, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: unknown type of array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;
  
  // Vector add
  E_Float vx = 0.; E_Float vy = 0.; E_Float vz = 1.;

  if (res == 1)
  {
    E_Int imjm, imjmkm, km1;
    imjm = im*jm;
    imjmkm = imjm*km;
    km1 = km+1;
    E_Int npts = imjm*km1;
    E_Int nfld = f->getNfld();
    PyObject* tpl = K_ARRAY::buildArray(nfld, varString, im, jm, km1);
    E_Float* nzp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF nz(npts, nfld, nzp, true);
    
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* newzonep = nz.begin(n);
      E_Float* fpn = f->begin(n);
#pragma omp parallel
      {
        E_Int ind, ind2;
#pragma omp for
        for (E_Int i = 0; i < imjmkm; i++)
        {
          newzonep[i] = fpn[i];
        }
      
#pragma omp for
        for (E_Int i = 0; i < im*jm; i++)
        {
          ind   = i + (km-1)*imjm;
          ind2  = i + imjmkm; 
          newzonep[ind2] = fpn[ind];
        }
      }
    }

    if (posx > 0 && posy > 0 && posz > 0)
    {
      E_Float* nfx = nz.begin(posx);
      E_Float* nfy = nz.begin(posy);
      E_Float* nfz = nz.begin(posz);
      E_Float* fx = f->begin(posx);
      E_Float* fy = f->begin(posy);
      E_Float* fz = f->begin(posz);
      for (E_Int j = 0; j < jm; j++)
#pragma omp parallel
      {
        E_Int ind, ind2;
#pragma omp for
        for (E_Int i = 0; i < im; i++)
        {
          ind = i + j*im + (km-1)*imjm;
          ind2  = i + j*im + imjmkm; 
          nfx[ind2] = fx[ind] + vx;
          nfy[ind2] = fy[ind] + vy;
          nfz[ind2] = fz[ind] + vz;
        }
      }
    }

    RELEASESHAREDS(array,f);
    return tpl;
  }
  else if (res == 2)
  {
    if (K_STRING::cmp(eltType, "BAR") !=0 && 
        K_STRING::cmp(eltType, "QUAD")!=0 && 
        K_STRING::cmp(eltType, "TRI") !=0 &&
        K_STRING::cmp(eltType,"NGON") !=0)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "addkplane: only for BAR, QUAD, TRI, NGON or struct-arrays.");
      return NULL;
    }
    
    // Duplication des coordonnees
    char elt[16]; 
    if (K_STRING::cmp(eltType, "BAR") == 0) // -> QUAD
      strcpy(elt, "QUAD");
    else if (K_STRING::cmp(eltType, "QUAD") == 0) // -> HEXA
      strcpy(elt, "HEXA");
    else if (K_STRING::cmp(eltType, "TRI") == 0) // -> PENTA
      strcpy(elt, "PENTA");
    else if (K_STRING::cmp(eltType, "NGON") == 0) // -> NGON
      strcpy(elt, "NGON");
    PyObject* tpl;
    if (K_STRING::cmp(eltType, "NGON") != 0) // Elements basiques
    { 
      E_Int np = f->getSize(); E_Int npts = 2*np;
      E_Int nfld = f->getNfld();
      E_Int ne = cn->getSize(); E_Int nvert = 2* cn->getNfld();
      tpl = K_ARRAY::buildArray(nfld, varString, npts, ne, -1, elt);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI connect(ne, nvert, cnnp, true);
      E_Float* nzp = K_ARRAY::getFieldPtr(tpl);
      FldArrayF nz(npts, nfld, nzp, true);
      
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* newzonep = nz.begin(n);
        E_Float* fp = f->begin(n);
        for (E_Int i = 0; i < np; i++) newzonep[i] = fp[i];
        for (E_Int i = 0; i < np; i++) newzonep[i+np] = fp[i];
      }
      if (posx > 0 && posy >0 && posz > 0)
      {
        E_Float* nfx = nz.begin(posx);
        E_Float* nfy = nz.begin(posy);
        E_Float* nfz = nz.begin(posz);
        E_Float* fx = f->begin(posx);
        E_Float* fy = f->begin(posy);
        E_Float* fz = f->begin(posz);
        for (E_Int i = 0; i < np; i++)
        {
          nfx[i+np] = fx[i] + vx;
          nfy[i+np] = fy[i] + vy;
          nfz[i+np] = fz[i] + vz;
        }
      }
      
      if (nvert == 4) // -> QUAD
      {
        E_Int* cn1 = connect.begin(1);
        E_Int* cn2 = connect.begin(2);
        E_Int* cn3 = connect.begin(3);
        E_Int* cn4 = connect.begin(4);
        E_Int* cnp1 = cn->begin(1);
        E_Int* cnp2 = cn->begin(2);

        for (E_Int i = 0; i < ne; i++)
        {
          cn1[i] = cnp1[i];
          cn2[i] = cnp2[i];
          cn3[i] = cnp2[i]+np;
          cn4[i] = cnp1[i]+np;
        }
      }
      else if (nvert == 8) // -> HEXA
      {
        E_Int* cn1 = connect.begin(1);
        E_Int* cn2 = connect.begin(2);
        E_Int* cn3 = connect.begin(3);
        E_Int* cn4 = connect.begin(4);
        E_Int* cn5 = connect.begin(5);
        E_Int* cn6 = connect.begin(6);
        E_Int* cn7 = connect.begin(7);
        E_Int* cn8 = connect.begin(8);
        E_Int* cnp1 = cn->begin(1);
        E_Int* cnp2 = cn->begin(2);
        E_Int* cnp3 = cn->begin(3);
        E_Int* cnp4 = cn->begin(4);
        
        for (E_Int i = 0; i < ne; i++)
        {
          cn1[i] = cnp1[i];
          cn2[i] = cnp2[i];
          cn3[i] = cnp3[i];
          cn4[i] = cnp4[i];
          cn5[i] = cnp1[i]+np;
          cn6[i] = cnp2[i]+np;
          cn7[i] = cnp3[i]+np;
          cn8[i] = cnp4[i]+np;
        }
      }
      else if (nvert == 6) // -> PENTA
      {
        E_Int* cn1 = connect.begin(1);
        E_Int* cn2 = connect.begin(2);
        E_Int* cn3 = connect.begin(3);
        E_Int* cn4 = connect.begin(4);
        E_Int* cn5 = connect.begin(5);
        E_Int* cn6 = connect.begin(6);
        E_Int* cnp1 = cn->begin(1);
        E_Int* cnp2 = cn->begin(2);
        E_Int* cnp3 = cn->begin(3);
        
        for (E_Int i = 0; i < ne; i++)
        {
          cn1[i] = cnp1[i];
          cn2[i] = cnp2[i];
          cn3[i] = cnp3[i];
          cn4[i] = cnp1[i]+np;
          cn5[i] = cnp2[i]+np;
          cn6[i] = cnp3[i]+np;
        }
      }
    }//elts basiques
    else  // NGONs
    {
      E_Int nps = f->getSize(); E_Int npv = 2*nps; E_Int nfld = f->getNfld();
      E_Int* cnsp = cn->begin(); // pointeur sur la connectivite NGON surfacique      
      E_Int nfs = cnsp[0];
      E_Int sizeFNs = cnsp[1]; //  taille de la connectivite Face/Noeuds
      E_Int sizeEFs = cnsp[3+sizeFNs]; //  taille de la connectivite Elts/Faces
      E_Int nes = cnsp[sizeFNs+2];  // nombre total d elements
      FldArrayI posEltsSurf; K_CONNECT::getPosElts(*cn, posEltsSurf);
      FldArrayI posFacesSurf; K_CONNECT::getPosFaces(*cn, posFacesSurf);
      // on verifie que le NGON est surfacique a partir de la premiere face
      if (cnsp[2] != 2) // la face a plus de 2 sommets ce n'est pas une arete
      {
        PyErr_SetString(PyExc_TypeError,
                        "addkplane: NGON array must be a surface.");
        RELEASESHAREDU(array, f, cn); return NULL;
      }
     
      E_Int sizeEFv = sizeEFs+nes*2;// (nfacess+2) faces dans le volume
      //E_Int nev = nes; // nb d elts dans le NGON volumique
      E_Int nfv = nfs + 2*nes;//nb de faces ds le NGON volumique
      E_Int sumFS = 0;// dimensionnement du tableau faces/noeuds pour les faces correspondant aux elts surfaciques
      E_Int* ptr = cnsp+sizeFNs+4;
      E_Int e = 0;
      while (e < nes)
      {
        E_Int nfloc = ptr[0];
        sumFS += nfloc+1;// pour chq face vol : nfacesloc vertex + 1 pour dimensionner
        ptr += nfloc+1;
        e++;
      }
      E_Int sizeFNv = nfs*(4+1) + 2*sumFS;// (nb de sommets + 1)
      E_Int csize = sizeEFv+sizeFNv+4;
      tpl = K_ARRAY::buildArray(nfld, varString, npv, nes, -1, "NGON", false, csize);
      E_Float* nzp = K_ARRAY::getFieldPtr(tpl);
      FldArrayF nz(npv, nfld, nzp, true);
      E_Int* cnvp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cnv(csize, 1, cnvp, true);
      cnvp[0] = nfv;
      cnvp[1] = sizeFNv;
      cnvp[sizeFNv+2] = nes;
      cnvp[sizeFNv+3] = sizeEFv;
      // duplication des champs, avec z = z+1 
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* newzonep = nz.begin(n);
        E_Float* fp = f->begin(n);
        for (E_Int i = 0; i < nps; i++) newzonep[i] = fp[i];
        for (E_Int i = 0; i < nps; i++) newzonep[i+nps] = fp[i];
      }
      if (posx > 0 && posy > 0 && posz > 0)
      {
        E_Float* nfx = nz.begin(posx);
        E_Float* nfy = nz.begin(posy);
        E_Float* nfz = nz.begin(posz);
        E_Float* fx = f->begin(posx);
        E_Float* fy = f->begin(posy);
        E_Float* fz = f->begin(posz);
        for (E_Int i = 0; i < nps; i++)
        {
          nfx[i+nps] = fx[i] + vx;
          nfy[i+nps] = fy[i] + vy;
          nfz[i+nps] = fz[i] + vz;
        }
      }
      //=======================================================================
      // connectivites
      //=======================================================================
      E_Int* ptrFNv = cnvp+2;//ptr cFN vol
      E_Int* ptrFNs = cnsp+2;//ptr cFN surf
      // a partir de chq face construction des faces laterales "quad" 
      // extrudee a partir des faces surfaciques
      E_Int nofv = 0;
      while (nofv < nfs)
      {
        ptrFNv[0] = 4;
        ptrFNv[1] = ptrFNs[1];
        ptrFNv[2] = ptrFNs[2];
        ptrFNv[3] = ptrFNs[2]+nps;
        ptrFNv[4] = ptrFNs[1]+nps;
        ptrFNv += 5; ptrFNs += 3;
        nofv++;
      }     
      // a partir des elts: recup des faces laterales: meme numerotation 
      // qu'en surfacique
      E_Int noe = 0;
      E_Int* ptrEFv = cnvp+sizeFNv+4;//ptr cEF vol
      E_Int* ptrEFs = cnsp+sizeFNs+4;//ptr cEF surf
      while (noe < nes)
      {
        E_Int nfacessloc = ptrEFs[0];
        ptrEFv[0] = nfacessloc+2;
        for (E_Int nof = 1; nof <= nfacessloc; nof++)
          ptrEFv[nof] = ptrEFs[nof];
        noe++;
        ptrEFs += nfacessloc+1;
        ptrEFv += nfacessloc+3;
      }
      
      // construction des faces NGons
      noe = 0;
      ptrEFv = cnvp+sizeFNv+4;//ptr cEF vol
      ptrEFs = cnsp+sizeFNs+4;//ptr cEF surf
      std::vector<E_Int> indices;
      while (noe < nes) 
      {
        // les vertex surfaciques sont dans l'ordre rotatif
        indices.clear();
        K_CONNECT::getVertexIndices(cn->begin(), posFacesSurf.begin(), posEltsSurf[noe], indices);
        E_Int nvert = indices.size();

        //creation de la face correspondant a l elt surfacique
        ptrFNv[0] = nvert;
        for (E_Int i = 0; i < nvert; i++) ptrFNv[i+1] = indices[i];        
        ptrFNv+= nvert+1; 
        //creation de la face shiftee en z+1
        ptrFNv[0] = nvert;
        for (E_Int i = 0; i < nvert; i++) ptrFNv[i+1] = indices[i]+nps;        
        ptrFNv+= nvert+1;
 
        //modif de l'elt: on remplit les 2 derniers
        E_Int nfacesV = ptrEFv[0];
        ptrEFv[nfacesV-1] = nofv+1;
        ptrEFv[nfacesV]   = nofv+2;
        noe++; nofv += 2; ptrEFv += nfacesV+1;
      }    
    }//NGONs
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: unknow type of array.");
    return NULL;
  }
}

//=============================================================================
/* copie les champs de arrayC dans les centres de arrayK 
   arrayC doit definir les champs en centres de l'array arrayK avant 
   son addkplane(N) 
   On fournit arrayK et pas seulement N car on a besoin de connaitre la
   nature de l array (a cause du passage de 2D a 3D si N>1) */
//=============================================================================
PyObject* K_TRANSFORM::addkplaneCenters(PyObject* self, PyObject* args)
{
  PyObject *arrayC, *arrayK;
  E_Int N;
  if (!PYPARSETUPLEI(args,
                    "OOl", "OOi",
                    &arrayC, &arrayK, &N))
  {
      return NULL;
  }

  // Check array of centers
  E_Int imc, jmc, kmc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray(arrayC, varStringc, fc, imc, jmc, kmc, 
                                     cnc, eltTypec, true);

  if (resc != 1 && resc != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: unknown type of array for centers.");
    return NULL;
  }
  // Check array of nodes already with N kplanes
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(arrayK, varString, f, im, jm, km, 
                                    cn, eltType, true);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: unknown type of array for nodes.");
    return NULL;
  }
  if (res != resc)
  {
    RELEASESHAREDB(resc, arrayC, fc, cnc);
    RELEASESHAREDB(res, arrayK, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: array of centers and nodes must be both structured or unstructured.");
    return NULL;
  }
  if (resc == 2) 
  {
    if (K_STRING::cmp(eltTypec,"BAR*") != 0 && 
        K_STRING::cmp(eltTypec,"TRI*") != 0 && 
        K_STRING::cmp(eltTypec,"QUAD*") != 0 && 
        K_STRING::cmp(eltTypec,"NGON*") != 0)
    {
      RELEASESHAREDB(resc, arrayC, fc, cnc);
      RELEASESHAREDB(res, arrayK, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "addkplane: array of centers and nodes must be both structured or unstructured.");
      return NULL;
    }
  }
  E_Int nfld = fc->getNfld();
  
  PyObject* tpl;
  E_Float val;
  if (resc == 1) 
  {
    E_Int size = imc*jmc*kmc;
    E_Int imcjmc = imc*jmc;
    tpl = K_ARRAY::buildArray(nfld, varStringc, imc, jmc, km-1);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF field(imc*jmc*(km-1),nfld,fp,true);
    if (km == 2) // cas 2D
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* ptrFc = fc->begin(n);
        E_Float* ptrF = field.begin(n);
        #pragma omp parallel for
        for (E_Int ind = 0; ind < size; ind++) ptrF[ind] = ptrFc[ind];
      }
    }
    else 
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* ptrFc = fc->begin(n);
        E_Float* ptrF = field.begin(n);
        for (E_Int noz = 0; noz < km-1; noz++)
          #pragma omp parallel for
          for (E_Int ind = 0; ind < imcjmc; ind++)
          {
            val = ptrFc[ind];
            ptrF[ind] = val;      
            ptrF[ind+noz*imcjmc] = val;
          }
      }
    }
  }
  else 
  {
    E_Int csize = cn->getSize()*cn->getNfld();    
    E_Int* cncp = cnc->begin();
    E_Int* cnp0 = cn->begin();
    E_Int nelts, nelts0;
    if (K_STRING::cmp(eltType,"NGON") != 0) nelts = cn->getSize();
    else nelts = cnp0[2+cnp0[1]];
    if (K_STRING::cmp(eltTypec,"NGON*") != 0) nelts0 = cnc->getSize();
    else nelts0 = cncp[2+cncp[1]];

    tpl = K_ARRAY::buildArray(nfld, varStringc, nelts, nelts,
                              -1, eltType, true, csize);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fn(nelts, nfld, fp, true); 
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayI cnn(cn->getSize(), cn->getNfld(), cnnp, true); cnn = *cn;

    if (K_STRING::cmp(eltTypec,"BAR*") == 0 || 
        K_STRING::cmp(eltTypec,"TRI*") == 0 || 
        K_STRING::cmp(eltTypec,"QUAD*") == 0) // addkplane: QUAD*
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* ptrFc = fc->begin(n);
        E_Float* ptrF = fn.begin(n);
        for (E_Int noz = 0; noz < N; noz++)
          for (E_Int ind = 0; ind < nelts0; ind++)         
            ptrF[ind+noz*nelts0] =  ptrFc[ind];      
      }
    }
    else if (K_STRING::cmp(eltTypec,"NGON*") == 0) 
    {
      E_Int* cncp = cnc->begin();
      // il faut verifier que le NGON* est 2D pour extruder
      if (cncp[2] != 2) // la face a plus de 2 sommets ce n'est pas une arete
      {
        PyErr_SetString(PyExc_TypeError,
                        "addkplane: NGON array must be a surface.");
        RELEASESHAREDU(arrayC, fc, cnc);
        RELEASESHAREDU(arrayK, f, cn);
        return NULL;
      }
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* ptrFc = fc->begin(n);
        E_Float* ptrF = fn.begin(n);
        for (E_Int noz = 0; noz < N; noz++)
          for (E_Int ind = 0; ind < nelts0; ind++)         
            ptrF[ind+noz*nelts0] =  ptrFc[ind];      
      }
    }    
  }
  RELEASESHAREDB(resc, arrayC, fc, cnc);
  RELEASESHAREDB(res, arrayK, f, cn);
  return tpl;
}
