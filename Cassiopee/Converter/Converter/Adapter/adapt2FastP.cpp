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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Adaptation connectivite pour FastP
   1. Tri des faces internes et externes
   2. ajoute un element plat a la fin dans NFACE
*/
//=============================================================================
PyObject* K_CONVERTER::adapt2FastP(PyObject* self, PyObject* args)
{
  PyObject* NGon; PyObject* NFace; PyObject* PE; PyObject* NGon_intext; PyObject* NFace_intext; PyObject* Ptlist_bc; PyObject* Ptlist_rac; PyObject* Ptlist_racD;
  E_Int nelts;
  if (!PYPARSETUPLE_(args, OOOO_ OOOO_ I_, &NGon, &NFace, &PE, &NGon_intext, &NFace_intext, &Ptlist_bc, &Ptlist_rac, &Ptlist_racD, &nelts))
    return NULL;

  // Check numpy NGon
  FldArrayI* cNGon;
  E_Int res = K_NUMPY::getFromNumpyArray(NGon, cNGon);
  if (res == 0) { PyErr_SetString(PyExc_TypeError, "adapt2FastP: NGon is invalid."); return NULL; }

  // Check numpy NFace
  FldArrayI* cNFace;
  res = K_NUMPY::getFromNumpyArray(NFace, cNFace);
  if (res == 0) { PyErr_SetString(PyExc_TypeError, "adapt2FastP: NFace is invalid."); return NULL; }

  // Check numpy PE (ParentElements)
  FldArrayI* cPE;
  res = K_NUMPY::getFromNumpyArray(PE, cPE);
  if (res == 0) { PyErr_SetString(PyExc_TypeError, "adapt2FastP: numpy is invalid."); return NULL; }

  // Check numpy  Nb face interne/extern (Ngon)
  FldArrayI* cNGon_intext;
  res = K_NUMPY::getFromNumpyArray(NGon_intext, cNGon_intext);
  if (res == 0) { PyErr_SetString(PyExc_TypeError, "adapt2FastP:  NGon IntExt is invalid."); return NULL; }

  // Check numpy  Nb elts interne/extern (Nface)
  FldArrayI* cNFace_intext;
  res = K_NUMPY::getFromNumpyArray(NFace_intext, cNFace_intext);
  if (res == 0) { PyErr_SetString(PyExc_TypeError, "adapt2FastP: NFace IntExt is invalid."); return NULL; }


  E_Int NbPtlist_bc  = PyList_Size(Ptlist_bc  );
  E_Int NbPtlist_rac = PyList_Size(Ptlist_rac );
  E_Int NbPtlist_racD= PyList_Size(Ptlist_racD);

  E_Int* size_ptlist_bc = new E_Int [  NbPtlist_bc  *3];
  E_Int* size_ptlist_rac = new E_Int [ NbPtlist_rac *3];
  E_Int* size_ptlist_racD = new E_Int [NbPtlist_racD*3];

  E_Int** ipt_ptlist_bc; E_Int** ipt_ptlist_rac; E_Int** ipt_ptlist_racD;
  ipt_ptlist_bc  = new E_Int*[NbPtlist_bc + NbPtlist_rac + NbPtlist_racD ];
  ipt_ptlist_rac = ipt_ptlist_bc  + NbPtlist_bc;
  ipt_ptlist_racD= ipt_ptlist_rac + NbPtlist_rac;

  for (E_Int i = 0; i < NbPtlist_bc; i++)
  {
    PyObject* ptlistArray  = PyList_GetItem(Ptlist_bc,i); 
    FldArrayI* Ptlist;
    K_NUMPY::getFromPointList(ptlistArray, Ptlist);
    ipt_ptlist_bc[i]   = Ptlist->begin();
    size_ptlist_bc[i]  = Ptlist->getSize();
    //printf(" fen bc= %d %d %d \n", size_ptlist_bc[i], Ptlist->getSize(),  i);
  }
  for (E_Int i = 0; i < NbPtlist_rac; i++)
  {
    PyObject* ptlistArray  = PyList_GetItem(Ptlist_rac,i); 
    FldArrayI* Ptlist;
    K_NUMPY::getFromPointList(ptlistArray, Ptlist); 
    ipt_ptlist_rac[i] = Ptlist->begin();

    size_ptlist_rac[i] = Ptlist->getSize();
     //printf(" fen rac= %d %d \n",  size_ptlist_rac[i], i);
  }
  for (E_Int i = 0; i < NbPtlist_rac; i++)
  {
    PyObject* ptlistArray  = PyList_GetItem(Ptlist_racD,i); 
    FldArrayI* Ptlist;
    K_NUMPY::getFromPointList(ptlistArray, Ptlist); 
    ipt_ptlist_racD[i] = Ptlist->begin();

    size_ptlist_racD[i]= Ptlist->getSize();
    // printf(" fen racD= %d %d \n", size_ptlist_racD[i], i);
  }

  // Compte les faces int et ext
  E_Int nfaces = cPE->getSize();

  E_Int* ptNGon = cNGon->begin();
  FldArrayI posFaces(nfaces);
  K_CONNECT::getPosFacets(ptNGon, 0, nfaces, posFaces);
  E_Int* posFacesp = posFaces.begin();
  E_Int* ptr;

  E_Int* ng_intext = cNGon_intext->begin(1);
  E_Int* nf_intext = cNFace_intext->begin(1);

  E_Int nfaces_int0 = ng_intext[0]; //nb face interne couche zero
  E_Int nfaces_rac0 = ng_intext[1]; //nb face raccord couche zero
  E_Int nfaces_bc0  = ng_intext[2]; //nb face bc      couche zero
  E_Int nfaces_int1 = ng_intext[3]; //nb face interne couche un
  E_Int nfaces_bc1  = ng_intext[4]; //nb face bc      couche un   
  E_Int nfaces_int2 = ng_intext[5]; //nb face         couche deux

  E_Int Nelt0  = nf_intext[0]; //nb element  couche zero

  //printf(" NFACETOT=  %d \n", nfaces); 
  //on determine nombre face couche 0 et 1 pour chaque BC
   for (E_Int i = 0; i < NbPtlist_bc; i++)
  {
     E_Int face = size_ptlist_bc[i] - 1;
     E_Int c    = 0;
     while ( ipt_ptlist_bc[i][face] > nfaces_int0 + nfaces_rac0 + nfaces_bc0 )
       { c   = c+1;
        face = face-1;
       }
    size_ptlist_bc[i+NbPtlist_bc  ] =  size_ptlist_bc[i] - c;  // nombre face couche 0
    size_ptlist_bc[i+NbPtlist_bc*2] =  c;                     //nombre de face couche 1
     //printf(" size bc= %d %d %d %d  \n", i, size_ptlist_bc[i], size_ptlist_bc[i+NbPtlist_bc],size_ptlist_bc[i+NbPtlist_bc*2]);
  } 

  //printf("Nfi0= %d, Nfrac0= %d, Nfbc0= %d, Nfi1=%d, Nfbc1= %d, Nfi2=%d  \n", nfaces_int0,nfaces_rac0,nfaces_bc0,nfaces_int1,nfaces_bc1,nfaces_int2);

  E_Int* PEG = cPE->begin(1);
  E_Int* PED = cPE->begin(2);

  // new interface position (face externe mise a la fin):  nn
  E_Int* nn  = new E_Int [nfaces];
  E_Int* nni = new E_Int [nfaces];
  
  E_Int ni0 = 0;
  E_Int nrac= nfaces_int0;
  E_Int ni1 = nrac  + nfaces_rac0;
  E_Int nbc1= ni1   + nfaces_int1;
  E_Int nbc0= nbc1  + nfaces_bc1;
  //E_Int ni2 = nbc0  + nfaces_bc0;

  //printf(" pt rac= %d, ptfi1= %d, ptbc0= %d, ptbc1= %d \n", nrac, ni1, nbc0,nbc1);

  //tri face couche zero
  E_Int  cbc, crac, ci;
  cbc =0; crac=0; ci=0;
  E_Int nfaces0 = nfaces_int0 + nfaces_rac0 + nfaces_bc0;

  //tri face intern couche zero
  //printf(" NFACE0=  %d \n", nfaces0); 
  for (E_Int i = 0; i < nfaces0; i++)
  {
    if ( (PED[i] <=  Nelt0) && (PEG[i] <= Nelt0) && (  PED[i]*PEG[i] != 0)) { nn[i] = ni0 ; ni0++;  ci++;} // printf(" Fint0: %d %d %d %d \n",nn[i],i, PED[i],PEG[i] ); } //face intern
    //printf(" coue0: %d %d %d %d \n",i,nn[i],PED[i],PEG[i] );
  }
  
  //tri face rac couche zero
  for (E_Int rac = 0; rac < NbPtlist_rac; rac++)
  {
     E_Int nface_rac = size_ptlist_rac[rac];
     for (E_Int l = 0; l < nface_rac; l++)
     { E_Int i =  ipt_ptlist_rac[ rac ][l]-1;
       nn[i] = nrac; nrac++; crac++;
       //printf("rac nni %d %d %d  \n", nn[i], i, rac); 
     }
  }
  //tri face bcc couche zero
  for (E_Int bc = 0; bc < NbPtlist_bc; bc++)
  {
     E_Int nface_bc = size_ptlist_bc[bc + NbPtlist_bc];
     for (E_Int l = 0; l < nface_bc; l++)
     { E_Int i =  ipt_ptlist_bc[ bc ][l]-1;
       nn[i] = nbc0; nbc0++; cbc++;
       //printf("bc0 nni %d %d %d  \n", nn[i], i, bc); 
     }
  }

  //tri face interne  couche un   
  cbc =0; crac=0; ci=0;
  E_Int nfaces1 = nfaces_int1 + nfaces_bc1;
  for (E_Int i = nfaces0 ; i < nfaces0 + nfaces1  ; i++)
  {
    if  (  PED[i]*PEG[i] != 0)    { nn[i] = ni1 ; ni1++;  ci++; }   //printf(" couche1: %d %d %d %d \n",nn[i],i,PED[i],PEG[i] );    } //face intern
  }
  //tri face bcc couche un
  for (E_Int bc = 0; bc < NbPtlist_bc; bc++)
  {
     for (E_Int l = size_ptlist_bc[bc + NbPtlist_bc] ; l < size_ptlist_bc[bc]; l++)
     { E_Int i =  ipt_ptlist_bc[ bc ][l]-1;
       nn[i] = nbc1; nbc1++; cbc++;
       //printf("bc1 nni %d %d %d  \n", nn[i], i, bc); 
     }
  }
  //printf("couche1 %d %d %d %d \n", cbc, crac, ci, nfaces1);
  //face couche deux: rien a trier   
  E_Int nfaces2 = nfaces_int2;
  for (E_Int i = nfaces0 + nfaces1; i < nfaces0 + nfaces1  + nfaces2 ; i++)  nn[i] = i;

  //for (E_Int i = 0; i < nfaces; i++) printf("i= %d,  nn=%d, nni= %d, PED=, %d, PEG= %d \n",i,nn[i],nni[i], PED[i] ,PEG[i]);
  E_Int e1, e2, nf, ints;

  for (E_Int i = 0; i < nfaces; i++) nni[ nn[i] ] = i;



  // modify PE
  FldArrayI temp(nfaces,2);
  E_Int* temp1 = temp.begin(1); E_Int* temp2 = temp.begin(2);
  for (E_Int i = 0; i < nfaces; i++)
  {
    e1 = PED[i]; e2 = PEG[i];
    temp1[nn[i]] = e1;
    temp2[nn[i]] = e2;
  }
  for (E_Int i = 0; i < nfaces; i++)
  {
    PED[i] = temp1[i]; PEG[i] = temp2[i];
  //printf("PEDG0 %d %d %d \n", PED[i] ,PEG[i] ,i);
  }
  temp.malloc(0);

  // modify ngon elementconnectivity
  temp.malloc(cNGon->getSize());
  E_Int* tt = temp.begin();
  for (E_Int i = 0; i < nfaces; i++)
  {
    ptr = ptNGon+posFacesp[nni[i]];
    nf = ptr[0];
    for (E_Int j = 0; j <= nf; j++) tt[j] = ptr[j];
    tt += nf+1;
  }

  tt = temp.begin();
  for (E_Int i = 0; i < cNGon->getSize(); i++)
  {
    ptNGon[i] = tt[i];
  }
  posFaces.malloc(0);

  // modify BC
  for (E_Int i = 0; i < NbPtlist_bc; i++)
  {
    for (E_Int j = 0; j < size_ptlist_bc[i]; j++)
    {
      E_Int face = ipt_ptlist_bc[i][j];
      ipt_ptlist_bc[i][j] = nn[face-1]+1;
      //printf(" pt rac= %d %d %d %d\n", nn[face-1]+1, face, j,i);
    }
  }
  // modify pointlist Match
  for (E_Int i = 0; i < NbPtlist_rac; i++)
  {
    for (E_Int j = 0; j < size_ptlist_rac[i]; j++)
    {
      E_Int face = ipt_ptlist_rac[i][j];
      ipt_ptlist_rac[i][j] = nn[face-1]+1;
    }
  }
  // modify pointlistdonor Match
  for (E_Int i = 0; i < NbPtlist_racD; i++)
  {
    for (E_Int j = 0; j < size_ptlist_racD[i]; j++)
    {
      E_Int face = ipt_ptlist_racD[i][j];
      ipt_ptlist_racD[i][j] = nn[face-1]+1;
    }
  }

  // modify NFACES (returned) + 
  // un element ajoute par face externe (ghost cells)
  // pas encore de tri des elements
  E_Int* ptNFace = cNFace->begin();
  E_Int size     = cNFace->getSize();
  FldArrayI posElts(nelts);
  K_CONNECT::getPosFacets(ptNFace, 0, nelts, posElts);
  E_Int* posEltsp = posElts.begin();

  //ajout de 1 element contenant 1 seule face pour chaque face BC
  //besoin de 2 stockage: nb de face=1  et No de face
  E_Int new_elts= nfaces_bc0+nfaces_bc1;
  PyObject* nnfo = K_NUMPY::buildNumpyArray(size+2*new_elts, 1, 1, 1); // + exts
  E_Int*     nnf = K_NUMPY::getNumpyPtrI(nnfo);                       // nouveau NFace

  ptr = cNFace->begin();
  for (E_Int i = 0; i < nelts; i++)
  {
    ptr = ptNFace+posEltsp[i];
    nf = ptr[0];
    //printf("%d %d\n",nf,nn[i]);
    nnf[0] = nf;
    for (E_Int j = 1; j <= nf; j++) { nnf[j] = nn[ptr[j]-1]+1; }
    nnf += nf+1;
  }
  
  // Ajout des ghost cells sur les faces BC
  ints = nfaces_int0 + nfaces_rac0 + nfaces_int1;

  for (E_Int i = 0; i < new_elts; i++)
  {
    nnf[0] = 1;        // 1 face
    nnf[1] = ints+i+1; // no de face externe
    nnf += 2;
  }
  
  // modify PE pour elements ghosts
  for (E_Int i = 0; i < new_elts; i++)
  {
    if (PEG[i+ints] == 0)
    {
      printf(" PEGGGGGGG nul " SF_D2_ "\n", ints, i);
      PEG[i+ints] = nelts+i+1;
    }

    if (PED[i+ints] == 0) 
    {
      //printf(" PED nul %d %d %d \n", PED[i+ints], ints, i);
      //PED[i+ints] = PEG[i+ints];
      //PEG[i+ints] = nelts+i+1;

      // Modif ivan
      PED[i+ints] = nelts+i+1;
    }
  }

//  for (E_Int i = 0; i < nfaces; i++)
//  {
//  printf("PEDG2 %d %d %d \n", PED[i] ,PEG[i] ,i);
//  }

  delete [] nn; delete [] nni; delete [] ipt_ptlist_bc;
  delete [] size_ptlist_bc;
  delete [] size_ptlist_rac;
  delete [] size_ptlist_racD;
  RELEASESHAREDN(NGon, cNGon);
  RELEASESHAREDN(PE, cPE);
  RELEASESHAREDN(NFace, cNFace);
  RELEASESHAREDN(NFace_intext, cNFace_intext);
  RELEASESHAREDN(NGon_intext, cNGon_intext);
  return nnfo;
}
