/*    
    Copyright 2013-2017 Onera.

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
# include <stdio.h>
#include <math.h>
# include "connector.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Calcule et stocke les coefficients d'interpolation 
   CAS SANS DOUBLE WALL 
   IN: receiverArray: points a interpoler définis sous forme de maillage
   IN: donorArrays: maillages donneurs. La localisation du donneur 
        (noeuds/centres/centres étendus) doit être effectuee au prealable
   IN: Order: ordre des interpolations (2, 3, 5)
   IN: Nature: 0: produit des cellN=0 -> donneur invalide; 
               1: cellN=0 ou 2 -> donneur invalide
   IN: PenaliseBorders: 1: penalite sur le volume des pts ou cellules frontieres
   IN: allHooks != Py_None: un hook par adt associé à un donor 
   OUT: [donorBlks,donorInd1D, donorType, coefs, extrapInd1D, orphanInd1D] 
        donorBlks: no du blk donneur, démarre à 0
        donorInd1D: indice global (structure), de l elt (NS) du donneur
        donorType: type d interpolation effectué localement
        coefs: coefficients d interpolation, stockés selon le type
        extrapInd1D: indices des pts extrapolés
        orphanInd1D: indices des pts orphelins
*/
//=============================================================================
PyObject* K_CONNECTOR::indiceToCoord2(PyObject* self, PyObject* args)
{
  PyObject *indiceslist; // domaines d interpolation
  PyObject *rangedonor;
  PyObject *transfo;
  PyObject *profondeur;
  PyObject *dirD;
  PyObject *typ;
  E_Int dir;
  E_Int nb_ind;
  E_Int ni;
  E_Int nj;
  E_Int nk;

  if (!PYPARSETUPLEI(args,
                    "OOOOOOlllll", "OOOOOOiiiii",
		     &indiceslist, &rangedonor,&transfo, &profondeur, &dirD, &typ, &dir, &nb_ind, &ni, &nj, &nk))
  {
      return NULL;
  }

  ni=ni-1;
  nj=nj-1;
  nk=nk-1;

  cout << ni <<" "<<nj<<" "<<nk << endl;

  E_Int i;
  E_Int j;
  E_Int k;
  E_Int imin=ni;
  E_Int imax=1;
  E_Int jmin=nj;
  E_Int jmax=1;
  E_Int kmin=nk;
  E_Int kmax=1;
  E_Int li;
  E_Int lj;
  E_Int lk;
  // E_Int dirD;
 
 /// Recuperation du tableau de stockage des valeurs
  //PyObject* indiceslist = PyList_GetItem(indiceslist,0); FldArrayF* drodm;
  //K_NUMPY::getFromNumpyArray(drodmArray, drodm, true); E_Float* iptdrodm = drodm->begin();

  FldArrayI* ind_list;
  K_NUMPY::getFromNumpyArray(indiceslist,ind_list, true); E_Int* ipt_ind_list = ind_list->begin();

  FldArrayI* rangedonor_;
  K_NUMPY::getFromNumpyArray(rangedonor,rangedonor_, true);  E_Int* ipt_rangedonor = rangedonor_->begin();

  FldArrayI* transfo_;
  K_NUMPY::getFromNumpyArray(transfo,transfo_, true);  E_Int* ipt_transfo = transfo_->begin();

  FldArrayI* profondeur_;
  K_NUMPY::getFromNumpyArray(profondeur,profondeur_, true);  E_Int* ipt_profondeur = profondeur_->begin();

  FldArrayI* dirD_;
  K_NUMPY::getFromNumpyArray(dirD,dirD_, true);  E_Int* ipt_dirD = dirD_->begin();

  FldArrayI* typ_;
  K_NUMPY::getFromNumpyArray(typ,typ_, true);  E_Int* ipt_typ = typ_->begin();

  //  for (E_Int i=0; i < 3; i++)
  //  {
  //   if ( abs(ipt_transfo[i]) == abs(dir))
  //	{
	  
  //	  ipt_dirD[0] = (i+1) * (-dir/abs(dir))*(ipt_transfo[i]/abs(ipt_transfo[i]));

  //	}

  // }



   
   for (E_Int ind = 0 ; ind < nb_ind ; ind++)
    {

      //cout << ipt_typ[ind]  << endl;
      if (ipt_typ[ind]==22)

	{

	  k = floor(ipt_ind_list[ind]/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind] - (k-1)*ni*nj;
	  j = floor(j/ni) + 1;
	  i = ipt_ind_list[ind] - (k-1)*ni*nj - (j-1)*ni + 1;

	  if (i < imin){ imin=i;}
	  if (i > imax){ imax=i;}
	  if (j < jmin){ jmin=j;}
	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+1)/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind]+1 - (k-1)*ni*nj;
	  j = floor(j/ni) + 1;
	  i = ipt_ind_list[ind]+1 - (k-1)*ni*nj - (j-1)*ni + 1;

	  if (i < imin){ imin=i;}
	  if (i > imax){ imax=i;}
	  if (j < jmin){ jmin=j;}
	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+ni)/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind]+ni - (k-1)*ni*nj;
	  j = floor(j/ni) + 1;
	  i = ipt_ind_list[ind]+ni - (k-1)*ni*nj - (j-1)*ni + 1;

	  if (i < imin){ imin=i;}
	  if (i > imax){ imax=i;}
	  if (j < jmin){ jmin=j;}
	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+ni+1)/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind]+ni+1 - (k-1)*ni*nj;
	  j = floor(j/ni) + 1;
	  i = ipt_ind_list[ind]+ni+1 - (k-1)*ni*nj - (j-1)*ni + 1;

	  if (i < imin){ imin=i;}
	  if (i > imax){ imax=i;}
	  if (j < jmin){ jmin=j;}
	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}


	}


      else if (ipt_typ[ind]==1)

	{

	  k = floor(ipt_ind_list[ind]/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind] - (k-1)*ni*nj;
 	  j = floor(j/ni) + 1;
	  i = ipt_ind_list[ind] - (k-1)*ni*nj - (j-1)*ni + 1;

	  if (i < imin){ imin=i;}
	  if (i > imax){ imax=i;}
	  if (j < jmin){ jmin=j;}
	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	}

      else if (ipt_typ[ind]==2)

	{  


	  k = floor(ipt_ind_list[ind]/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind] - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i = ipt_ind_list[ind] - (k-1)*ni*nj - (j-1)*ni + 1;

   	  if (i < imin){ imin=i;}
   	  if (i > imax){ imax=i;}
   	  if (j < jmin){ jmin=j;}
   	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+1)/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind]+1 - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i =  ipt_ind_list[ind]+1 - (k-1)*ni*nj - (j-1)*ni + 1;

   	  if (i < imin){ imin=i;}
   	  if (i > imax){ imax=i;}
   	  if (j < jmin){ jmin=j;}
   	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+ni)/(ni*nj)) + 1;
   	  j = ipt_ind_list[ind]+ni - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i = ipt_ind_list[ind]+ni - (k-1)*ni*nj -  (j-1)*ni + 1;

   	  if (i < imin){ imin=i;}
   	  if (i > imax){ imax=i;}
   	  if (j < jmin){ jmin=j;}
   	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+ni+1)/(ni*nj)) + 1;
   	  j = ipt_ind_list[ind]+ni+1 - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i = ipt_ind_list[ind]+ni+1 - (k-1)*ni*nj - (j-1)*ni + 1;

   	  if (i < imin){ imin=i;}
   	  if (i > imax){ imax=i;}
   	  if (j < jmin){ jmin=j;}
   	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+ni*nj)/(ni*nj)) + 1;
   	  j = ipt_ind_list[ind]+ni*nj - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i = ipt_ind_list[ind]+ni*nj - (k-1)*ni*nj - (j-1)*ni + 1;

   	  if (i < imin){ imin=i;}
   	  if (i > imax){ imax=i;}
   	  if (j < jmin){ jmin=j;}
   	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+ni*nj+1)/(ni*nj)) + 1;
   	  j = ipt_ind_list[ind]+ni*nj+1 - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i = ipt_ind_list[ind]+ni*nj+1 - (k-1)*ni*nj -  (j-1)*ni + 1;

   	  if (i < imin){ imin=i;}
   	  if (i > imax){ imax=i;}
   	  if (j < jmin){ jmin=j;}
   	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+ni*nj+ni)/(ni*nj)) + 1;
   	  j = ipt_ind_list[ind]+ni*nj+ni - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i = ipt_ind_list[ind]+ni*nj+ni - (k-1)*ni*nj - (j-1)*ni + 1;

   	  if (i < imin){ imin=i;}
   	  if (i > imax){ imax=i;}
   	  if (j < jmin){ jmin=j;}
   	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}

	  k = floor((ipt_ind_list[ind]+ni*nj+ni+1)/(ni*nj)) + 1;
   	  j = ipt_ind_list[ind]+ni*nj+ni+1 - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i = ipt_ind_list[ind]+ni*nj+ni+1 - (k-1)*ni*nj - (j-1)*ni + 1;

   	  if (i < imin){ imin=i;}
   	  if (i > imax){ imax=i;}
   	  if (j < jmin){ jmin=j;}
   	  if (j > jmax){ jmax=j;}
   	  if (k < kmin){ kmin=k;}
   	  if (k > kmax){ kmax=k;}
 
	}

    }


   cout << imin <<" "<< imax <<" "<<jmin<<" "<<jmax<<" "<<kmin<<" "<<kmax<<" "<<  endl;


   li = imax - imin;
   lj = jmax - jmin;
   lk = kmax - kmin;


   // Correction des indices
   kmin=kmin-2;
   if (kmin == -1 or kmin== 0)
     {
       kmin = 1;
     }
   kmax = kmin + lk;
   if (kmax > nk-4)
     {
       kmax = nk -4;
     } 
   if (kmax <= 0)
     {
       kmax = 1;
     } 


   jmin=jmin-2;
   if (jmin == -1 or jmin== 0)
     {
       jmin = 1;
     }
   jmax = jmin + lj;
   if (jmax > nj - 4)
     {
       jmax = nj -4;
     } 
   if (jmax <= 0)
     {
       jmax = 1;
     } 

   imin=imin-2;
   if (imin == -1 or imin== 0)
     {
       imin = 1;
     }
   imax = imin + li;
   if (imax > ni-4)
     {
       imax = ni -4;
     } 
   if (imax <= 0)
     {
       imax = 1;
     } 

   cout << imin <<" "<< imax <<" "<<jmin<<" "<<jmax<<" "<<kmin<<" "<<kmax<<" "<<  endl;

   if (nk == 1) // 2D

     {

       if (lj >= li and imax > (ni-1)/2)  // dirR = +1
	{
	  ipt_dirD[0] = 1;
	}
      else if (lj >= li and imin < (ni-1)/2)  // dirR = -1
	{
	  ipt_dirD[0] = -1;
	}
      else if (li >= lj and jmax > (nj-1)/2)  // dirR = +2
	{
	  ipt_dirD[0] = 2;
	}
      else if (li >= lj and jmin < (nj-1)/2)  // dirR = -2
	{
	  ipt_dirD[0] = -2;
	}

     }

   else //3D

     {

      if (lk == min(min(li,lj),lk) and kmax > (nk-1)/2)  // dirR = -3
	{
	  ipt_dirD[0] = -3;
	}
      else if (lk == min(min(li,lj),lk) and kmin < (nk-1)/2)  // dirR = 3
	{
	  ipt_dirD[0] = 3;
	}
      else if (lj == min(min(li,lj),lk) and jmax > (nj-1)/2)  // dirR = +2
	{
	  ipt_dirD[0] = 2;
	}
      else if (lj == min(min(li,lj),lk) and jmin < (nj-1)/2)  // dirR = -2
	{
	  ipt_dirD[0] = -2;
	}
      else if (li == min(min(li,lj),lk) and imax > (ni-1)/2)  // dirR = +2
	{
	  ipt_dirD[0] = 1;
	}
      else if (li == min(min(li,lj),lk) and imin < (ni-1)/2)  // dirR = -2
	{
	  ipt_dirD[0] = -1;
	}



     }
      
   

 if (ipt_dirD[0]==1)
   {
     ipt_profondeur[0] = imax - imin ;
     imin = imax;
   }
 else if (ipt_dirD[0]==-1)
   {
     ipt_profondeur[0] = imax - imin ;
     imax = imin;
   }
 else if (ipt_dirD[0]==2)
   {
     ipt_profondeur[0] = jmax - jmin ;
     jmin = jmax;
   }
 else if (ipt_dirD[0]==-2)
   {
     ipt_profondeur[0] = jmax - jmin ;
     jmax = jmin;
   }
 else if (ipt_dirD[0]==3)
   {
     ipt_profondeur[0] = kmax - kmin ;
     kmin = kmax;
   }
 else if (ipt_dirD[0]==-3)
   {
     ipt_profondeur[0] = kmax - kmin ;
     kmax = kmin;
   }

 //cout << "profondeur= " << ipt_profondeur[0] << endl;
 /*
   kmax=kmax-2;
   if (kmax == nk-2 or kmax == nk-3)
     {
       kmax = nk -4;
     }
   if (kmax <= 0)
     {
       kmax = 1;
     }

   kmin=kmin-2;
   if (kmin == -1 or kmin== 0)
     {
       kmin = 1;
     }

   jmax=jmax-2;
   if (jmax == nj-2 or jmax== nj-3)
     {
       jmax = nj -4;
     }
   if (jmax <= 0)
     {
       jmax = 1;
     }


   jmin=jmin-2;
   if (jmin == -1 or jmin== 0)
     {
       jmin = 1;
     }

   imax=imax-2;
   if (imax == ni-2 or imax== ni-3)
     {
       imax = ni -4;
     }
   if (imax <= 0)
     {
       imax = 1;
     }

   imin=imin-2;
   if (imin == -1 or imin== 0)
     {
       imin = 1;
     }

 */
      
  ipt_rangedonor[0]=imin;
  ipt_rangedonor[1]=imax;
  ipt_rangedonor[2]=jmin;
  ipt_rangedonor[3]=jmax;
  ipt_rangedonor[4]=kmin;
  ipt_rangedonor[5]=kmax;



 RELEASESHAREDN( indiceslist  , ind_list );
 RELEASESHAREDN( rangedonor  , rangedonor_ );
 RELEASESHAREDN( transfo  , transfo_ );
 RELEASESHAREDN( profondeur  , profondeur_ );
 RELEASESHAREDN( dirD  , dirD_ );

 Py_INCREF(Py_None);
 return Py_None;

}
