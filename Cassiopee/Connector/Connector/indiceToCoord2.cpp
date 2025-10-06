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
  E_Int ni, nj, nk;

  if (!PYPARSETUPLE_(args, OOOO_ OO_ IIII_ I_,
		     &indiceslist, &rangedonor,&transfo, &profondeur, &dirD, &typ, &dir, &nb_ind, &ni, &nj, &nk))
  {
      return NULL;
  }

  //cout <<"ni= "<< ni <<" "<<"nj= "<<nj<<" "<<"nk= " << nk << endl;

  E_Int i, j, k;
  E_Int imin=ni;
  E_Int imax=1;
  E_Int jmin=nj;
  E_Int jmax=1;
  E_Int kmin=nk;
  E_Int kmax=1;
  E_Int li;
  E_Int lj;
  E_Int lk;
 
 /// Recuperation du tableau de stockage des valeurs
  FldArrayI* ind_list;
  K_NUMPY::getFromNumpyArray(indiceslist, ind_list); E_Int* ipt_ind_list = ind_list->begin();

  FldArrayI* rangedonor_;
  K_NUMPY::getFromNumpyArray(rangedonor, rangedonor_); E_Int* ipt_rangedonor = rangedonor_->begin();

  FldArrayI* transfo_;
  K_NUMPY::getFromNumpyArray(transfo, transfo_); //E_Int* ipt_transfo = transfo_->begin();

  FldArrayI* profondeur_;
  K_NUMPY::getFromNumpyArray(profondeur, profondeur_); E_Int* ipt_profondeur = profondeur_->begin();

  FldArrayI* dirD_;
  K_NUMPY::getFromNumpyArray(dirD, dirD_); E_Int* ipt_dirD = dirD_->begin();

  FldArrayI* typ_;
  K_NUMPY::getFromNumpyArray(typ, typ_); E_Int* ipt_typ = typ_->begin();

   //
   //Determination fenetre donneur
   //
   for (E_Int ind = 0 ; ind < nb_ind ; ind++)
    {
      if (ipt_typ[ind]==22) // interp ordre 2 2D
	{
	  k = floor(ipt_ind_list[ind]/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind] - (k-1)*ni*nj;
	  j = floor(j/ni) + 1;
	  i = ipt_ind_list[ind] - (k-1)*ni*nj - (j-1)*ni + 1;

	  if (i   < imin){ imin=i;  }
	  if (i+1 > imax){ imax=i+1;}
	  if (j   < jmin){ jmin=j;  }
	  if (j+1 > jmax){ jmax=j+1;}
   	  if (k   < kmin){ kmin=k;  }
   	  if (k   > kmax){ kmax=k;  }
	}
      else if (ipt_typ[ind]==1) //abutting
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
      else if (ipt_typ[ind]==2) // interp ordre 2 3D
	{  
	  k = floor(ipt_ind_list[ind]/(ni*nj)) + 1;
   	  j =  ipt_ind_list[ind] - (k-1)*ni*nj;
	  j = floor(j/ni)+1;
    	  i = ipt_ind_list[ind] - (k-1)*ni*nj - (j-1)*ni + 1;

   	  if (i   < imin){ imin=i;  }
   	  if (i+1 > imax){ imax=i+1;}
   	  if (j   < jmin){ jmin=j;  }
   	  if (j+1 > jmax){ jmax=j+1;}
   	  if (k   < kmin){ kmin=k;  }
   	  if (k+1 > kmax){ kmax=k+1;}
	}
    }

   if (imin < 3 or jmin < 3){cout << "danger_min" << endl;}
   if (jmax > nj-2 or imax > ni-2){cout << "danger_max" << endl;}

   //epaisseur couche donneuse
   li = imax - imin;
   lj = jmax - jmin;
   lk = kmax - kmin;

   // Correction des indices: c --> fortran Fast??
   kmin=kmin-2;
   if (kmin == -1 or kmin== 0) { kmin = 1; }
   
   ///proposition
   kmax = kmax - 2;
   if (kmax > nk-4) { kmax = nk-4;} 
   if (kmax <= 0  ) { kmax = 1; } 

   jmin=jmin-2;
   if (jmin == -1 or jmin== 0) { jmin = 1; }
   ///proposition
   jmax = jmax - 2;
   if (jmax > nj - 4) { jmax = nj -4;} 
   if (jmax <= 0    ) { jmax = 1; } 

   imin=imin-2;
   if (imin == -1 or imin== 0) { imin = 1; }

   ///proposition
   imax = imax - 2;
   if (imax > ni-4) { imax = ni -4;} 
   if (imax <= 0  ) { imax = 1; } 

   cout <<"indices corriges= "<< imin <<" "<< imax <<" "<<jmin<<" "<<jmax<<" "<<kmin<<" "<<kmax<<  endl;

   //on determine si fenetre donneuse est en imin , imax,.....
   if (nk == 1) // 2D
     {
      if      (lj >= li and imax > (ni-1)/2)
        { ipt_dirD[0] = 1;
          ipt_profondeur[0] = imax - imin ;
          imin = imax;
        } 
      else if (lj >= li and imin < (ni-1)/2)
        { ipt_dirD[0] =-1;
          ipt_profondeur[0] = imax - imin ;
          imax = imin;
        } 
      else if (li >= lj and jmax > (nj-1)/2)
        { ipt_dirD[0] = 2;
          ipt_profondeur[0] = jmax - jmin ;
          jmin = jmax;
        } 
      else if (li >= lj and jmin < (nj-1)/2)
        { ipt_dirD[0] =-2;
          ipt_profondeur[0] = jmax - jmin ;
          jmax = jmin;
        } 
     }
   else //3D
     {
      if      (lk == min(min(li,lj),lk) and kmax > (nk-1)/2)
        { ipt_dirD[0] = 3;
          ipt_profondeur[0] = kmax - kmin ;
          kmin = kmax;
        } 
      else if (lk == min(min(li,lj),lk) and kmin < (nk-1)/2)
        { ipt_dirD[0] =-3;
          ipt_profondeur[0] = kmax - kmin ;
          kmax = kmin;
        } 
      else if (lj == min(min(li,lj),lk) and jmax > (nj-1)/2)
        { ipt_dirD[0] = 2;
          ipt_profondeur[0] = jmax - jmin ;
          jmin = jmax;
        } 
      else if (lj == min(min(li,lj),lk) and jmin < (nj-1)/2)
        { ipt_dirD[0] =-2;
          ipt_profondeur[0] = jmax - jmin ;
          jmax = jmin;
        } 
      else if (li == min(min(li,lj),lk) and imax > (ni-1)/2)
        { ipt_dirD[0] = 1;
          ipt_profondeur[0] = imax - imin ;
          imin = imax;
        } 
      else if (li == min(min(li,lj),lk) and imin < (ni-1)/2)
        { ipt_dirD[0] =-1;
          ipt_profondeur[0] = imax - imin ;
          imax = imin;
        } 
     }
      
  cout <<"dirD= "<< ipt_dirD[0] <<  endl;
  cout << "profondeur= " << ipt_profondeur[0] << endl;
      
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
 RELEASESHAREDN( typ, typ_ );

 Py_INCREF(Py_None);


 return Py_None;


}
