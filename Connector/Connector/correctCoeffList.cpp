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
PyObject* K_CONNECTOR::correctCoeffList(PyObject* self, PyObject* args)
{
  PyObject *indiceslist; // domaines d interpolation
  PyObject *coefflist;
  PyObject *typ;
  E_Int dim;
  E_Int nb_ind;
  E_Int ni;
  E_Int nj;
  E_Int nk;

  if (!PYPARSETUPLEI(args,
                    "OOOllll", "OOOiiii",
		     &indiceslist, &coefflist, &typ, &nb_ind,  &ni, &nj, &nk))
  {
      return NULL;
  }

  ni=ni-1;
  nj=nj-1;
  nk=nk-1;
  E_Int i;
  E_Int j;
  E_Int k;
  E_Int compt=0;
  E_Int imin=ni;
  E_Int imax=1;
  E_Int jmin=nj;
  E_Int jmax=1;
  E_Int kmin=nk;
  E_Int kmax=1;
  E_Int indCoef=0;
  E_Float a;
  E_Float b;
  E_Float xc;
  E_Float yc;
  E_Float zc;
  E_Float S=0.0;
  E_Float tab[8][3];
  E_Int tab_[8];
  E_Int dirD;
 
 /// Recuperation du tableau de stockage des valeurs
  //PyObject* indiceslist = PyList_GetItem(indiceslist,0); FldArrayF* drodm;
  //K_NUMPY::getFromNumpyArray(drodmArray, drodm, true); E_Float* iptdrodm = drodm->begin();

  FldArrayI* ind_list;
  K_NUMPY::getFromNumpyArray(indiceslist,ind_list, true); E_Int* ipt_ind_list = ind_list->begin();

  FldArrayF* coefflist_;
  K_NUMPY::getFromNumpyArray(coefflist,coefflist_, true);  E_Float* ipt_coefflist = coefflist_->begin();

  FldArrayI* typ_;
  K_NUMPY::getFromNumpyArray(typ,typ_, true);  E_Int* ipt_typ = typ_->begin();

  //cout << "nb_ind= " << nb_ind << endl;

 for (E_Int ind = 0 ; ind < nb_ind ; ind++)

   {

     //cout << "type= " << int(ipt_typ[ind]) << endl;

   if (ipt_typ[ind]==22)

     {
      
      	  k = floor(ipt_ind_list[ind]/(ni*nj)) + 1;
	  j = ipt_ind_list[ind] - (k-1)*ni*nj;
	  j = floor(j/ni) + 1;
	  i = ipt_ind_list[ind] -  (k-1)*ni*nj - (j-1)*ni + 1;



	  if (i == 2 and j > 2 and j < nj- 2) // bande imin 
	    {
	      a = ipt_coefflist[indCoef]   ;
	      b = ipt_coefflist[indCoef+2] ;
	      ipt_coefflist[indCoef+1] = ipt_coefflist[indCoef+1] + a/2.0 + b/2.0 ;	      
	      ipt_coefflist[indCoef+3] = ipt_coefflist[indCoef+3] + a/2.0 + b/2.0 ;
	      ipt_coefflist[indCoef] = 0.0   ;
	      ipt_coefflist[indCoef+2] = 0.0 ;
	    }


	  if (i == ni-2 and j > 2 and j < nj- 2) // bande imax
	    {
	      a = ipt_coefflist[indCoef+1]   ;
	      b = ipt_coefflist[indCoef+3] ;
	      ipt_coefflist[indCoef] = ipt_coefflist[indCoef] + a/2.0 + b/2.0 ;	      
	      ipt_coefflist[indCoef+2] = ipt_coefflist[indCoef+2] + a/2.0 + b/2.0;
	      ipt_coefflist[indCoef+1] = 0.0   ;
	      ipt_coefflist[indCoef+3] = 0.0 ;
	    }


	  if (j == 2 and i > 2 and i < ni- 2) // bande jmin 
	    {
	      a = ipt_coefflist[indCoef]   ;
	      b = ipt_coefflist[indCoef+1] ;
	      ipt_coefflist[indCoef+2] = ipt_coefflist[indCoef+2] + a/2.0 + b/2.0 ;	      
	      ipt_coefflist[indCoef+3] = ipt_coefflist[indCoef+3] + a/2.0 + b/2.0;
	      ipt_coefflist[indCoef] = 0.0   ;
	      ipt_coefflist[indCoef+1] = 0.0 ;


	    }


	  if (j == nj-2 and i > 2 and i < ni- 2) // bande jmax
	    {
	      a = ipt_coefflist[indCoef+2] ; 
	      b = ipt_coefflist[indCoef+3] ;
	      ipt_coefflist[indCoef] = ipt_coefflist[indCoef] + a/2.0 + b/2.0 ;	      
	      ipt_coefflist[indCoef+1] = ipt_coefflist[indCoef+1] + a/2.0 + b/2.0;
	      ipt_coefflist[indCoef+2] = 0.0 ;
	      ipt_coefflist[indCoef+3] = 0.0 ;
	    }


	  if (i == 2 and j == 2) // coin bas gauche 
	    {      
	      ipt_coefflist[indCoef+3] = 1.0 ;
	      ipt_coefflist[indCoef+1] = 0.0 ;
	      ipt_coefflist[indCoef]   = 0.0 ;
	      ipt_coefflist[indCoef+2] = 0.0 ;
	    }


	  if (i == 2 and j ==  nj- 2) // coin haut gauche
	    {
	      ipt_coefflist[indCoef+3] = 0.0 ;
	      ipt_coefflist[indCoef+1] = 1.0 ;
	      ipt_coefflist[indCoef]   = 0.0 ;
	      ipt_coefflist[indCoef+2] = 0.0 ;

	    }


	  if (i == ni-2 and j == 2) // coin bas droite 
	    {
	      ipt_coefflist[indCoef+3] = 0.0 ;
	      ipt_coefflist[indCoef+1] = 0.0 ;
	      ipt_coefflist[indCoef]   = 0.0 ;
	      ipt_coefflist[indCoef+2] = 1.0 ;

	    }


	  if (i == ni-2 and j == nj- 2) // coin haut droite
	    {
	      ipt_coefflist[indCoef+3] = 0.0 ;
	      ipt_coefflist[indCoef+1] = 0.0 ;
	      ipt_coefflist[indCoef]   = 1.0 ;
	      ipt_coefflist[indCoef+2] = 0.0 ;

	    }

	  indCoef += 4 ;

     }

   else if (ipt_typ[ind]==1)

     {
	  indCoef += 1;
     }

   else if (ipt_typ[ind]==2)
 
     { 
   
 	   k = floor(ipt_ind_list[ind]/(ni*nj)) + 1;
	   j = ipt_ind_list[ind] - (k-1)*ni*nj;
	   j = floor(j/ni)+1;
	   i = ipt_ind_list[ind] - (k-1)*ni*nj - (j-1)*ni + 1;

 
	   // Centre du cube 0
	   xc = float(2*i + 1)/2.0;
	   yc = float(2*j + 1)/2.0;
	   zc = float(2*k + 1)/2.0;
	   tab[0][0]=xc;
	   tab[0][1]=yc;
	   tab[0][2]=zc;

	   // Centre du cube 1
	   xc = float(2*i + 3)/2.0;
	   yc = float(2*j + 1)/2.0;
	   zc = float(2*k + 1)/2.0;
	   tab[1][0]=xc;
	   tab[1][1]=yc;
	   tab[1][2]=zc;

	   // Centre du cube 2
	   xc = float(2*i + 1)/2.0;
	   yc = float(2*j + 3)/2.0;
	   zc = float(2*k + 1)/2.0;
	   tab[2][0]=xc;
	   tab[2][1]=yc;
	   tab[2][2]=zc;

	   // Centre du cube 3
	   xc = float(2*i + 3)/2.0;
	   yc = float(2*j + 3)/2.0;
	   zc = float(2*k + 1)/2.0;
	   tab[3][0]=xc;
	   tab[3][1]=yc;
	   tab[3][2]=zc;

	   // Centre du cube 4
	   xc = float(2*i + 1)/2.0;
	   yc = float(2*j + 1)/2.0;
	   zc = float(2*k + 3)/2.0;
	   tab[4][0]=xc;
	   tab[4][1]=yc;
	   tab[4][2]=zc;

	   // Centre du cube 5
	   xc = float(2*i + 3)/2.0;
	   yc = float(2*j + 1)/2.0;
	   zc = float(2*k + 3)/2.0;
	   tab[5][0]=xc;
	   tab[5][1]=yc;
	   tab[5][2]=zc;

	   // Centre du cube 6
	   xc = float(2*i + 1)/2.0;
	   yc = float(2*j + 3)/2.0;
	   zc = float(2*k + 3)/2.0;
	   tab[6][0]=xc;
	   tab[6][1]=yc;
	   tab[6][2]=zc;

	   // Centre du cube 7
	   xc = float(2*i + 3)/2.0;
	   yc = float(2*j + 3)/2.0;
	   zc = float(2*k + 3)/2.0;
	   tab[7][0]=xc;
	   tab[7][1]=yc;
	   tab[7][2]=zc;

	   for (E_Int z = 0; z < 8; z++)
	     {
	       if (tab[z][0] > 3.0 and tab[z][0] < float(ni-1) and tab[z][1] > 3.0 and tab[z][1] < float(nj-1) and tab[z][2] > 3.0 and tab[z][2] < float(nk-1))
	   	 {
	   	   // cube ok
	   	   tab_[z]=-1;
	   	   // compte le nombre de cubes ok
	   	   compt = compt + 1;
	   	 } 
	       else
	   	 {
   	   	   // cube pas ok
	   	   tab_[z]=0;
	   	 }
	       //  cout << tab_[z] << endl;

	     }	  

	   // if (compt==0)
	   // {
	   //   cout <<"Danger"<< endl;
	   // }


	   if (compt == 1 or compt == 2 or compt == 3 or compt == 4 or compt == 5 or compt == 6 or compt == 7) // On corrige les coeffs des molecules qui contiennent les deux types de cellules
	     {
	       for ( j=0; j < 8; j++)
	   	 {
	   	   if (tab_[j]==0) // On tombe sur un cube pas ok
	   	     {
	   	       for ( E_Int z=0; z < 8; z++)// On repartit le coeff du cube pas ok sur les cubes ok 
	   		 {
	   		   if (tab_[z]==-1) 
	   		     {
	   		       ipt_coefflist[indCoef+z] = ipt_coefflist[indCoef+z] + ipt_coefflist[indCoef+j]/float(compt) ;

	   		     }


	   		 }
   	   	       ipt_coefflist[indCoef+j] = 0.0 ;


	   	     }


	   	 }



	     }

	   //  for ( j=0; j < 8; j++)
	   //	{

	   // S = S + ipt_coefflist[indCoef+j]; 
	   //	}
	   //cout << i <<" "<<j<<" "<<k<<" "<< ipt_coefflist[indCoef] << endl;
	   //cout << i+1 <<" "<<j<<" "<<k<<" "<< ipt_coefflist[indCoef+1] << endl;
	   //cout << i <<" "<<j+1<<" "<<k<<" "<< ipt_coefflist[indCoef+2] << endl;
	   //cout << i+1 <<" "<<j+1<<" "<<k<<" "<< ipt_coefflist[indCoef+3] << endl;
	   //cout << i <<" "<<j<<" "<<k+1<<" "<< ipt_coefflist[indCoef+4] << endl;
	   //cout << i+1 <<" "<<j<<" "<<k+1<<" "<< ipt_coefflist[indCoef+5] << endl;
	   //cout << i <<" "<<j+1<<" "<<k+1<<" "<< ipt_coefflist[indCoef+6] << endl;
	   //cout << i+1 <<" "<<j+1<<" "<<k+1<<" "<< ipt_coefflist[indCoef+7] << endl;


	   //   cout << S<< endl;
	   //   S=0;
	   indCoef += 8;
	   compt=0;
	  
        
     }

   

   }





 RELEASESHAREDN( indiceslist  , ind_list );
 RELEASESHAREDN( coefflist  , coefflist_ );

 Py_INCREF(Py_None);
 return Py_None;

}
