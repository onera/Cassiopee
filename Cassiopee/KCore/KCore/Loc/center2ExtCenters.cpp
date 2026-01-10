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

# include "loc.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Passage des données de centres aux centres étendus
   Les valeurs aux bords (centres des faces externes) sont extrapolées 
   des valeurs des centres adjacents
   Retourne 1 en cas de succès, 0 en cas d'échec*/
//=============================================================================
E_Int K_LOC::center2ExtCenterStruct(E_Int im, E_Int jm, E_Int km,
                                    FldArrayF& FCenter,
                                    E_Int ime, E_Int jme, E_Int kme, 
                                    FldArrayF& FExtCenter)
{
  E_Int i, j, k, indc, inde, ic, jc, kc;
  E_Int nfld = FCenter.getNfld();
  E_Int size = ime*jme*kme;

  if (FExtCenter.getSize() != size || FExtCenter.getNfld() != nfld) 
    FExtCenter.malloc(size, nfld);
  FExtCenter.setAllValuesAtNull();
  E_Int imejme = ime*jme;
  E_Int imjm = im*jm;
  if ( ime > 1 && jme > 1 && kme > 1 ) 
  {
    for (E_Int v = 1; v <= nfld; v++)
    {
      E_Float* fc = FCenter.begin(v);
      E_Float* fe = FExtCenter.begin(v);
      
      //centres
      for (k = 1; k < kme-1; k++) 
        for (j = 1; j < jme-1; j++) 
          for (i = 1; i < ime-1; i++) 
          {
            inde = i + j * ime + k*imejme;
            indc = (i-1) + (j-1) * im + (k-1)*imjm;
            fe[inde] = fc[indc];
          }
      /* plan i=0 : on extrapole les centres au dessus*/
      i = 0; ic = 0;
       for (k = 0; k < kme; k++) 
         for (j = 0; j < jme; j++) 
         {
           if (j == 0) jc = 0;
           else if (j == jme-1) jc = jm-1;
           else jc = j-1;

           if (k == 0) kc = 0;
           else if (k == kme-1) kc = km-1;
           else kc = k-1;

           inde = i + j * ime + k*imejme;
           indc = ic + jc * im + kc*imjm;
           fe[inde] = fc[indc];
         }       
       
       /* plan i=ime-1 : on extrapole les centres au dessus*/
       i = ime-1; ic = im-1;
       for (k = 0; k < kme; k++) 
         for (j = 0; j < jme; j++) 
         {
           if (j == 0) jc = 0;
           else if (j == jme-1) jc = jm-1;
           else jc = j-1;

           if (k == 0) kc = 0;
           else if (k == kme-1) kc = km-1;
           else kc = k-1;

           inde = i + j * ime + k*imejme;
           indc = ic + jc * im + kc*imjm;
           fe[inde] = fc[indc];
         } 
       /* plan j=0 : on extrapole les centres au dessus */
       j = 0; jc = 0;
       for (k = 0; k < kme; k++) 
         for (i = 0; i < ime; i++) 
         {
           if (i == 0) ic = 0;
           else if (i == ime-1) ic = im-1;
           else ic = i-1;

           if (k == 0) kc = 0;
           else if (k == kme-1) kc = km-1;
           else kc = k-1;
           inde = i + j * ime + k*imejme;
           indc = ic + jc * im + kc*imjm;
           fe[inde] = fc[indc];
         }       
       
       /* plan j=jme-1 : on extrapole les centres au dessus*/
       j = jme-1; jc = jm-1;
       for (k = 0; k < kme; k++) 
         for (i = 0; i < ime; i++) 
         {
           if (i == 0) ic = 0;
           else if (i == ime-1) ic = im-1;
           else ic = i-1;

           if (k == 0) kc = 0;
           else if (k == kme-1) kc = km-1;
           else kc = k-1;

           inde = i + j * ime + k*imejme;
           indc = ic + jc * im + kc*imjm;
           fe[inde] = fc[indc];
         }        
       /* plan k=0 : on extrapole les centres au dessus*/
       k = 0; kc = 0;
       for (j = 0; j < jme; j++) 
         for (i = 0; i < ime; i++) 
         {
           if (i == 0) ic = 0;
           else if (i == ime-1) ic = im-1;
           else ic = i-1;

           if (j == 0) jc = 0;
           else if (j == jme-1) jc = jm-1;
           else jc = j-1;

           inde = i + j * ime + k*imejme;
           indc = ic + jc * im + kc*imjm;
           fe[inde] = fc[indc];
         }       
       
       /* plan k=kme-1 : on extrapole les centres au dessus*/
       k = kme-1; kc = km-1;
       for (j = 0; j < jme; j++) 
         for (i = 0; i < ime; i++) 
         {
           if (i == 0) ic = 0;
           else if (i == ime-1) ic = im-1;
           else ic = i-1;

           if (j == 0) jc = 0;
           else if (j == jme-1) jc = jm-1;
           else jc = j-1;

           inde = i + j * ime + k*imejme;
           indc = ic + jc * im + kc*imjm;
           fe[inde] = fc[indc];
         }       
    }
  }
  return 1;
}
