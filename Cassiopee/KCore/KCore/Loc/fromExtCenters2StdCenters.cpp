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
# include <stdio.h>

using namespace std;
using namespace K_FLD;
//=============================================================================
/* Retourne indTab sous la forme [i,i+1, ..., j,j+1] des centres 
   a partir de la cellule indExt en centres etendus
   Si une des cellules en centres est au bord extrapB = 1
   indTab est alloue ici */
//=============================================================================
short K_LOC::fromExtCenters2StdCenters( 
  E_Int ime, E_Int jme, E_Int kme, 
  E_Int indExt, E_Int type, FldArrayI& indTab,
  E_Int& extrapB)
{
  E_Int imr = ime-2;//maillage en noeuds
  E_Int jmr = jme-2;//maillage en noeuds
  E_Int kmr = kme-2;//maillage en noeuds
  E_Int imejme = ime*jme;
  E_Int kc = indExt/imejme; 
  E_Int jc = (indExt-kc*imejme)/ime;
  E_Int ic = indExt-jc*ime-kc*imejme;
  ic += 1; jc+= 1; kc+=1;

  E_Int inci1, inci2, incj1, incj2, inck1, inck2;
  E_Int inci3, inci4, incj3, incj4, inck3, inck4;
  extrapB = 0;

  switch (type)
  {
    case 1:
      indTab.malloc(3);
      ic = ic-2; jc = jc-2; kc = kc-2;
      inci1 = 1; incj1 = 1; inck1 = 1;
      // Dim  = 2
      if ( imr == 1 ) 
      {ic = 0; inci1 = 0;}    
      else if ( jmr == 1 )
      {jc = 0;incj1 = 0; }
      else if ( kmr == 1 )
      {kc = 0; inck1 = 0;}
      indTab[0] = ic;
      indTab[1] = jc;
      indTab[2] = kc;
      //direction i : bords
      if (ic == -1)
      {indTab[0] = 0; extrapB = 1;}
       
      // direction j : bords
      if (jc == -1)
      {indTab[1] = 0; extrapB = 1;}

      //direction k : bords
      if (kc == -1)
      {indTab[2] = 0; extrapB = 1;}
      break;

    case 2:
      indTab.malloc(6);
      ic = ic-2; jc = jc-2; kc = kc-2;
      inci1 = 1; incj1 = 1; inck1 = 1;
      // Dim  = 2
      if ( imr == 1 ) 
      {ic = 0; inci1 = 0;}    
      else if ( jmr == 1 )
      {jc = 0;incj1 = 0; }
      else if ( kmr == 1 )
      {kc = 0; inck1 = 0;}
      indTab[0] = ic;
      indTab[1] = ic+inci1;
      indTab[2] = jc;
      indTab[3] = jc+incj1;
      indTab[4] = kc;
      indTab[5] = kc+inck1;

      //direction i : bords
      if (ic == -1)
      {indTab[0] = 0; extrapB = 1;}
      
      else if (ic == imr-1)
      {indTab[1] = ic; extrapB = 1;}
 
      // direction j : bords
      if (jc == -1)
      {indTab[2] = 0; extrapB = 1;}

      else if (jc == jmr-1)
      {indTab[3] = jc; extrapB = 1;}

      //direction k : bords
      if (kc == -1)
      {indTab[4] = 0; extrapB = 1;}
      else if (kc == kmr-1)
      {indTab[5] = kc; extrapB = 1;} 

      break;
    
    case 3:
    indTab.malloc(9);

    ic = ic-2; jc = jc-2; kc = kc-2;
    inci1 = 1; inci2 = 2;
    incj1 = 1; incj2 = 2;
    inck1 = 1; inck2 = 2;
    
    // Dim  = 2
    if ( imr == 1 ) 
    {ic = 0; inci1 = 0; inci2 = 0;}    
    else if ( jmr == 1 )
    {jc = 0;incj1 = 0; incj2 = 0;}
    else if ( kmr == 1 )
    {kc = 0; inck1 = 0; inck2 = 0;}
   
    indTab[0] = ic;
    indTab[1] = ic+inci1;
    indTab[2] = ic+inci2;
    indTab[3] = jc;
    indTab[4] = jc+incj1;
    indTab[5] = jc+incj2;
    indTab[6] = kc;
    indTab[7] = kc+inck1;
    indTab[8] = kc+inck2; 

    //direction i : bords 
    if (ic == -1)
    {indTab[0] = 0; extrapB = 1;}
    
    else if (ic == imr-2)
    {indTab[2] = ic+1; extrapB = 1;}
 
    else if ( ic == imr-1)
    {
      indTab[1] = ic;
      indTab[2] = ic;
      extrapB = 1;
    }

    //direction j : bords 
    if (jc == -1)
    {indTab[3] = 0; extrapB = 1;}
    
    else if (jc == jmr-2)
    {indTab[5] = jc+1; extrapB = 1;}
    
    else if ( jc == jmr-1)
    {
      indTab[4] = jc;
      indTab[5] = jc;
      extrapB = 1;
    }

    //direction k : bords 
    if (kc == -1)
    {indTab[6] = 0; extrapB = 1;}
    
    else if (kc == kmr-2)
    {indTab[8] = kc+1; extrapB = 1;}
    
    else if ( kc == kmr-1)
    {
      indTab[7] = kc;
      indTab[8] = kc;
      extrapB = 1;
    }
    break;

    case 5: 
      indTab.malloc(15);
      ic = ic-2;  jc = jc-2; kc = kc-2;
      inci1 = 1; incj1 = 1; inck1 = 1;
      inci2 = 2; incj2 = 2; inck2 = 2;
      inci3 = 3; incj3 = 3; inck3 = 3;
      inci4 = 4; incj4 = 4; inck4 = 4;

      // Dimension 2
      if ( imr == 1 ) 
      {
        ic = 0;
        inci1 = 0;
        inci2 = 0;
        inci3 = 0;
        inci4 = 0;
      }
     
      else if ( jmr == 1 )
      {
        jc = 0;
        incj1 = 0;
        incj2 = 0;
        incj3 = 0;
        incj4 = 0;
      }
      else if ( kmr == 1 )
      {
        kc = 0;   
        inck1 = 0;
        inck2 = 0;  
        inck3 = 0;
        inck4 = 0;
      }

      indTab[0] = ic;
      indTab[1] = ic+inci1;
      indTab[2] = ic+inci2;
      indTab[3] = ic+inci3;
      indTab[4] = ic+inci4;
      indTab[5] = jc;
      indTab[6] = jc+incj1;
      indTab[7] = jc+incj2;
      indTab[8] = jc+incj3;
      indTab[9] = jc+incj4;
      indTab[10] = kc;
      indTab[11] = kc+inck1;
      indTab[12] = kc+inck2; 
      indTab[13] = kc+inck3;
      indTab[14] = kc+inck4;
 
      //direction i : bords 
      if (ic == -1)
      {indTab[0] = 0; extrapB = 1;}
      
      else if ( ic == imr-4)
      {indTab[4] = imr-1; extrapB = 1;}

      else if ( ic == imr-3)
      {     
        indTab[3] = imr-1;
        indTab[4] = imr-1;
        extrapB = 1;
      }
      else if ( ic == imr-2)
      {
        indTab[2] = imr-1;
        indTab[3] = imr-1;
        indTab[4] = imr-1;
        extrapB = 1;
      }
      else if ( ic == imr-1)
      {
        indTab[1] = ic;
        indTab[2] = ic;
        indTab[3] = ic;
        indTab[4] = ic;
        extrapB = 1;
      }
      
      //direction j : bords 
      if (jc == -1)
      {indTab[5] = 0;extrapB = 1;}
      
      else if ( jc == jmr-4)
      {indTab[9] = jmr-1;extrapB = 1;}

      else if ( jc == jmr-3)
      {     
        indTab[8] = jmr-1;
        indTab[9] = jmr-1;
        extrapB = 1;
      }
      else if ( jc == jmr-2)
      {
        indTab[7] = jmr-1;
        indTab[8] = jmr-1;
        indTab[9] = jmr-1;
        extrapB = 1;
      }
      else if ( jc == jmr-1)
      {
        indTab[6] = jc;
        indTab[7] = jc;
        indTab[8] = jc;
        indTab[9] = jc;
        extrapB = 1;
      }
      
      //direction k : bords 
      if (kc == -1)
      {indTab[10] = 0; extrapB = 1;}
      
      else if ( kc == kmr-4)
      {indTab[14] = kmr-1; extrapB = 1;}

      else if ( kc == kmr-3)
      {     
        indTab[13] = kmr-1;
        indTab[14] = kmr-1;
        extrapB = 1;
      }
      else if ( kc == kmr-2)
      {
        indTab[12] = kmr-1;
        indTab[13] = kmr-1;
        indTab[14] = kmr-1;
        extrapB = 1;
      }
      else if ( kc == kmr-1)
      {
        indTab[11] = kc;
        indTab[12] = kc;
        indTab[13] = kc;
        indTab[14] = kc;
        extrapB = 1;
      }
      break;

    default:
      printf("Interpolation type is unknown.\n");
      return -1;
  }
  return 1;
}
