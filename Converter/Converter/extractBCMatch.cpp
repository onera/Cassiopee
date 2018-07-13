/*    
    Copyright 2013-2018 Onera.

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
#include "converter.h"

#define signVar(a) (a < 0 ? -1 : 1)

using namespace std;
using namespace K_FLD;


//=============================================================================
//=============================================================================
PyObject* K_CONVERTER::extractBCMatch(PyObject* self, PyObject* args )
{
  PyObject *fields, *indBC, *fldBC;

  E_Int niD, njD, nkD;       // dim zone donneuse
  E_Int niR, njR, nkR;       // dim zone receveuse 

  E_Int iminD, jminD, kminD; // indices fenetre donneuse
  E_Int imaxD, jmaxD, kmaxD; // indices fenetre donneuse
  E_Int iminR, jminR, kminR; // indices fenetre receveuse
  E_Int imaxR, jmaxR, kmaxR; // indices fenetre receveuse 

  E_Int triI, triJ, triK ;   // transform (issu du GC de la zone "receveuse")

  if (!PYPARSETUPLEI(args, "O(llllll)(llllll)(lll)(lll)", "O(iiiiii)(iiiiii)(iii)(iii)", 
                     &fields, &iminD, &jminD, &kminD, &imaxD, &jmaxD, &kmaxD,
		              &iminR, &jminR, &kminR, &imaxR, &jmaxR, &kmaxR, 
		              &niR, &njR, &nkR, 
                              &triI, &triJ, &triK )) return NULL;

  // Check array
  // ===========
  FldArrayF* FCenter; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(fields, varString, FCenter, niD, njD, nkD, 
                                    cn, eltType); 

  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError, "extractBCMatch: array must be structured."); 
    if (res == 2) RELEASESHAREDS(fields, FCenter);
    return NULL; 
  }

  // 
  E_Int dim       = 3;
  E_Int noindint  = 0;
  E_Int ind ;
  E_Int nbIntID   = (niD+1)*njD*K_FUNC::E_max(1,nkD) ;
  E_Int nbIntJD   = (njD+1)*niD*K_FUNC::E_max(1,nkD) ;
  E_Int nbIntIR   = (niR+1)*njR*K_FUNC::E_max(1,nkR) ;
  E_Int nbIntJR   = (njR+1)*niR*K_FUNC::E_max(1,nkR) ;

  E_Int ifaceR ;
  E_Int jfaceR ;
  E_Int kfaceR ;
  E_Int shift  ;

  // compute dim 
  if ((niD == 1) or (njD == 1) or (nkD ==1)) 
  { 
    dim = 2; 
  }

  // printf("niD = %d ; njD = %d ; nkD = %d \n", niD, njD, nkD);

  // build output arrays 
  // ===================
  E_Int nfld = FCenter->getNfld();
  E_Int nint = max(1,(imaxD-iminD))*max(1,(jmaxD-jminD))*max(1,(kmaxD-kminD)); 

  // 1. tableau des indices 
  // ---------------------
  PyObject* indFace = K_NUMPY::buildNumpyArray(nint,1,1);
  E_Int* ptrIndFace = K_NUMPY::getNumpyPtrI(indFace);

  // 2. tableau des champs 
  // ---------------------
  PyObject* tpl = K_ARRAY::buildArray2(nfld,varString,nint,1,1,2); 

  FldArrayF*  fBC; 
  FldArrayI* cnBC;
  E_Int ni2, nj2, nk2;
  K_ARRAY::getFromArray2(tpl, varString, fBC, ni2, nj2, nk2, cnBC, eltType);

  // printf("trirac : %d %d %d \n",triI,triJ,triK);

  // Cas 2D
  // ======
  if (dim == 2)
  {
    // Face donneuse en i 
    // ******************
    if (iminD == imaxD) 
    { 
      // printf("Frontiere en i \n");
      // 1. tableau des indices
      // ----------------------

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) // A supprimer
        {
          E_Int indFaceD   = iminD - 1 + jface*(niD+1) ; // A supprimer
 
          // printf("indD : %d \n", indFaceD);// A supprimer
          noindint++;// A supprimer
        }// A supprimer

      // 1.a. face receveuse en i
      if (iminR == imaxR) 
      {
        // printf("Frontiere receveuse en i \n");
	noindint = 0;

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
        { 
          if (triJ > 0)
	  {
	    jfaceR =  jface + jminR - jminD ; // OK ! 
	  }
	  else
	  {
	    jfaceR = (jmaxR-2) - (jface-jminD+1); // OK !
	  }

          ptrIndFace[noindint] = iminR - 1 + jfaceR*(niR+1) ;
          // printf("** indR : %d \n", ptrIndFace[noindint]);
          noindint++;
	}
      }
      // 1.b. face receveuse en j
      if (jminR == jmaxR) 
      {
        // printf("Frontiere receveuse en j \n");
	noindint = 0;

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
        { 
          if (triI > 0)
	  {
	    ifaceR = jface + iminR - jminD ;
	  }
	  else
	  {
	    ifaceR = (imaxR-2) - (jface-jminD+1)  ; 
	  }

          shift  = (jminR-1)*niR + nbIntIR ; 
          ptrIndFace[noindint] = shift + ifaceR ;
          // printf("indR   : %d \n", ptrIndFace[noindint]);
          noindint++;
	}
      } // jminR=jmaxR

      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0 ;

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
        {
          if (iminD==1) { ind = jface*niD         ; } 
          else          { ind = jface*niD + niD-1 ; } 

          fld[noindint] = fce[ind] ; 
          noindint++;
        } 
      } // var loop 
    }
    // Si frontiere en j
    // *****************
    if (jminD == jmaxD) 
    {      
      // printf("Frontiere en j \n");
      // 1. tableau des indices 
      E_Int shift = (jminD-1)*niD + nbIntID ;

      for (E_Int iface = iminD - 1 ; iface < imaxD-1 ; iface ++)  // A supprimer
      {
        E_Int indFaceD = shift + iface ; // A supprimer
        // printf("indD : %d \n", indFaceD);// A supprimer
        noindint++; // A supprimer
      }   // A supprimer

      // printf("indR : ");
      // 1.a. face receveuse en i
      if (iminR == imaxR) 
      {
        // printf("Frontiere receveuse en i \n");
	noindint = 0;

	for (E_Int iface = iminD - 1 ; iface < imaxD-1 ; iface ++)
        { 
          if (triJ > 0)
	  {
	    jfaceR =  iface + jminR - iminD ; 
	  }
	  else
	  {
	    jfaceR = (jmaxR-2) - (iface-iminD+1);
	  }

          ptrIndFace[noindint] = iminR - 1 + jfaceR*(niR+1) ;
          // printf(" %d ", ptrIndFace[noindint]);
          noindint++;
	}
      } // iminR==imaxR 

      // 1.b. face receveuse en j
      if (jminR == jmaxR) 
      {
        // printf("Frontiere receveuse en j \n");
	noindint = 0;

	for (E_Int iface = iminD - 1 ; iface < imaxD-1 ; iface ++)
        { 
          E_Int shift  = (jminR-1)*niR + nbIntIR ;

          if (triI > 0)
	  {
	    ifaceR = iface + iminR - iminD ;
	  }
	  else
	  {
	    ifaceR = (imaxR-2) - (iface-iminD+1) ;
	  }

          ptrIndFace[noindint] = shift + ifaceR ;
          // printf(" %d ", ptrIndFace[noindint]);
          noindint++;
	}
      }
  
      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0 ;

        for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
        {
          if (jminD==1) { ind = iface               ; } 
          else          { ind = iface + niD*(njD-1) ; } 

          fld[noindint] = fce[ind] ; 
          noindint++;
        }
      } // var loop
    }// (si frontiere en j)
 
  } //(si dim=2)

  // Cas 3D
  // ======
  else if (dim == 3)
  {
    // printf("kminR : %d \n", kminR); 
    // ********************
    // Frontiere donneuse i 
    // ********************
    if (iminD == imaxD)
    {
      // printf("Frontiere donneuse en i \n");
      noindint = 0 ;
      // printf("indD : ");// A supprimer

      // 1. tableau des indices  // A supprimer 
      for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)  // A supprimer 
      {
        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)  // A supprimer 
        {
          E_Int indFaceD = iminD - 1 + jface*(niD+1) + kface*(niD+1)*njD ; // A supprimer 
          // printf("%d ", indFaceD);// A supprimer
          noindint++; // A supprimer 
        } 
      } 
      // printf("\n ");// A supprimer
 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en i 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (iminR == imaxR)
      {
        // printf("Frontiere receveuse en i \n");
        noindint = 0 ;

        if (abs(triJ)==2) // kD <-> kR et  jD <-> jR
	{
          // printf("Face receveuse jD=jR, kD=kR) \n");
          // printf("indR : ");
          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triK > 0) { kfaceR = kface + kminR - kminD ;     }
	    else          { kfaceR = (kmaxR-2)-(kface-kminD+1) ; }

            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
            {
              if (triJ > 0){ jfaceR = jface + jminR - jminD ;    }
	      else         { jfaceR = (jmaxR-2)-(jface-jminD+1); }

              ptrIndFace[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triJ=2

        if (abs(triJ)==3) // kD <-> jR et  jD <-> kR
	{
          // printf("Face receveuse jD=kR, kD=jR) \n");
          // printf("indR : ");

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)  
            { 
              if (triK > 0) { kfaceR = jface + kminR - jminD ;     }
	      else          { kfaceR = (kmaxR-2)-(jface-jminD+1) ; }

              if (triJ > 0) { jfaceR = kface + jminR - kminD ;     }
	      else          { jfaceR = (jmaxR-2)-(kface-kminD+1);  }
	    
            ptrIndFace[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
            // printf("%d ", ptrIndFace[noindint]);
            noindint++;
	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }
     
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en j 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (jminR == jmaxR)
      {
        // printf("Frontiere receveuse en j \n");
        noindint = 0 ;

        if (abs(triI)==2) // jD <-> iR et  kD <-> kR
	{
          // printf("Face receveuse jD=iR, kD=kR) \n");
          // printf("indR : ");
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triK > 0) { kfaceR = kface + kminR - kminD ;     }
	    else          { kfaceR = (kmaxR-2)-(kface-kminD+1) ; }

            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
            {
              if (triI > 0){ ifaceR = jface + iminR - jminD ;    }
	      else         { ifaceR = (imaxR-2)-(jface-jminD+1); }

              ptrIndFace[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triI=1

        if (abs(triI)==3) // jD <-> kR et  kD <-> iR
	{
          // printf("Face receveuse jD=kR, kD=iR) \n");
          // printf("indR : ");
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)  
            { 
              if (triK > 0) { kfaceR = jface + kminR - jminD ;     }
	      else          { kfaceR = (kmaxR-2)-(jface-jminD+1) ; }

              if (triI > 0) { ifaceR = kface + iminR - kminD ;     }
	      else          { ifaceR = (imaxR-2)-(kface-kminD+1);  }
	    
            ptrIndFace[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
            // printf("%d ", ptrIndFace[noindint]);
            noindint++;
	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }

    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en k
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (kminR == kmaxR)
      {
        // printf("Frontiere receveuse en k \n");
        noindint = 0 ;

        if (abs(triI)==2) // jD <-> iR et  kD <-> jR
      	{
          // printf("Face receveuse jD=iR, kD=jR) \n");
          // printf("indR : ");
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triJ > 0) { jfaceR = kface + jminR - kminD ;     }
      	    else          { jfaceR = (jmaxR-2)-(kface-kminD+1) ; }

            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
            {
              if (triI > 0){ ifaceR = jface + iminR - jminD ;    }
      	      else         { ifaceR = (imaxR-2)-(jface-jminD+1); }

              ptrIndFace[noindint] = shift + ifaceR + jfaceR*niR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
      	    }
          } 
        } // triI=1
        // printf("\n");

        if (abs(triI)==3) // jD <-> jR et  kD <-> iR
      	{
          // printf("Face receveuse jD=jR, kD=iR) \n");
          // printf("indR : ");
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            if (triI > 0) { ifaceR = kface + iminR - kminD ;     }
      	    else          { ifaceR = (imaxR-2)-(kface-kminD+1);  }

            for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)  
            { 
              if (triJ > 0) { jfaceR = jface + jminR - jminD ;     }
      	      else          { jfaceR = (jmaxR-2)-(jface-jminD+1) ; }
	    
            ptrIndFace[noindint] = shift + ifaceR + jfaceR*niR ;
            // printf("%d ", ptrIndFace[noindint]);
            noindint++;
      	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }

      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0 ;

        for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
        {
          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
          {
            if (iminD==1) { ind = iminD-1 + jface*niD + kface*njD*niD; } 
            else          { ind = niD-1    + jface*niD + kface*njD*niD; } 

            fld[noindint] = fce[ind] ; 
            noindint++;
          }
	}

      } // var loop 
    }
    // Si frontiere en j
    // *****************
    if (jminD == jmaxD)
    {
      // printf("Frontiere donneuse en j \n");
      noindint = 0 ;
      // printf("indD : ");// A supprimer

      // 1. tableau des indices 
      E_Int shift = (jminD-1)*niD + nbIntID ;

      for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) // A supprimer
      {
        for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) // A supprimer
        {
          E_Int indFaceD  = shift + iface + kface*niD*(njD+1) ;// A supprimer
          // printf("%d ", indFaceD);// A supprimer
          noindint++;// A supprimer
        }// A supprimer
      }
      // printf("\n ");// A supprimer

    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en i 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (iminR == imaxR)
      {
        // printf("Frontiere receveuse en i \n");
        noindint = 0 ;

        if (abs(triJ)==1) // kD <-> kR et  iD <-> jR
	{
          // printf("Face receveuse jD=jR, kD=kR) \n");
          // printf("indR : ");
          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triK > 0) { kfaceR = kface + kminR - kminD ;     }
	    else          { kfaceR = (kmaxR-2)-(kface-kminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triJ > 0){ jfaceR = iface + jminR - iminD ;    }
	      else         { jfaceR = (jmaxR-2)-(iface-iminD+1); }

              ptrIndFace[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triJ=2

        if (abs(triJ)==3) // kD <-> jR et  iD <-> kR
	{
          // printf("Face receveuse iD=kR, kD=jR) \n");
          // printf("indR : ");

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triK > 0) { kfaceR = iface + kminR - iminD ;     }
	      else          { kfaceR = (kmaxR-2)-(iface-iminD+1) ; }

              if (triJ > 0) { jfaceR = kface + jminR - kminD ;     }
	      else          { jfaceR = (jmaxR-2)-(kface-kminD+1);  }
	    
            ptrIndFace[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
            // printf("%d ", ptrIndFace[noindint]);
            noindint++;
	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }

    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en j 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (jminR == jmaxR)
      {
        // printf("Frontiere receveuse en j \n");
        noindint = 0 ;

        if (abs(triI)==1) // iD <-> iR et  kD <-> kR
	{
          // printf("Face receveuse iD=iR, kD=kR) \n");
          // printf("indR : ");
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triK > 0) { kfaceR = kface + kminR - kminD ;     }
	    else          { kfaceR = (kmaxR-2)-(kface-kminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triI > 0){ ifaceR = iface + iminR - iminD ;    }
	      else         { ifaceR = (imaxR-2)-(iface-iminD+1); }

              ptrIndFace[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triI=1

        if (abs(triI)==3) // iD <-> kR et  kD <-> iR
	{
          // printf("Face receveuse iD=kR, kD=iR) \n");
          // printf("indR : ");
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triK > 0) { kfaceR = iface + kminR - iminD ;     }
	      else          { kfaceR = (kmaxR-2)-(iface-iminD+1) ; }

              if (triI > 0) { ifaceR = kface + iminR - kminD ;     }
	      else          { ifaceR = (imaxR-2)-(kface-kminD+1);  }
	    
            ptrIndFace[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
            // printf("%d ", ptrIndFace[noindint]);
            noindint++;
	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }

    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en k 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (kminR == kmaxR)
      {
        // printf("Frontiere receveuse en k \n");
        noindint = 0 ;

        if (abs(triI)==1) // iD <-> iR et  kD <-> jR
	{
          // printf("Face receveuse iD=iR, kD=jR) \n");
          // printf("indR : ");
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
          {
            if (triJ > 0) { jfaceR = kface + jminR - kminD ;     }
	    else          { jfaceR = (jmaxR-2)-(kface-kminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triI > 0){ ifaceR = iface + iminR - iminD ;    }
	      else         { ifaceR = (imaxR-2)-(iface-iminD+1); }

              ptrIndFace[noindint] = shift + ifaceR + jfaceR*niR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triI=1

        if (abs(triI)==3) // iD <-> jR et  kD <-> iR
	{
          // printf("Face receveuse iD=jR, kD=iR) \n");
          // printf("indR : ");
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++)
          {
            if (triI > 0) { ifaceR = kface + iminR - kminD ;     }
	    else          { ifaceR = (imaxR-2)-(kface-kminD+1);  }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triJ > 0) { jfaceR = iface + jminR - iminD ;     }
	      else          { jfaceR = (jmaxR-2)-(iface-iminD+1) ; }
	    
              ptrIndFace[noindint] = shift + ifaceR + jfaceR*niR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }

      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0 ;

        for (E_Int kface = kminD-1 ; kface < kmaxD-1 ; kface ++) 
        {
          for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
          {
            if (jminD==1) { ind = iface + kface*njD*niD            ; } 
            else         { ind = iface + kface*njD*niD + (njD-1)*niD; } 

            fld[noindint] = fce[ind] ; 
            noindint++;
          }
	}

      } // var loop
    }
    // Si frontiere en k
    // *****************
    if (kminD == kmaxD)
    {
      // printf("Frontiere donneuse en k \n");
      // printf("indD : ");
      E_Int shift = (kminD-1)*niD*njD + nbIntID + nbIntJD ;

      for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
      {
        for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
        {
          ptrIndFace[noindint] = shift + iface + jface*niD ;
          // printf("%d ", ptrIndFace[noindint]);
          noindint++;
        }
      }
      // printf("\n ");

    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en i 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (iminR == imaxR)
      {
        // printf("Frontiere receveuse en i \n");
        noindint = 0 ;

        if (abs(triJ)==2) // iD <-> kR et  jD <-> jR
	{
          // printf("Face receveuse iD=kR, jD=jR) \n");
          // printf("indR : ");
          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
          {
            if (triJ > 0){ jfaceR = jface + jminR - jminD ;    }
	    else         { jfaceR = (jmaxR-2)-(jface-jminD+1); }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triK > 0) { kfaceR = iface + kminR - iminD ;     }
	      else          { kfaceR = (kmaxR-2)-(iface-iminD+1) ; }

              ptrIndFace[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triJ=2

        if (abs(triJ)==1) // iD <-> jR et  jD <-> kR
	{
          // printf("Face receveuse iD=jR, jD=kR) \n");
          // printf("indR : ");

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)
          {
            if (triK > 0) { kfaceR = jface + kminR - jminD ;     }
	    else          { kfaceR = (kmaxR-2)-(jface-jminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triJ > 0) { jfaceR = iface + jminR - iminD ;     }
	      else          { jfaceR = (jmaxR-2)-(iface-iminD+1);  }

              ptrIndFace[noindint] = iminR - 1 + jfaceR*(niR+1) + kfaceR*(niR+1)*njR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }

    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en j 
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (jminR == jmaxR)
      {
        // printf("Frontiere receveuse en j \n");
        noindint = 0 ;

        if (abs(triI)==2) // iD <-> kR et  jD <-> iR
	{
          // printf("Face receveuse iD=kR, jD=iR) \n");
          // printf("indR : ");
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
          {
            if (triI > 0){ ifaceR = jface + iminR - jminD ;    }
	    else         { ifaceR = (imaxR-2)-(jface-jminD+1); }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triK > 0) { kfaceR = iface + kminR - iminD ;     }
	      else          { kfaceR = (kmaxR-2)-(iface-iminD+1) ; }

              ptrIndFace[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triJ=2

        if (abs(triI)==1) // iD <-> iR et  jD <-> kR
	{
          // printf("Face receveuse iD=iR, jD=kR) \n");
          // printf("indR : ");
          E_Int shift = (jminR-1)*niR + nbIntIR ;

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)
          {
            if (triK > 0) { kfaceR = jface + kminR - jminD ;     }
	    else          { kfaceR = (kmaxR-2)-(jface-jminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triI > 0) { ifaceR = iface + iminR - iminD ;     }
	      else          { ifaceR = (imaxR-2)-(iface-iminD+1);  }

              ptrIndFace[noindint] = shift + ifaceR + kfaceR*niR*(njR+1) ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }
     
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Frontiere receveuse en k
    // ~~~~~~~~~~~~~~~~~~~~~~~~
      if (kminR == kmaxR)
      {
        // printf("Frontiere receveuse en k \n");
        noindint = 0 ;

        if (abs(triI)==1) // iD <-> iR et  jD <-> jR
	{
          // printf("Face receveuse iD=iR, jD=jR) \n");
          // printf("indR : ");
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
          {
            if (triJ > 0){ jfaceR = jface + jminR - jminD ;    }
	    else         { jfaceR = (jmaxR-2)-(jface-jminD+1); }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
            {
              if (triI > 0) { ifaceR = iface + iminR - iminD ;     }
	      else          { ifaceR = (imaxR-2)-(iface-iminD+1) ; }

              ptrIndFace[noindint] = shift + ifaceR + jfaceR*niR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triI=1

        if (abs(triI)==2) // iD <-> jR et  jD <-> iR
	{
          // printf("Face receveuse iD=jR, jD=iR) \n");
          // printf("indR : ");
          E_Int shift = (kminR-1)*niR*njR + nbIntIR + nbIntJR ;

          for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++)
          {
            if (triI > 0) { ifaceR = jface + iminR - jminD ;     }
	    else          { ifaceR = (imaxR-2)-(jface-jminD+1) ; }

            for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++)  
            { 
              if (triJ > 0) { jfaceR = iface + jminR - iminD ;     }
	      else          { jfaceR = (jmaxR-2)-(iface-iminD+1);  }

              ptrIndFace[noindint] = shift + ifaceR + jfaceR*niR ;
              // printf("%d ", ptrIndFace[noindint]);
              noindint++;
	    }
          } 
        } // triJ=3
 
        // printf("\n ");
      }

      // 2. tableau des champs 
      for (E_Int var = 1; var <= nfld; var++)
      {  
        E_Float* fld = fBC->begin(var);
        E_Float* fce = FCenter->begin(var);

        noindint = 0 ;

        for (E_Int jface = jminD-1 ; jface < jmaxD-1 ; jface ++) 
        {
          for (E_Int iface = iminD-1 ; iface < imaxD-1 ; iface ++) 
          {
            if (kminD==1) { ind = iface + jface*niD            ; } 
            else          { ind = iface + jface*niD + (nkD-1)*niD*njD; } 

            fld[noindint] = fce[ind] ; 
            noindint++;
          }
	}

      } // var loop
    }
  }


  RELEASESHAREDS(fields, FCenter);

  PyObject* tplOut ;
  tplOut = Py_BuildValue("[OO]",indFace,tpl);

  return tplOut; 
}


