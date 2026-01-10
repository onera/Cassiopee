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

# include "converter.h"

using namespace K_FLD;
using namespace std;

# define BLOCK                                                     \
  vdir1[0] = xt[inddir1]-xt[ind]; vdir2[0] = xt[inddir2]-xt[ind];       \
  vdir1[1] = yt[inddir1]-yt[ind]; vdir2[1] = yt[inddir2]-yt[ind];       \
  vdir1[2] = zt[inddir1]-zt[ind]; vdir2[2] = zt[inddir2]-zt[ind];       \
  vx = vdir1[0]+vdir2[0]; vy = vdir1[1]+vdir2[1]; vz = vdir1[2]+vdir2[2]; \
  if (K_FUNC::fEqualZero(vdir1[0]-vdir2[0]) == true &&                    \
  K_FUNC::fEqualZero(vdir1[1]-vdir2[1]) == true &&                    \
      K_FUNC::fEqualZero(vdir1[2]-vdir2[2]) == true )                     \
  {vx = vx/2.; vy = vy/2.; vz = vz/2.;}                                   \
  xt[indgc] = xt[ind]+vx; yt[indgc] = yt[ind]+vy; zt[indgc] = zt[ind]+vz;

//=============================================================================
/**/
//=============================================================================
PyObject* K_CONVERTER::fillCornerGhostCells2(PyObject* self, PyObject* args)
{
  //E_Int err = 0;
  PyObject* zone;
  E_Int ngc;
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  if (!PYPARSETUPLE_(args, O_ I_ SSS_, &zone, &ngc, 
                     &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
    return NULL;
  /* zone a modifier */
  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  E_Int res = K_PYTREE::getFromZone(zone, 1, 0, varString,
                                    fields, locs, im, jm, km,
                                    cn, cnSize, cnNfld, eltType, hook,
                                    GridCoordinates, 
                                    FlowSolutionNodes, FlowSolutionCenters);
  if (res != 1) 
  {
    PyErr_SetString(PyExc_TypeError,"fillCornerGhostCells2: zone must be structured.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  E_Int dim = 3;
  if ( im > 1 && jm > 1 ) 
  {
    if ( km > 1 ) dim = 3;
    else dim = 2;
  }
  else // on ne fait rien
  {
    RELEASESHAREDZ(hook, varString, eltType);
    Py_INCREF(Py_None);
    return Py_None;
  }

   E_Int posx = K_ARRAY::isCoordinateXPresent(varString);      
   E_Int posy = K_ARRAY::isCoordinateYPresent(varString);      
   E_Int posz = K_ARRAY::isCoordinateZPresent(varString);    
   if (posx == -1 || posy == -1 || posz == -1) 
   {
     PyErr_SetString(PyExc_TypeError,"fillCornerGhostCells2: coordinates not found.");
     RELEASESHAREDZ(hook, varString, eltType);
     return NULL;
   }
   E_Float* xt = fields[posx];
   E_Float* yt = fields[posy];
   E_Float* zt = fields[posz];
   // 8 aretes a traiter en premier, ensuite les coins
   E_Int indgc, ind, inddir1, inddir2;
   E_Int imjm = im*jm;
   E_Float vdir1[3]; E_Float vdir2[3];
   E_Float vx, vy, vz;
   
   if ( dim == 3 ) 
   {
     /* ARETES AD, EH, BC, FG */
     for (E_Int j = ngc; j < jm-ngc; j++)
     {
       // 1. arete AD (i=1;k=1) 
       for (E_Int k = ngc; k > 0; k--)
         for (E_Int i = ngc; i > 0; i--)
         {
           ind = i + j*im + k*imjm; inddir1 = ind-1; inddir2 = ind-imjm;
           indgc = (i-1) + j*im + (k-1)*imjm;
           BLOCK;
         }
     
       // 2. arete EH (i=1;k=kmax)
       for (E_Int k = km-ngc; k < km; k++)
         for (E_Int i = ngc; i > 0; i--)
         {
           ind = i + j*im + (k-1)*imjm; inddir1 = ind-1; inddir2 = ind+imjm;
           indgc = (i-1) + j*im + k*imjm;
           BLOCK;
         }

       // 3. arete BC(i=imax,k=1)
       for (E_Int k = ngc; k > 0; k--)
         for (E_Int i = im-ngc; i < im; i++)
         {
           ind = i-1 + j*im + k*imjm; inddir1 = ind+1; inddir2 = ind-imjm;
           indgc = i + j*im + (k-1)*imjm;
           BLOCK;
         }

       // 4. arete FG (i=imax,k=kmax)
       for (E_Int k = km-ngc; k < km; k++)
         for (E_Int i = im-ngc; i < im; i++)
         {
           ind = i-1 + j*im + (k-1)*imjm; inddir1 = ind+1; inddir2 = ind+imjm;
           indgc = i + j*im + k*imjm;
           BLOCK;         
         }     
     }
     /* ARETES AB, CD, EF, GH */
     for (E_Int i = ngc; i < im-ngc; i++)
     {
       //1. arete AB (j=1,k=1)
       for (E_Int k = ngc; k > 0; k--)
         for (E_Int j = ngc; j > 0; j--)
         {
           ind = i + j*im + k*imjm; inddir1 = ind-im; inddir2 = ind-imjm;       
           indgc = i + (j-1)*im + (k-1)*imjm;
           BLOCK;        
         }
       //2. arete DC (j=1,k=km)
       for (E_Int k = km-ngc; k < km; k++)
         for (E_Int j = ngc; j > 0; j--)
         {
           ind = i + j*im + (k-1)*imjm; inddir1 = ind-im; inddir2 = ind+imjm;
           indgc = i + (j-1)*im + k*imjm;
           BLOCK;
         }
       //3. arete EF(j=jm,k=1)
       for (E_Int k = ngc; k > 0; k--)
         for (E_Int j = jm-ngc; j < jm; j++)
         {
           ind = i + (j-1)*im + k*imjm; inddir1 = ind+im; inddir2 = ind-imjm;         
           indgc = i + j*im + (k-1)*imjm;
           BLOCK;
         }
       //4. arete HG (j=jm,k=km)
       for (E_Int k = km-ngc; k < km; k++)
         for (E_Int j = jm-ngc; j < jm; j++)
         {
           ind = i + (j-1)*im + (k-1)*imjm; inddir1 = ind+im; inddir2 = ind+imjm;
           indgc = i + j*im + k*imjm;
           BLOCK;
         }
     }
     /* ARETES AE, BF, DH, CG */
     for (E_Int k = ngc; k < km-ngc; k++)
     {
       // 1. arete AE (i=1,j=1)
       for (E_Int j = ngc; j > 0; j--)
         for (E_Int i = ngc; i > 0; i--)
         {
           ind = i + j*im + k*imjm; inddir1 = ind-1; inddir2 = ind-im;         
           indgc = (i-1) + (j-1)*im + k*imjm;
           BLOCK;
         }
       // 2. arete BF (i=imax,j=1)
       for (E_Int j = ngc; j > 0; j--)
         for (E_Int i = im-ngc; i < im; i++)
         {
           ind = i-1 + j*im + k*imjm; inddir1 = ind+1; inddir2 = ind-im;        
           indgc = i + (j-1)*im + k*imjm;
           BLOCK;
         }
       //3. arete DH (i=1,j=jmax)
       for (E_Int j = jm-ngc; j < jm; j++)
         for (E_Int i = ngc; i > 0; i--)
         {
           ind = i + (j-1)*im + k*imjm; inddir1 = ind-1; inddir2 = ind+im;
           indgc = (i-1) + j*im + k*imjm;
           BLOCK;
         }
       //4. arete CG (i=imax,j=jmax)
       for (E_Int j = jm-ngc; j < jm; j++)
         for (E_Int i = im-ngc; i < im; i++)
         {
           ind = (i-1) + (j-1)*im + k*imjm; inddir1 = ind+1; inddir2 = ind+im;
           indgc = i + j*im + k*imjm;
           BLOCK;
         }
     }
     /* 8 COINS */
     //E_Int i0, j0, k0;

     // 1. coin (1,1,1)
     for (E_Int j = 0; j < ngc; j++)
     {
       for (E_Int k = ngc; k > 0; k--)
         for (E_Int i = ngc; i > 0; i--)
         {
           ind = i + j*im + k*imjm; inddir1 = ind-1; inddir2 = ind-imjm;
           indgc = (i-1) + j*im + (k-1)*imjm;
           BLOCK;
         }
       //2. coin (1,1,km)
       for (E_Int k = km-ngc; k <km; k++)
         for (E_Int i = ngc; i > 0; i--)
         {
           ind = i + j*im + (k-1)*imjm; inddir1 = ind-1; inddir2 = ind+imjm;
           indgc = (i-1) + j*im + k*imjm;
           BLOCK;
         }
       //3. coin (im,1,1)
       for (E_Int k = ngc; k > 0; k--)
         for (E_Int i = im-ngc; i < im; i++)
         {
           ind = i-1 + j*im + k*imjm; inddir1 = ind+1; inddir2 = ind-imjm;
           indgc = i + j*im + (k-1)*imjm;
           BLOCK;
         }
       //4. coin (im,1,km)
       for (E_Int k = km-ngc; k < km; k++)
         for (E_Int i = im-ngc; i < im; i++)
         {
           ind = i-1 + j*im + (k-1)*imjm; inddir1 = ind+1; inddir2 = ind+imjm;
           indgc = i + j*im + k*imjm;
           BLOCK;         
         }     
     }
     for (E_Int j = jm-ngc; j < jm; j++)
     {
       // 5. coin (1,jm,1)
       for (E_Int k = ngc; k > 0; k--)
         for (E_Int i = ngc; i > 0; i--)
         {
           ind = i + j*im + k*imjm; inddir1 = ind-1; inddir2 = ind-imjm;
           indgc = (i-1) + j*im + (k-1)*imjm;
           BLOCK;
         }
       //6. coin (1,jm,km)
       for (E_Int k = km-ngc; k <km; k++)
         for (E_Int i = ngc; i > 0; i--)
         {
           ind = i + j*im + (k-1)*imjm; inddir1 = ind-1; inddir2 = ind+imjm;
           indgc = (i-1) + j*im + k*imjm;
           BLOCK;
         }
       //7. coin (im,jm,1)
       for (E_Int k = ngc; k > 0; k--)
         for (E_Int i = im-ngc; i < im; i++)
         {
           ind = i-1 + j*im + k*imjm; inddir1 = ind+1; inddir2 = ind-imjm;
           indgc = i + j*im + (k-1)*imjm;
           BLOCK;
         }
       //8. coin (im,jm,km)
       for (E_Int k = km-ngc; k < km; k++)
         for (E_Int i = im-ngc; i < im; i++)
         {
           ind = i-1 + j*im + (k-1)*imjm; inddir1 = ind+1; inddir2 = ind+imjm;
           indgc = i + j*im + k*imjm;
           BLOCK;         
         } 
     }
   }
   else  // cas 2D
   {
     //1. coin A (1,1,1)
     for (E_Int j = ngc; j > 0; j--)
       for (E_Int i = ngc; i > 0; i--)
       {
         ind = i + j*im; inddir1 = ind-1; inddir2 = ind-im;
         indgc = (i-1) + (j-1)*im;
         BLOCK;
       } 
     //2. coin B (im,1,1)
     for (E_Int j = ngc; j > 0; j--)
       for (E_Int i = im-ngc; i < im; i++)
       {
         ind = i-1 + j*im; inddir1 = ind+1; inddir2 = ind-im;
         indgc = i + (j-1)*im;
         BLOCK;
       } 
     //3. coin C (1,jm,1)
     for (E_Int j = jm-ngc; j < jm; j++)
       for (E_Int i = ngc; i > 0; i--)
       {
         ind = i + (j-1)*im; inddir1 = ind-1; inddir2 = ind+im;
         indgc = (i-1) + j*im;
         BLOCK;
       } 
     //4. coin D (im,jm,1)
     for (E_Int j = jm-ngc; j < jm; j++)
       for (E_Int i = im-ngc; i < im; i++)
       {
         ind = i-1 + (j-1)*im; inddir1 = ind+1; inddir2 = ind+im;
         indgc = i + j*im;
         BLOCK;
       }
   }
   RELEASESHAREDZ(hook, varString, eltType);
   Py_INCREF(Py_None);
   return Py_None;
}
