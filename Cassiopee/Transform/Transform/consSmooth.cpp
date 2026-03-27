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

# include "transform.h"
using namespace K_FLD;
using namespace K_SEARCH;

//=============================================================================
/* Local functions */
//=============================================================================

inline bool relativeDist(E_Float p1x, E_Float p1y, E_Float p1z, E_Float nbg1x, E_Float nbg1y, E_Float nbg1z, E_Float p2x, E_Float p2y, E_Float p2z) 
{
  E_Float refL2 = (p1x - nbg1x)*(p1x - nbg1x) + (p1y - nbg1y)*(p1y - nbg1y) + (p1z - nbg1z)*(p1z - nbg1z);
  E_Float dist2 = (p1x - p2x)*(p1x - p2x) + (p1y - p2y)*(p1y - p2y) + (p1z - p2z)*(p1z - p2z);
  return (dist2 >= 1e-12 * refL2);
}

inline E_Int stepCorrection(E_Int isOpen, E_Int nUnique) 
{
  // Case of a CLOSED curve with an ODD number of points
  if ((isOpen == 0) && (nUnique % 2 != 0)) 
  {
    printf("Warning: consSmooth: step=2 is invalid for closed curve with odd number of points. Forcing step=1.\n");
    return 1;
  }
  
  // Case of an OPEN curve with an EVEN number of points
  if ((isOpen == 1) && (nUnique % 2 == 0)) 
  {
    printf("Warning: consSmooth: step=2 is invalid for open curve with even number of points. Forcing step=1.\n");
    return 1;
  }

  return 2;
}

// ============================================================================
/* Conservative smoothing (Kuprat)*/
// IN: array: mesh (Structured i-array ou unstructured 1D in the (x,y) plane)
// IN: sweeps: number of smoothing sweeps
// IN: twoWays: 0 = one way smoothing, 1 = one way and back (default 0, may lead to assymetry) 
// IN: step: way step (1 ou 2) default 1
// OUT: smoothed mesh (same as input)
// ============================================================================
PyObject* K_TRANSFORM::consSmooth(PyObject* self, PyObject* args)
{
  PyObject* array;

  E_Int sweeps = 1;
  E_Int step = 1;
  E_Int twoWays = 0;
  
  if (!PYPARSETUPLE_(args, O_ III_, &array, &sweeps,  &twoWays, &step))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =  K_ARRAY::getFromArray3(array, varString,
                                      f, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: array is invalid.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  if (twoWays != 0 && twoWays != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: int twoWays is invalid.");
    return NULL;
  }

  if (step < 1 && step > 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: int step is invalid.");
    return NULL;
  }

  // Build output
  PyObject* tpl;
  FldArrayF* fo;
  FldArrayI* cno = NULL;

  if (res == 1) // structured
  {
    tpl = K_ARRAY::buildArray3(*f, varString, im, jm, km, f->getApi());
    K_ARRAY::getFromArray3(tpl, fo);

    E_Float* x = fo->begin(posx);
    E_Float* y = fo->begin(posy);
    E_Float* z = fo->begin(posz);

    E_Int npts = im;
    
    E_Int idx0;
    E_Int idx1;
    E_Int idx2;
    E_Int idx3;

    // open or closed ?
    E_Int isOpen = relativeDist(x[0], y[0], z[0], x[1], y[1], z[1], x[npts-1], y[npts-1], z[npts-1]);

    if (isOpen)
    {
      printf("Info: consSmooth: open geometry: fixed nodes %d and %d.\n", 0, npts-1);
    }
    {
      printf("Info: consSmooth: closed geometry with double points.\n");
    }


    E_Int nUnique = (isOpen == 0) ? npts - 1 : npts; 
    E_Int start = 0;
    E_Int end = (isOpen == 0) ? nUnique : nUnique - 3; 

    // step 2 invalid for some cases, correction by taking step = 1
    if (step == 2)
    {
      step = stepCorrection (isOpen, nUnique);
    }

    for (E_Int k = 0; k < sweeps; k++)
    {

      for (E_Int j = 0; j < twoWays + 1; j++)
      {
        
        for (E_Int i = start; i < end; i = i + step)
        {
          
          if (isOpen == 0)
          {
            // Cycle, always in [0, nUnique-1]
            idx0 = (j == 0) ? i % nUnique : (nUnique - i % nUnique) % nUnique;
            idx1 = (j == 0) ? (i + 1) % nUnique : (nUnique - ((i + 1)   % nUnique)) % nUnique;
            idx2 = (j == 0) ? (i + 2) % nUnique : (nUnique - ((i + 2)   % nUnique)) % nUnique;
            idx3 = (j == 0) ? (i + 3) % nUnique : (nUnique - ((i + 3)   % nUnique)) % nUnique;
          }
          else 
          {
            idx0 = (j == 0) ? i  : (nUnique-1) - i;
            idx1 = (j == 0) ? (i + 1) : (nUnique-1) - (i + 1);
            idx2 = (j == 0) ? (i + 2)  : (nUnique-1) - (i + 2);
            idx3 = (j == 0) ? (i + 3)  : (nUnique-1) - (i + 3) ;
          }

          /* Get points i, i+1, i+2, i+3 */
          E_Float xi = x[idx0], yi = y[idx0], zi = z[idx0]; 
          E_Float xip1 = x[idx1], yip1 = y[idx1], zip1 = z[idx1]; 
          E_Float xip2 = x[idx2], yip2 = y[idx2], zip2 = z[idx2]; 
          E_Float xip3 = x[idx3], yip3 = y[idx3], zip3 = z[idx3];

          /* Compute deltas (i+3 - i), (i+2 - i) et (i+1 - i) */
          E_Float dv1x = xip1 - xi, dv1y = yip1 - yi, dv1z = zip1 - zi; /* xi+1 - xi */
          E_Float dv2x = xip2 - xi, dv2y = yip2 - yi, dv2z = zip2 - zi; /* xi+2 - xi */
          E_Float dv3x = xip3 - xi, dv3y = yip3 - yi, dv3z = zip3 - zi; /* xi+3 - xi */

          /* Set uNormal = unit normal to baseline (i+3;i) */
          E_Float normeDv3 = (dv3x * dv3x) + (dv3y * dv3y) + (dv3z * dv3z); /*  ||{xi+3-xi}**||² = ||xi+3-xi||² */
          
          if (normeDv3 < 1e-12) continue;
          E_Float divNorme = 1.0 / normeDv3 ;/*  1 / ||xi+3-xi||² */

          E_Float uNormalx = divNorme * (-dv3y), uNormaly = divNorme * (dv3x), uNormalz = 0;  /* {xi+3-xi}** / ||xi+3-xi||² */

          /* Compute signed area */
          E_Float aire = 0.5 * (dv3x * dv2y - dv3y * dv2x) + \
                        0.5 * (dv2x * dv1y - dv2y * dv1x) ;

          E_Float h = 1.5 * aire;

          x[idx1] = 2.0/3.0 * xi + 1.0 /3.0 * xip3 + h * uNormalx; y[idx1] = 2.0/3.0 * yi + 1.0 /3.0 * yip3 + h * uNormaly; z[idx1] = 2.0/3.0 * zi + 1.0 /3.0 * zip3 + h * uNormalz;
          x[idx2] = 1.0/3.0 * xi + 2.0 /3.0 * xip3 + h * uNormalx; y[idx2] = 1.0/3.0 * yi + 2.0 /3.0 * yip3 + h * uNormaly; z[idx2] = 1.0/3.0 * zi + 2.0 /3.0 * zip3 + h * uNormalz;
          
        }
      }
    
    }

    if (isOpen == 0)
    {
      x[nUnique] = x[0]; 
      y[nUnique] = y[0]; 
      z[nUnique] = z[0];
    }

    RELEASESHAREDS(tpl, fo); // free memory
  
  }
  else if (res == 2) // unstructured
  {  
    tpl = K_ARRAY::buildArray3(f->getNfld(), varString,
                               f->getSize(), cn->getSize(),
                               eltType, false, f->getApi());
    K_ARRAY::getFromArray3(tpl, fo, cno);

    E_Float* x = fo->begin(posx);
    E_Float* y = fo->begin(posy);
    E_Float* z = fo->begin(posz);

    for (E_Int v = 1; v <= f->getNfld(); v++) 
    {
      E_Float* fp1 = f->begin(v);
      E_Float* fp2 = fo->begin(v);
      for (E_Int i = 0; i < f->getSize(); i++) fp2[i] = fp1[i];
    }
    
    for (E_Int n = 1; n <= cno->getNfld(); n++) 
    {
      for (E_Int i = 0; i < cno->getSize(); i++) 
      {
        (*cno)(i,n) = (*cn)(i,n);
      }
    }

    E_Int npts = f->getSize(); 

    std::vector< std::vector<E_Int> > nodeAdj(npts);// vertex/elt adjacents
    K_CONNECT::connectEV2VE(*cn, nodeAdj);
  
    // Open or closed, find boundary nodes
    E_Int n1 = -1, n2 = -1;
    E_Int isOpen = 0;

    for (E_Int p = 0; p < npts; p++) 
    {
      if ((E_Int)nodeAdj[p].size() == 1) 
      {
        if (n1 == -1) n1 = p;
        else n2 = p;
      } 
    }

    if (n1 != -1 && n2 != -1)
    {
      E_Int elt = nodeAdj[n1][0]; // ID de l'élément
      E_Int nA = (*cn)(elt, 1) - 1;
      E_Int nB = (*cn)(elt, 2) - 1;
      E_Int n1Ngb = (nA == n1) ? nB : nA;

      isOpen = relativeDist(x[n1], y[n1], z[n1], x[n1Ngb], y[n1Ngb], z[n1Ngb], x[n2], y[n2], z[n2]);

      if (isOpen == 0)
      {
        // CASE 2: mesh is topologically OPEN, but geometrically CLOSED.
        printf("Info: consSmooth: closed geometry with double points: %d and %d.\n", n1, n2);
      }
      else // CASE 3: mesh is OPEN (extremities will stay fixed)

      {
        printf("Info: consSmooth: open geometry: fixed nodes %d and %d.\n", n1, n2);
      } 
    }
    else 
    {
      printf("Info: consSmooth: closed geometry.\n");
    }

    // Create a chained of ordered nodes
    std::vector<E_Int> bar;
    
    E_Int firstNode = (isOpen == 1) ? n1
                    : (n1 != -1) ? n1 // CASE 2: start from first double point
                    : ((*cn)(0, 1) - 1); // CASE 1: closed, start from any node 
    E_Int cur = firstNode, prev = -1;
    bar.push_back(cur);

    while (true)
    {
      E_Int next = -1;
      
      for (E_Int elt : nodeAdj[cur])  // elt = element index
      {
        E_Int nA = (*cn)(elt, 1) - 1;
        E_Int nB = (*cn)(elt, 2) - 1;
        E_Int neighbor = (nA == cur) ? nB : nA; // neighbour is the other node of element
        if (neighbor != prev) { next = neighbor; break; }
      }

      if (next == -1) break;                              // CASE 3: node with only one node reached = end of open chain
      if (isOpen == 0 && next == firstNode) break;        // CASE 1: return to start
      if (isOpen == 0 && n1 != -1 && next == n2) break;   // CASE 2: we dont keep double point
      
      bar.push_back(next);
      prev = cur;
      cur = next;
      if (isOpen == 1 && cur == n2) break; // we reach the end of opened curve
    }

    E_Int nUnique = (E_Int)bar.size();

    E_Int start = 0;
    E_Int end = (isOpen == 0) ? nUnique : nUnique - 3;

    // step 2 invalid for some cases, correction by taking step = 1
    if (step == 2)
    {
      step = stepCorrection(isOpen, nUnique);
    }
  
    for (E_Int k = 0; k < sweeps; k++)
    {
      for (E_Int j = 0; j < twoWays + 1; j++)
      {
        for (E_Int i = start; i < end; i += step)
        {
          E_Int c0, c1, c2, c3; // indexes in BAR
          if (isOpen == 0)
          {
            // Closed case
            c0 = (j == 0) ? i % nUnique : (nUnique - i % nUnique) % nUnique;
            c1 = (j == 0) ? (i + 1) % nUnique : (nUnique - ((i + 1)   % nUnique)) % nUnique;
            c2 = (j == 0) ? (i + 2) % nUnique : (nUnique - ((i + 2)   % nUnique)) % nUnique;
            c3 = (j == 0) ? (i + 3) % nUnique : (nUnique - ((i + 3)   % nUnique)) % nUnique;
          }
          else
          {
            // Open case
            c0 = (j == 0) ? i : (nUnique - 1) - i;
            c1 = (j == 0) ? (i + 1) : (nUnique - 1) - (i + 1);
            c2 = (j == 0) ? (i + 2) : (nUnique - 1) - (i + 2);
            c3 = (j == 0) ? (i + 3) : (nUnique - 1) - (i + 3);
          }

          // Nodes indices
          E_Int idx0 = bar[c0];
          E_Int idx1 = bar[c1];
          E_Int idx2 = bar[c2];
          E_Int idx3 = bar[c3];

          E_Float xi = x[idx0], yi = y[idx0], zi = z[idx0]; 
          E_Float xip1 = x[idx1], yip1 = y[idx1], zip1 = z[idx1]; 
          E_Float xip2 = x[idx2], yip2 = y[idx2], zip2 = z[idx2]; 
          E_Float xip3 = x[idx3], yip3 = y[idx3], zip3 = z[idx3];

          /* Compute deltas (i+3 - i), (i+2 - i) et (i+1 - i) */
          E_Float dv1x = xip1 - xi, dv1y = yip1 - yi, dv1z = zip1 - zi; /* xi+1 - xi */
          E_Float dv2x = xip2 - xi, dv2y = yip2 - yi, dv2z = zip2 - zi; /* xi+2 - xi */
          E_Float dv3x = xip3 - xi, dv3y = yip3 - yi, dv3z = zip3 - zi; /* xi+3 - xi */

          /* Set uNormal = unit normal to baseline (i+3;i) */
          E_Float normeDv3 = (dv3x * dv3x) + (dv3y * dv3y) + (dv3z * dv3z); /*  ||{xi+3-xi}**||² = ||xi+3-xi||² */
          
          if (normeDv3 < 1e-12) continue;
          E_Float divNorme = 1.0 / normeDv3 ;/*  1 / ||xi+3-xi||² */

          E_Float uNormalx = divNorme * (-dv3y), uNormaly = divNorme * (dv3x), uNormalz = 0;  /* {xi+3-xi}** / ||xi+3-xi||² */

          /* Compute signed area */
          E_Float aire = 0.5 * (dv3x * dv2y - dv3y * dv2x) + \
                        0.5 * (dv2x * dv1y - dv2y * dv1x) ;

          E_Float h = 1.5 * aire;

          x[idx1] = 2.0/3.0 * xi + 1.0 /3.0 * xip3 + h * uNormalx; y[idx1] = 2.0/3.0 * yi + 1.0 /3.0 * yip3 + h * uNormaly; z[idx1] = 2.0/3.0 * zi + 1.0 /3.0 * zip3 + h * uNormalz;
          x[idx2] = 1.0/3.0 * xi + 2.0 /3.0 * xip3 + h * uNormalx; y[idx2] = 1.0/3.0 * yi + 2.0 /3.0 * yip3 + h * uNormaly; z[idx2] = 1.0/3.0 * zi + 2.0 /3.0 * zip3 + h * uNormalz;
          
        }
      }
    }

    if (isOpen == 0 && n1 != -1 && n2 != -1)
    {
      x[n2] = x[n1];
      y[n2] = y[n1];
      z[n2] = z[n1];
    }

    RELEASESHAREDU(tpl, fo, cno); // free memory
  }

  RELEASESHAREDB(res, array, f, cn); 
  return tpl;

}



