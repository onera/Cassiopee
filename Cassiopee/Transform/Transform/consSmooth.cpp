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
# include <algorithm>
# include <random>
# include <array>
using namespace K_FLD;
using namespace K_SEARCH;

//=============================================================================
/* Local functions */
//=============================================================================

// compute relative distance between p1 and p2, ref taken from p1 and nbg1
inline bool relativeDist(E_Float p1x, E_Float p1y, E_Float p1z, E_Float nbg1x, E_Float nbg1y, E_Float nbg1z, E_Float p2x, E_Float p2y, E_Float p2z) 
{
  E_Float refL2 = (p1x - nbg1x)*(p1x - nbg1x) + (p1y - nbg1y)*(p1y - nbg1y) + (p1z - nbg1z)*(p1z - nbg1z);
  E_Float dist2 = (p1x - p2x)*(p1x - p2x) + (p1y - p2y)*(p1y - p2y) + (p1z - p2z)*(p1z - p2z);
  return (dist2 >= 1e-12 * refL2);
}

inline E_Int pbStep (E_Int isOpen, E_Int nUnique) 
{
  // Case of a CLOSED curve with an ODD number of points
  if ((isOpen == 0) && (nUnique % 2 != 0)) 
  {
    printf("Warning: consSmooth: step=2 will not give a perfect symetric smoothing for closed curve with an odd number of points.\n");
    return 1;
  }
  
  // Case of an OPEN curve with an EVEN number of points
  if ((isOpen == 1) && (nUnique % 2 == 0)) 
  {
    printf("Warning: consSmooth: step=2 will not give a perfect symetric smoothing for open curve with an even number of points.\n");
    return 1;
  }
  
  return 0;
}

// compute the diff between v1 and v2, return ej
inline void calDiff(E_Float v1x, E_Float v1y, E_Float v1z,
  E_Float v2x, E_Float v2y, E_Float v2z,
  E_Float& ejx, E_Float& ejy, E_Float& ejz)
{
  ejx = v2x - v1x;
  ejy = v2y - v1y;
  ejz = v2z - v1z;
}

// compute the cross product of e1 and e2, return v
inline void calProdVec(E_Float e1x, E_Float e1y, E_Float e1z,
  E_Float e2x, E_Float e2y, E_Float e2z,
  E_Float& vx, E_Float& vy, E_Float& vz)
{
  vx = e1y*e2z - e1z*e2y;
  vy = e1z*e2x - e1x*e2z;
  vz = e1x*e2y - e1y*e2x;
}

inline void calArea(E_Float* x, E_Float* y, E_Float* z,
  const std::vector<E_Int>& ring, E_Int nNbrs, E_Int p,
  E_Float& Ax, E_Float& Ay, E_Float& Az, bool& inverted, 
  E_Float refx = 0.0, E_Float refy = 0.0, E_Float refz = 0.0)
{
  Ax = 0.0; Ay = 0.0; Az = 0.0;
  inverted = false;

  E_Float ejx, ejy, ejz, epx, epy, epz;
  E_Int jNext;
  // calculate the total area (reference orientation)
  std::vector<E_Float> axs(nNbrs), ays(nNbrs), azs(nNbrs);
  for (E_Int j = 0; j < nNbrs; j++)
  {
    jNext = (j+1) % nNbrs;
    ejx = x[ring[j]]-x[p], ejy = y[ring[j]]-y[p], ejz = z[ring[j]]-z[p];
    epx = x[ring[jNext]]-x[p], epy = y[ring[jNext]]-y[p], epz = z[ring[jNext]]-z[p];
    calProdVec(ejx, ejy, ejz, epx, epy, epz, axs[j], ays[j], azs[j]);
    Ax += axs[j]; Ay += ays[j]; Az += azs[j];
  }

  // Reference either giver or uses global orientation
  E_Float bx = (refx != 0.0 || refy != 0.0 || refz != 0.0) ? refx : Ax;
  E_Float by = (refx != 0.0 || refy != 0.0 || refz != 0.0) ? refy : Ay;
  E_Float bz = (refx != 0.0 || refy != 0.0 || refz != 0.0) ? refz : Az;

  // Global inversion: total area opposite the reference
  if (Ax*bx + Ay*by + Az*bz < 0.0) 
  { 
    inverted = true; 
    return; 
  }

  // Local inversion: an individual triangle opposite the reference
  for (E_Int j = 0; j < nNbrs; j++)
  {
    if (axs[j]*bx + ays[j]*by + azs[j]*bz < 0.0)
    { 
      inverted = true; 
      break; 
    }
  }
}

// find neighbours of p (starts 0) starting with pStart and in the direct sense
bool buildRingAroundP(E_Int p, E_Int pStart, const std::vector<E_Int>& eltsOfP, FldArrayI* cn, std::vector<E_Int>& ring)
{
  E_Int nNbrs = eltsOfP.size();
  std::vector<std::array<E_Int,3>> pNbrs(nNbrs);

  // Extract the 2 nodes != p for each triangle.
  E_Int k = 0;
  for (E_Int elt : eltsOfP)
  {
    E_Int n0 = (*cn)(elt, 1) - 1;
    E_Int n1 = (*cn)(elt, 2) - 1;
    E_Int n2 = (*cn)(elt, 3) - 1;
    E_Int idx = 0;
    if (n0 != p) pNbrs[k][idx++] = n0;
    if (n1 != p) pNbrs[k][idx++] = n1;
    if (n2 != p) pNbrs[k][idx++] = n2;
    k++;
  }

  // Starting triangle: CCW sequence (p → pStart → cur)
  E_Int startK = -1, cur = -1;
  for (E_Int j = 0; j < nNbrs; j++)
  {
    E_Int elt = eltsOfP[j];
    E_Int n0 = (*cn)(elt, 1) - 1;
    E_Int n1 = (*cn)(elt, 2) - 1;
    E_Int n2 = (*cn)(elt, 3) - 1;
    if      (n0 == p && n1 == pStart) { cur = n2; startK = j; break; }
    else if (n1 == p && n2 == pStart) { cur = n0; startK = j; break; }
    else if (n2 == p && n0 == pStart) { cur = n1; startK = j; break; }
  }
  if (startK == -1) return false;

  ring.push_back(pStart);
  ring.push_back(cur);

  pNbrs[startK][0] = pNbrs[startK][1] = -1;

  // Chaining of the remaining nodes
  for (E_Int i = 2; i < nNbrs; i++)
  {
    bool found = false;
    E_Int pPrev = ring[i-1];
    for (E_Int j = 0; j < nNbrs; j++)
    {
      if (pNbrs[j][0] == -1) continue;
      E_Int pA = pNbrs[j][0];
      E_Int pB = pNbrs[j][1];
      if (pA == pPrev)
      {
        ring.push_back(pB);
        pNbrs[j][0] = pNbrs[j][1] = -1;
        found = true; break;
      }
      if (pB == pPrev)
      {
        ring.push_back(pA);
        pNbrs[j][0] = pNbrs[j][1] = -1;
        found = true; break;
      }
    }
    if (!found) return false; // boundary node
  }
          

  if (ring.size() != static_cast<std::size_t>(nNbrs)) return false;

  // Check: ensure the last node properly closes the loop at pStart.
  E_Int last = ring.back();
  for (E_Int elt : eltsOfP)
  {
    E_Int n0 = (*cn)(elt, 1) - 1;
    E_Int n1 = (*cn)(elt, 2) - 1;
    E_Int n2 = (*cn)(elt, 3) - 1;
    if ((n0==p && n1==last && n2==pStart) ||
        (n1==p && n2==last && n0==pStart) ||
        (n2==p && n0==last && n1==pStart))
      return true;
  }
  return false; // unclosed loop = boundary node

}

// ============================================================================
/* Conservative smoothing (Kuprat) */
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
  E_Float omega;
  if (!PYPARSETUPLE_(args, O_ III_ R_, &array, &sweeps,  &twoWays, &step, &omega))
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
  PyObject* tpl = NULL;
  FldArrayF* fo;
  FldArrayI* cno = NULL;
  
  // Get dimension A VERIFIER PLUS INTENSEMENT
  E_Int dim = 3;
  if (res == 1) // Structured
  {
    E_Int im1 = im - 1;
    E_Int jm1 = jm - 1;
    E_Int km1 = km - 1;
    if (im1 * jm1 * km1 == 0)
    {
      if ((im1 && (jm1 || km1)) || (jm1 && km1)) dim = 2;
      else dim = 1;
    }
  }
  else if (res == 2 && eltType != nullptr) // Unstructured
  {
    dim = K_CONNECT::getDimME(eltType);
  }
  
  if (dim == 1) // 1d curves
  {
    //printf ("consSmooth: on curves.\n");
    if (res == 1) // structured
    {
      tpl = K_ARRAY::buildArray3(*f, varString, im, jm, km, f->getApi());
      K_ARRAY::getFromArray3(tpl, fo);

      E_Float* x = fo->begin(posx);
      E_Float* y = fo->begin(posy);
      E_Float* z = fo->begin(posz);

      E_Int npts = im;      
      E_Int idx0, idx1, idx2, idx3;

      // open or closed ?
      E_Int isOpen = relativeDist(x[0], y[0], z[0], x[1], y[1], z[1], x[npts-1], y[npts-1], z[npts-1]);

      if (isOpen)
      {
        printf("Info: consSmooth: open geometry: fixed nodes " SF_D_ " and " SF_D_ ".\n", 0, npts-1);
      }
      else
      {
        printf("Info: consSmooth: closed geometry with double points.\n");
      }

      E_Int nUnique = (isOpen == 0) ? npts - 1 : npts; 
      E_Int start = 0;
      E_Int end = (isOpen == 0) ? nUnique : nUnique - 3; 


      E_Int isPbStep = 0;
      
      if (step == 2)
      {
        isPbStep = pbStep (isOpen, nUnique); 
      }

      E_Int currStart = start;
      E_Int currEnd = end;
      
      for (E_Int k = 0; k < sweeps; k++)
      {

        if (isPbStep == 1 && isOpen == 0)
        {
          currStart = start + k;
          currEnd = end + k;
        }

        for (E_Int j = 0; j < twoWays + 1; j++)
        {

          for (E_Int i = currStart; i < currEnd; i = i + step)
          {
            
            if (isOpen == 0)
            {
              // Cycle, always in [0, nUnique-1]
              if (j == 0) 
              {
                idx0 = i % nUnique;
                idx1 = (i + 1) % nUnique;
                idx2 = (i + 2) % nUnique;
                idx3 = (i + 3) % nUnique;
              } 
              else 
              {
                E_Int offset = 0;
                if (isPbStep == 1)
                {
                  offset = (j % 2 != 0) ? 1 : 0;
                }
                idx0 = (nUnique - (i % nUnique) + offset) % nUnique;
                idx1 = (nUnique - ((i + 1) % nUnique) + offset) % nUnique;
                idx2 = (nUnique - ((i + 2) % nUnique) + offset) % nUnique;
                idx3 = (nUnique - ((i + 3) % nUnique) + offset) % nUnique;
              }
            }
            else 
            {
              if (j == 0) 
              {
                idx0 = i;
                idx1 = (i + 1);
                idx2 = (i + 2);
                idx3 = (i + 3);
              } 
              else 
              {
                E_Int offset = 0;
                if (isPbStep == 1)
                {
                  offset = (j % 2 != 0) ? 1 : 0;
                }
                
                idx0 = (nUnique-1) - (i + offset);
                idx1 = (nUnique-1) - (i + 1 + offset);
                idx2 = (nUnique-1) - (i + 2 + offset);
                idx3 = (nUnique-1) - (i + 3 + offset);

                if (idx0 < 4)
                {
                  continue;
                }
              }
            }

            // printf("j == %d: on regarde les points : %d, %d, %d, %d ET on change les points : %d and %d.\n", j, idx0, idx1, idx2, idx3, idx1, idx2);

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

            /* Determines new xi+1 and xi+2 */
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

      // to check if open cal relative dist of closure  (get ngb)
      if (n1 != -1 && n2 != -1)
      {
        E_Int elt = nodeAdj[n1][0]; // ID of element
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

      E_Int isPbStep = 0;
      E_Int currStart = start;
      E_Int currEnd = end;
      
      if (step == 2)
      {
        isPbStep = pbStep (isOpen, nUnique); 
      }
    
      for (E_Int k = 0; k < sweeps; k++)
      {
        if (isPbStep == 1 && isOpen == 0)
        {
          currStart = start + k;
          currEnd = end + k;
        }
        for (E_Int j = 0; j < twoWays + 1; j++)
        {
          for (E_Int i = currStart; i < currEnd; i += step)
          {
            E_Int c0, c1, c2, c3; // indexes in BAR
            if (isOpen == 0)
            {
              // Closed case
              // Cycle, always in [0, nUnique-1]
              if (j == 0) 
              {
                c0 = i % nUnique;
                c1 = (i + 1) % nUnique;
                c2 = (i + 2) % nUnique;
                c3 = (i + 3) % nUnique;
              } 
              else 
              {
                E_Int offset = 0;
                if (isPbStep == 1)
                {
                  offset = (j % 2 != 0) ? 1 : 0;
                }
                c0 = (nUnique - (i % nUnique) + offset) % nUnique;
                c1 = (nUnique - ((i + 1) % nUnique) + offset) % nUnique;
                c2 = (nUnique - ((i + 2) % nUnique) + offset) % nUnique;
                c3 = (nUnique - ((i + 3) % nUnique) + offset) % nUnique;
              }
            }
            else
            {
              // Open case
              if (j == 0) 
              {
                c0 = i;
                c1 = (i + 1);
                c2 = (i + 2);
                c3 = (i + 3);
              } 
              else 
              {
                E_Int offset = 0;
                if (isPbStep == 1)
                {
                  offset = (j % 2 != 0) ? 1 : 0;
                }
                
                c0 = (nUnique-1) - (i + offset);
                c1 = (nUnique-1) - (i + 1 + offset);
                c2 = (nUnique-1) - (i + 2 + offset);
                c3 = (nUnique-1) - (i + 3 + offset);

                if (c0 < 4)
                {
                  continue;
                }
              }
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

            /* Determines new xi+1 and xi+2 */
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
  } 
  else if (dim == 2) // TRI surface
  {
    //printf ("consSmooth: on surfaces.\n");
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

    E_Int nelts = cn->getSize(); // number of rows = number of elements (triangles)
    E_Int npts = f->getSize(); // total number of nodes
    E_Int ns = cn->getNfld(); // = 3 for TRI

    std::vector< std::vector<E_Int> > cVE(npts);// vertex/elt adjacents
    K_CONNECT::connectEV2VE(*cn, cVE);

    // retrieval of internal edges (count == 2)
    std::vector<std::pair<E_Int, E_Int>> activeEdges;
    E_Float maxLen2 = 0.0;
    for (E_Int t = 0; t < nelts; t++)
    {
      for (E_Int edgeT = 0; edgeT < ns; edgeT++)
      {
        E_Int p1 = (*cn)(t, edgeT + 1) - 1;
        E_Int p2 = (*cn)(t, (edgeT + 1) % ns + 1) - 1;
        
        if (p1 > p2) continue; // Uniqueness

        // reference length (initial maximum edge)
        E_Float d2 = (x[p1]-x[p2])*(x[p1]-x[p2]) + (y[p1]-y[p2])*(y[p1]-y[p2]) + (z[p1]-z[p2])*(z[p1]-z[p2]);
        if (d2 > maxLen2) maxLen2 = d2;

        // Boundary check
        E_Int count = 0;
        for (E_Int a : cVE[p1])
          for (E_Int b : cVE[p2])
            if (a == b) count++;
        
        if (count == 2) 
        {
          activeEdges.push_back({p1, p2});
          //printf("active p1=%d p2=%d\n", p1, p2);
        }
      }
    }

    // initial normals
    std::vector<E_Float> refNx(npts, 0.0), refNy(npts, 0.0), refNz(npts, 0.0);
    for (E_Int p = 0; p < npts; p++)
    {
      if (cVE[p].empty()) continue;
      // Sum of the areas of all triangles adjacent to p
      for (E_Int t : cVE[p])
      {
        E_Int n0 = (*cn)(t, 1) - 1;
        E_Int n1 = (*cn)(t, 2) - 1;
        E_Int n2 = (*cn)(t, 3) - 1;
        E_Float ex = x[n1]-x[n0], ey = y[n1]-y[n0], ez = z[n1]-z[n0];
        E_Float fx = x[n2]-x[n0], fy = y[n2]-y[n0], fz = z[n2]-z[n0];
        E_Float ax, ay, az;
        calProdVec(ex, ey, ez, fx, fy, fz, ax, ay, az);
        refNx[p] += ax; refNy[p] += ay; refNz[p] += az;
      }
      // Normalisation
      E_Float norm = sqrt(refNx[p]*refNx[p] + refNy[p]*refNy[p] + refNz[p]*refNz[p]);
      if (norm > 1e-12) { refNx[p] /= norm; refNy[p] /= norm; refNz[p] /= norm; }
    }

    //E_Float lRef = (maxLen2 > 1e-28) ? sqrt(maxLen2) : 1.0;
    E_Int incoherence = 0;
    E_Float normAMin = 1e-28;

    for (E_Int sweep = 0; sweep < sweeps; sweep++)
    {
      std::mt19937 rng(sweep);
      std::shuffle(activeEdges.begin(), activeEdges.end(), rng);

      // each node processed at most once per pass
      E_Int nActiveEdges = (E_Int)activeEdges.size();
      // printf("There are %d edges to modify \n", nActiveEdges);

      for (E_Int edgeIdx = 0; edgeIdx < nActiveEdges; edgeIdx++)
      {    
        E_Int p1 = activeEdges[edgeIdx].first; // x1 de Kuprat
        E_Int p2 = activeEdges[edgeIdx].second; // x2 de Kuprat

        //printf("traitement edge %d: " SF_D_ " -> " SF_D_ "\n", edgeIdx, p1, p2);
        // printf("[triangle %d/%d] [%d/%d] p1=%d, p2=%d \n", t, nelts, edgeT, ns, p1, p2);

        // Ring around p1
        std::vector<E_Int> ring1;
        E_Int closed1 = buildRingAroundP(p1, p2, cVE[p1], cn, ring1);
        //printf("ring " SF_D_ " = ", p1);
        //for (size_t i = 0; i < ring1.size(); i++) printf("%d ", ring1[i]);
        //printf("\n");
        //printf("closed1=%d\n", closed1);
          
        //  Ring around p2
        std::vector<E_Int> ring2;
        E_Int closed2 = buildRingAroundP(p2, p1, cVE[p2], cn, ring2);
        //printf("ring %d = ", p2);
        //for (size_t i = 0; i < ring2.size(); i++) printf("%d ", ring2[i]);
        //printf("\n");
        //printf("closed2=%d\n", closed2);
        E_Int n1 = ring1.size();
        E_Int n2 = ring2.size();

        // check
        if (ring1[0] != p2) printf("PB1\n");
        if (ring2[0] != p1) printf("PB1\n");
        if (ring1[1] != ring2[n2-1]) printf("PB1\n");
        if (ring1[n1-1] != ring2[1]) printf("PB1\n");
          
        if (!closed1 || !closed2) 
        {
          continue;  // untreatable edge → skip
        }

        E_Int nNbrs1 = (E_Int)ring1.size();
        E_Int nNbrs2 = (E_Int)ring2.size();

        E_Float Ax1, Ay1, Az1;
        bool inverted1;
        calArea(x, y, z, ring1, nNbrs1, p1, Ax1, Ay1, Az1, inverted1, refNx[p1], refNy[p1], refNz[p1]);
        //printf("area1=%g %g %g\n", Ax1, Ay1, Az1);
        E_Float Ax2, Ay2, Az2;
        bool inverted2;
        calArea(x, y, z, ring2, nNbrs2, p2, Ax2, Ay2, Az2, inverted2, refNx[p2], refNy[p2], refNz[p2]);
        //printf("area2=%g %g %g\n", Ax2, Ay2, Az2);
          
        if (inverted1 || inverted2) incoherence += 1;
        
        E_Float vx = x[ring2[nNbrs2-1]] - x[ring2[1]];
        E_Float vy = y[ring2[nNbrs2-1]] - y[ring2[1]];
        E_Float vz = z[ring2[nNbrs2-1]] - z[ring2[1]];

        // check
        if (vx != x[ring1[1]] - x[ring1[n1-1]]) printf("PB2\n");
        if (vy != y[ring1[1]] - y[ring1[n1-1]]) printf("PB2\n");
        if (vz != z[ring1[1]] - z[ring1[n1-1]]) printf("PB2\n");

        E_Float sum1x = 0.0, sum1y = 0.0, sum1z = 0.0;
        for (E_Int s = 1; s < nNbrs1; s++)
        {
          sum1x += x[ring1[s]];
          sum1y += y[ring1[s]];
          sum1z += z[ring1[s]];
        }

        E_Float sum2x = 0.0, sum2y = 0.0, sum2z = 0.0;
        for (E_Int s = 1; s < nNbrs2; s++)
        {
          sum2x += x[ring2[s]];
          sum2y += y[ring2[s]];
          sum2z += z[ring2[s]];
        }

        // Equation (16) : x1^s
        E_Float denom = 1.0 / (E_Float)(nNbrs1 * nNbrs2 - 1);
        E_Float p1sx = (sum2x + (E_Float)nNbrs2 * sum1x) * denom;
        E_Float p1sy = (sum2y + (E_Float)nNbrs2 * sum1y) * denom;
        E_Float p1sz = (sum2z + (E_Float)nNbrs2 * sum1z) * denom;

        // Equation (15) : x2^s
        E_Float denom2 = 1.0 / (E_Float)nNbrs2;
        E_Float p2sx = (p1sx + sum2x) * denom2;
        E_Float p2sy = (p1sy + sum2y) * denom2;
        E_Float p2sz = (p1sz + sum2z) * denom2;
        
        E_Float dp1sx = omega * (p1sx - x[p1]), dp1sy = omega * (p1sy - y[p1]), dp1sz = omega * (p1sz - z[p1]); 
        E_Float dp2sx = omega * (p2sx - x[p2]), dp2sy = omega * (p2sy - y[p2]), dp2sz = omega * (p2sz - z[p2]); 

        // A = A1 + A2 + v × (dx1s - dx2s)
        E_Float diffx = dp1sx - dp2sx, diffy = dp1sy - dp2sy, diffz = dp1sz - dp2sz;
        E_Float crossx, crossy, crossz;
        calProdVec(vx, vy, vz, diffx, diffy, diffz, crossx, crossy, crossz);

        E_Float Ax = Ax1 + Ax2 + crossx;
        E_Float Ay = Ay1 + Ay2 + crossy;
        E_Float Az = Az1 + Az2 + crossz;
        E_Float normA = sqrt(Ax*Ax + Ay*Ay + Az*Az);

        if (normA > normAMin)
        {
          E_Float divnormA = 1.0 / normA;
          E_Float nx = Ax*divnormA, ny = Ay*divnormA, nz = Az*divnormA;
          //printf("nx=%g %g %g\n", nx, ny, nz);
              
          // dp1s · A1
          E_Float dot1 = dp1sx*Ax1 + dp1sy*Ay1 + dp1sz*Az1;

          // dp2s · A2
          E_Float dot2 = dp2sx*Ax2 + dp2sy*Ay2 + dp2sz*Az2;
              
          // dp2s · (v × dp1s)
          E_Float vxdp1x, vxdp1y, vxdp1z;
          calProdVec(vx, vy, vz, dp1sx, dp1sy, dp1sz, vxdp1x, vxdp1y, vxdp1z);
          E_Float dot3 = dp2sx*vxdp1x + dp2sy*vxdp1y + dp2sz*vxdp1z;

          E_Float h = -(dot1 + dot2 + dot3) / normA;
              
          //printf("h=%g\n", h);
          //printf("dp1s=%g %g %g\n", dp1sx, dp1sy, dp1sz);
          //printf("dp2s=%g %g %g\n", dp2sx, dp2sy, dp2sz);
          // final update
          x[p1] += dp1sx + h*nx;  y[p1] += dp1sy + h*ny;  z[p1] += dp1sz + h*nz;
          x[p2] += dp2sx + h*nx;  y[p2] += dp2sy + h*ny;  z[p2] += dp2sz + h*nz;

          // Check that the A1 A2 have not changed sign
          //E_Float Ax1n, Ay1n, Az1n; bool inv1n = false;
          //calArea(x, y, z, ring1, nNbrs1, p1, Ax1n, Ay1n, Az1n, inv1n, refNx[p1], refNy[p1], refNz[p1]);
          //E_Float Ax2n, Ay2n, Az2n; bool inv2n = false;
          //calArea(x, y, z, ring2, nNbrs2, p2, Ax2n, Ay2n, Az2n, inv2n, refNx[p2], refNy[p2], refNz[p2]);
          //if (inv1n || inv2n)
          //  printf("Warning: consSmooth: [%d Sweep][%d/%d] Inversion %d %d try omega = %f \n",sweep, countE, nActiveEdges,inv1n,inv2n, omega);
        }
      }
      if (incoherence > 0)
        printf("Warning: consSmooth: " SF_D_ " normal inversions\n", incoherence);
    }
    RELEASESHAREDU(tpl, fo, cno);
  }

  RELEASESHAREDB(res, array, f, cn); 
  return tpl;

}
