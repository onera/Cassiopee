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

// convertit un maillage quad en maillage tri en coupant suivant l'angle max

# include "converter.h"
# include "kcore.h"
# include <string.h>
# include <stdio.h>

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Conversion du maillage quad en maillage tri.
   Si l'input contient un champ indic, decoupe suivant
   ce champ, sinon
   on coupe suivant l'angle alpha maximum. */
//=============================================================================
PyObject* K_CONVERTER::convertQuad2Tri(PyObject* self, PyObject* args)
{
  PyObject* pArray;
  if (!PYPARSETUPLE_(args, O_, &pArray)) return NULL;
 
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(pArray, varString,
                                     f, ni, nj, nk, cn, eltType);

  // Test non structure ?
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(pArray,f);
    PyErr_SetString(PyExc_TypeError,
                    "convertQuad2tri: input array must be unstructured.");
    return NULL;
  }

  // Test QUAD
  if (strcmp(eltType, "QUAD") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertQuad2Tri: unstructured array must be QUAD.");
    RELEASESHAREDU(pArray, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertQuad2Tri: coord must be present in array.");
    RELEASESHAREDU(pArray, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  E_Int posi = K_ARRAY::isNamePresent("indic", varString); posi++;

  // Connectivite
  E_Int ne = cn->getSize();

  // TRI de sortie
  E_Int api = f->getApi();
  E_Int nt = 2*ne;
  FldArrayI ct(nt, 3);
  E_Int ntr = 0; // nbre reel de triangles
  E_Int* ct1 = ct.begin(1);
  E_Int* ct2 = ct.begin(2);
  E_Int* ct3 = ct.begin(3);
  
  E_Int ind1, ind2, ind3, ind4;
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  FldArrayI& cnp = *cn;
  E_Int* cnp1 = cnp.begin(1);
  E_Int* cnp2 = cnp.begin(2);
  E_Int* cnp3 = cnp.begin(3);
  E_Int* cnp4 = cnp.begin(4);

  if (posi == 0) // decoupage sans indic
  {
    E_Float alpha1, alpha2, ndir;
    E_Float ptA[3], ptB[3], ptC[3], dir[3];

    for (E_Int i = 0; i < ne; i++)
    {
      ind1 = cnp1[i]-1; ind2 = cnp2[i]-1; ind3 = cnp3[i]-1; ind4 = cnp4[i]-1;
      // Angles 2D - 124
      ptA[0] = x[ind1]; ptA[1] = y[ind1]; ptA[2] = z[ind1];
      ptB[0] = x[ind2]; ptB[1] = y[ind2]; ptB[2] = z[ind2];
      ptC[0] = x[ind4]; ptC[1] = y[ind4]; ptC[2] = z[ind4];
      dir[0] = (ptB[1]-ptA[1])*(ptC[2]-ptA[2])-(ptB[2]-ptA[2])*(ptC[1]-ptA[1]);
      dir[1] = (ptB[2]-ptA[2])*(ptC[0]-ptA[0])-(ptB[0]-ptA[0])*(ptC[2]-ptA[2]);
      dir[2] = (ptB[0]-ptA[0])*(ptC[1]-ptA[1])-(ptB[1]-ptA[1])*(ptC[0]-ptA[0]);
      ndir = sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
      ndir = 1./K_FUNC::E_max(ndir, 1.e-10);
      dir[0] = dir[0]*ndir; dir[1] = dir[1]*ndir; dir[2] = dir[2]*ndir;
      alpha1 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA, ptB, ptA, ptC, dir);

      ptA[0] = x[ind2]; ptA[1] = y[ind2]; ptA[2] = z[ind2];
      ptB[0] = x[ind3]; ptB[1] = y[ind3]; ptB[2] = z[ind3];
      ptC[0] = x[ind1]; ptC[1] = y[ind1]; ptC[2] = z[ind1];
      dir[0] = (ptB[1]-ptA[1])*(ptC[2]-ptA[2])-(ptB[2]-ptA[2])*(ptC[1]-ptA[1]);
      dir[1] = (ptB[2]-ptA[2])*(ptC[0]-ptA[0])-(ptB[0]-ptA[0])*(ptC[2]-ptA[2]);
      dir[2] = (ptB[0]-ptA[0])*(ptC[1]-ptA[1])-(ptB[1]-ptA[1])*(ptC[0]-ptA[0]);
      ndir = sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
      ndir = 1./K_FUNC::E_max(ndir, 1.e-10);
      dir[0] = dir[0]*ndir; dir[1] = dir[1]*ndir; dir[2] = dir[2]*ndir;
      alpha2 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA, ptB, ptA, ptC, dir);

      //printf("alpha= %f %f\n", alpha1, alpha2);

      if (alpha1 > alpha2) // coupe le plus obtu
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind3+1;
        ct1[2*i+1] = ind1+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;      
      }
      else // decoupage suivant 2
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind4+1;
        ct1[2*i+1] = ind2+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;      
      }
      ntr += 2;
    }
  }
  else if (1 == 0)// decoupage suivant indic (old)
  {
    E_Float* indic = f->begin(posi);
    E_Int t1, t2, t3, t4;
    for (E_Int i = 0; i < ne; i++)
    {
      ind1 = cnp1[i]-1; ind2 = cnp2[i]-1; ind3 = cnp3[i]-1; ind4 = cnp4[i]-1;
      t1 = floor(indic[ind1]+0.5); t2 = floor(indic[ind2]+0.5); 
      t3 = floor(indic[ind3]+0.5); t4 = floor(indic[ind4]+0.5);
      if (t1 == 1 && t2 == 0 && t3 == 1 && t4 == 0) // 2
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind3+1;
        ct1[2*i+1] = ind1+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;  
      }
      else if (t1 == 0 && t2 == 1 && t3 == 0 && t4 == 1) // 2
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind4+1;
        ct1[2*i+1] = ind2+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;
      }
      else if (t1 == 1 && t2 == 0 && t3 == 1 && t4 == 1) // 3
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind3+1;
        ct1[2*i+1] = ind1+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;
      }
      else if (t1 == 1 && t2 == 1 && t3 == 1 && t4 == 0) // 3
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind3+1;
        ct1[2*i+1] = ind1+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;
      }
      else if (t1 == 0 && t2 == 1 && t3 == 1 && t4 == 1) // 3
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind4+1;
        ct1[2*i+1] = ind2+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;
      }
      else if (t1 == 1 && t2 == 1 && t3 == 0 && t4 == 1) // 3
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind4+1;
        ct1[2*i+1] = ind2+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;
      }
      else
      {
        ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind3+1;
        ct1[2*i+1] = ind1+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;      
      }  
    }
    ntr += 2;
  }
  else // decoupage suivant indic (2) + angle non obtu
  {
    //E_Float alpha1, alpha2, alpha3, alpha4;
    E_Float ndir1, ndir2, ndir3, ndir4, ndirl;
    E_Float ptA1[3], ptB1[3], ptC1[3], dir1[3];
    E_Float ptA2[3], ptB2[3], ptC2[3], dir2[3];
    E_Float ptA3[3], ptB3[3], ptC3[3], dir3[3];
    E_Float ptA4[3], ptB4[3], ptC4[3], dir4[3];
    
    E_Float* indic = f->begin(posi);
    E_Int t1, t2, t3, t4;
    bool deg1, deg2; E_Float inverse1, inverse2;

    // Premiere passe : elimination des edges contractes
    for (E_Int i = 0; i < ne; i++)
    {
      ind1 = cnp1[i]-1; ind2 = cnp2[i]-1; ind3 = cnp3[i]-1; ind4 = cnp4[i]-1;

      // Verification des aretes (non nulles)
      // si arete nulle -> triangle direct
      ptA1[0] = x[ind2]-x[ind1]; 
      ptA1[1] = y[ind2]-y[ind1];
      ptA1[2] = z[ind2]-z[ind1];
      ptA2[0] = x[ind3]-x[ind2];
      ptA2[1] = y[ind3]-y[ind2];
      ptA2[2] = z[ind3]-z[ind2];
      ptA3[0] = x[ind4]-x[ind3]; 
      ptA3[1] = y[ind4]-y[ind3];
      ptA3[2] = z[ind4]-z[ind3];
      ptA4[0] = x[ind1]-x[ind4];
      ptA4[1] = y[ind1]-y[ind4];
      ptA4[2] = z[ind1]-z[ind4];
      ndir1 = sqrt(ptA1[0]*ptA1[0]+ptA1[1]*ptA1[1]+ptA1[2]*ptA1[2]);
      ndir2 = sqrt(ptA2[0]*ptA2[0]+ptA2[1]*ptA2[1]+ptA2[2]*ptA2[2]);
      ndir3 = sqrt(ptA3[0]*ptA3[0]+ptA3[1]*ptA3[1]+ptA3[2]*ptA3[2]);
      ndir4 = sqrt(ptA4[0]*ptA4[0]+ptA4[1]*ptA4[1]+ptA4[2]*ptA4[2]);

      if (ndir1 < 1.e-12 && ndir2 < 1.e-12 && ndir3 < 1.e-12 && ndir4 < 1.e-12) goto end0;

      if (ndir1 < 1.e-12)
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind3+1; ct3[ntr] = ind4+1;
        ptA4[0] = x[ind4]; ptA4[1] = y[ind4]; ptA4[2] = z[ind4];
        ptB4[0] = x[ind1]; ptB4[1] = y[ind1]; ptB4[2] = z[ind1];
        ptC4[0] = x[ind3]; ptC4[1] = y[ind3]; ptC4[2] = z[ind3];
        dir4[0] = (ptB4[1]-ptA4[1])*(ptC4[2]-ptA4[2])-(ptB4[2]-ptA4[2])*(ptC4[1]-ptA4[1]);
        dir4[1] = (ptB4[2]-ptA4[2])*(ptC4[0]-ptA4[0])-(ptB4[0]-ptA4[0])*(ptC4[2]-ptA4[2]);
        dir4[2] = (ptB4[0]-ptA4[0])*(ptC4[1]-ptA4[1])-(ptB4[1]-ptA4[1])*(ptC4[0]-ptA4[0]);
        ndirl = sqrt(dir4[0]*dir4[0]+dir4[1]*dir4[1]+dir4[2]*dir4[2]);
        if (ndirl < 1.e-11) printf("convertTri2Quad: " SF_D_ ": quad=" SF_D_ " is tri degen with edge contracted\n", ntr, i);
        ntr++;
        goto end0; 
      }
      else if (ndir2 < 1.e-12) 
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind4+1; 
        ptA1[0] = x[ind1]; ptA1[1] = y[ind1]; ptA1[2] = z[ind1];
        ptB1[0] = x[ind2]; ptB1[1] = y[ind2]; ptB1[2] = z[ind2];
        ptC1[0] = x[ind4]; ptC1[1] = y[ind4]; ptC1[2] = z[ind4];
        dir1[0] = (ptB1[1]-ptA1[1])*(ptC1[2]-ptA1[2])-(ptB1[2]-ptA1[2])*(ptC1[1]-ptA1[1]);
        dir1[1] = (ptB1[2]-ptA1[2])*(ptC1[0]-ptA1[0])-(ptB1[0]-ptA1[0])*(ptC1[2]-ptA1[2]);
        dir1[2] = (ptB1[0]-ptA1[0])*(ptC1[1]-ptA1[1])-(ptB1[1]-ptA1[1])*(ptC1[0]-ptA1[0]);
        ndirl = sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
        if (ndirl < 1.e-11) printf("convertTri2Quad: " SF_D_ ": quad=" SF_D_ " is 2-tri degen with edge contracted\n", ntr, i);
        ntr++;
        goto end0;
      }
      else if (ndir3 < 1.e-12) 
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind3+1;
        ptA2[0] = x[ind2]; ptA2[1] = y[ind2]; ptA2[2] = z[ind2];
        ptB2[0] = x[ind3]; ptB2[1] = y[ind3]; ptB2[2] = z[ind3];
        ptC2[0] = x[ind1]; ptC2[1] = y[ind1]; ptC2[2] = z[ind1];
        dir2[0] = (ptB2[1]-ptA2[1])*(ptC2[2]-ptA2[2])-(ptB2[2]-ptA2[2])*(ptC2[1]-ptA2[1]);
        dir2[1] = (ptB2[2]-ptA2[2])*(ptC2[0]-ptA2[0])-(ptB2[0]-ptA2[0])*(ptC2[2]-ptA2[2]);
        dir2[2] = (ptB2[0]-ptA2[0])*(ptC2[1]-ptA2[1])-(ptB2[1]-ptA2[1])*(ptC2[0]-ptA2[0]);
        ndirl = sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2]);
        if (ndirl < 1.e-11) printf("convertTri2Quad:" SF_D_ ": quad=" SF_D_ " is 3-tri degen with edge contracted\n", ntr, i);
        ntr++;
        goto end0; 
      }
      else if (ndir4 < 1.e-12) 
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind3+1;
        ptA2[0] = x[ind2]; ptA2[1] = y[ind2]; ptA2[2] = z[ind2];
        ptB2[0] = x[ind3]; ptB2[1] = y[ind3]; ptB2[2] = z[ind3];
        ptC2[0] = x[ind1]; ptC2[1] = y[ind1]; ptC2[2] = z[ind1];
        dir2[0] = (ptB2[1]-ptA2[1])*(ptC2[2]-ptA2[2])-(ptB2[2]-ptA2[2])*(ptC2[1]-ptA2[1]);
        dir2[1] = (ptB2[2]-ptA2[2])*(ptC2[0]-ptA2[0])-(ptB2[0]-ptA2[0])*(ptC2[2]-ptA2[2]);
        dir2[2] = (ptB2[0]-ptA2[0])*(ptC2[1]-ptA2[1])-(ptB2[1]-ptA2[1])*(ptC2[0]-ptA2[0]);
        ndirl = sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2]);
        if (ndirl < 1.e-11) printf("convertTri2Quad:" SF_D_ ": quad=" SF_D_ " is 4-tri degen with edge contracted\n", ntr, i); 
        ntr++; 
        goto end0; 
      }
      else // decoupage
      {
        t1 = floor(indic[ind1]+0.5); t2 = floor(indic[ind2]+0.5); 
        t3 = floor(indic[ind3]+0.5); t4 = floor(indic[ind4]+0.5);
        if (t1 == 1 && t2 == 0 && t3 == 1 && t4 == 0) // 2
        {
          ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind3+1;
          ct1[ntr+1] = ind1+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
          ntr += 2; 
          ptA2[0] = x[ind2]; ptA2[1] = y[ind2]; ptA2[2] = z[ind2];
          ptB2[0] = x[ind3]; ptB2[1] = y[ind3]; ptB2[2] = z[ind3];
          ptC2[0] = x[ind1]; ptC2[1] = y[ind1]; ptC2[2] = z[ind1];
          dir2[0] = (ptB2[1]-ptA2[1])*(ptC2[2]-ptA2[2])-(ptB2[2]-ptA2[2])*(ptC2[1]-ptA2[1]);
          dir2[1] = (ptB2[2]-ptA2[2])*(ptC2[0]-ptA2[0])-(ptB2[0]-ptA2[0])*(ptC2[2]-ptA2[2]);
          dir2[2] = (ptB2[0]-ptA2[0])*(ptC2[1]-ptA2[1])-(ptB2[1]-ptA2[1])*(ptC2[0]-ptA2[0]);
          ndir2 = sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2]);
          ptA4[0] = x[ind4]; ptA4[1] = y[ind4]; ptA4[2] = z[ind4];
          ptB4[0] = x[ind1]; ptB4[1] = y[ind1]; ptB4[2] = z[ind1];
          ptC4[0] = x[ind3]; ptC4[1] = y[ind3]; ptC4[2] = z[ind3];
          dir4[0] = (ptB4[1]-ptA4[1])*(ptC4[2]-ptA4[2])-(ptB4[2]-ptA4[2])*(ptC4[1]-ptA4[1]);
          dir4[1] = (ptB4[2]-ptA4[2])*(ptC4[0]-ptA4[0])-(ptB4[0]-ptA4[0])*(ptC4[2]-ptA4[2]);
          dir4[2] = (ptB4[0]-ptA4[0])*(ptC4[1]-ptA4[1])-(ptB4[1]-ptA4[1])*(ptC4[0]-ptA4[0]);
          ndir4 = sqrt(dir4[0]*dir4[0]+dir4[1]*dir4[1]+dir4[2]*dir4[2]);
          if (ndir2 < 1.e-11 && ndir4 < 1.e-11) printf("convertTri2Quad: generate tri degen when forcing indic\n");
          goto end0;
        }
        else if (t1 == 0 && t2 == 1 && t3 == 0 && t4 == 1) // 2
        {
          ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind4+1;
          ct1[ntr+1] = ind2+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
          ntr += 2; 
          ptA1[0] = x[ind1]; ptA1[1] = y[ind1]; ptA1[2] = z[ind1];
          ptB1[0] = x[ind2]; ptB1[1] = y[ind2]; ptB1[2] = z[ind2];
          ptC1[0] = x[ind4]; ptC1[1] = y[ind4]; ptC1[2] = z[ind4];
          dir1[0] = (ptB1[1]-ptA1[1])*(ptC1[2]-ptA1[2])-(ptB1[2]-ptA1[2])*(ptC1[1]-ptA1[1]);
          dir1[1] = (ptB1[2]-ptA1[2])*(ptC1[0]-ptA1[0])-(ptB1[0]-ptA1[0])*(ptC1[2]-ptA1[2]);
          dir1[2] = (ptB1[0]-ptA1[0])*(ptC1[1]-ptA1[1])-(ptB1[1]-ptA1[1])*(ptC1[0]-ptA1[0]);
          ndir1 = sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
          ptA3[0] = x[ind3]; ptA3[1] = y[ind3]; ptA3[2] = z[ind3];
          ptB3[0] = x[ind4]; ptB3[1] = y[ind4]; ptB3[2] = z[ind4];
          ptC3[0] = x[ind2]; ptC3[1] = y[ind2]; ptC3[2] = z[ind2];
          dir3[0] = (ptB3[1]-ptA3[1])*(ptC3[2]-ptA3[2])-(ptB3[2]-ptA3[2])*(ptC3[1]-ptA3[1]);
          dir3[1] = (ptB3[2]-ptA3[2])*(ptC3[0]-ptA3[0])-(ptB3[0]-ptA3[0])*(ptC3[2]-ptA3[2]);
          dir3[2] = (ptB3[0]-ptA3[0])*(ptC3[1]-ptA3[1])-(ptB3[1]-ptA3[1])*(ptC3[0]-ptA3[0]);
          ndir3 = sqrt(dir3[0]*dir3[0]+dir3[1]*dir3[1]+dir3[2]*dir3[2]);
          if (ndir1 < 1.e-11 && ndir3 < 1.e-11) printf("convertTri2Quad: generate A tri degen when forcing indic\n");
          goto end0;
        }
        else
        {
          ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind4+1;
          ct1[ntr+1] = ind2+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
          ntr += 2;
        }
      }
      end0:;
    }
    if (ntr != nt)
    {
      ct.reAllocMat(ntr, 3);
    }

    // Deuxieme passe
    //K_TRANSFORM::flipEdges(ct, f->getSize(), x, y, z, NULL);

    // Troisieme passe : decoupage suivant indic si possible
    for (E_Int i = 0; i < 0; i++)
    {
      ind1 = cnp1[i]-1; ind2 = cnp2[i]-1; ind3 = cnp3[i]-1; ind4 = cnp4[i]-1;

      t1 = floor(indic[ind1]+0.5); t2 = floor(indic[ind2]+0.5); 
      t3 = floor(indic[ind3]+0.5); t4 = floor(indic[ind4]+0.5);
      if (t1 == 1 && t2 == 0 && t3 == 1 && t4 == 0) // 2
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind3+1;
        ct1[ntr+1] = ind1+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
        ntr += 2; 
        ptA2[0] = x[ind2]; ptA2[1] = y[ind2]; ptA2[2] = z[ind2];
        ptB2[0] = x[ind3]; ptB2[1] = y[ind3]; ptB2[2] = z[ind3];
        ptC2[0] = x[ind1]; ptC2[1] = y[ind1]; ptC2[2] = z[ind1];
        dir2[0] = (ptB2[1]-ptA2[1])*(ptC2[2]-ptA2[2])-(ptB2[2]-ptA2[2])*(ptC2[1]-ptA2[1]);
        dir2[1] = (ptB2[2]-ptA2[2])*(ptC2[0]-ptA2[0])-(ptB2[0]-ptA2[0])*(ptC2[2]-ptA2[2]);
        dir2[2] = (ptB2[0]-ptA2[0])*(ptC2[1]-ptA2[1])-(ptB2[1]-ptA2[1])*(ptC2[0]-ptA2[0]);
        ndir2 = sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2]);
        ptA4[0] = x[ind4]; ptA4[1] = y[ind4]; ptA4[2] = z[ind4];
        ptB4[0] = x[ind1]; ptB4[1] = y[ind1]; ptB4[2] = z[ind1];
        ptC4[0] = x[ind3]; ptC4[1] = y[ind3]; ptC4[2] = z[ind3];
        dir4[0] = (ptB4[1]-ptA4[1])*(ptC4[2]-ptA4[2])-(ptB4[2]-ptA4[2])*(ptC4[1]-ptA4[1]);
        dir4[1] = (ptB4[2]-ptA4[2])*(ptC4[0]-ptA4[0])-(ptB4[0]-ptA4[0])*(ptC4[2]-ptA4[2]);
        dir4[2] = (ptB4[0]-ptA4[0])*(ptC4[1]-ptA4[1])-(ptB4[1]-ptA4[1])*(ptC4[0]-ptA4[0]);
        ndir4 = sqrt(dir4[0]*dir4[0]+dir4[1]*dir4[1]+dir4[2]*dir4[2]);
        if (ndir2 < 1.e-11 && ndir4 < 1.e-11) printf("B tri degen\n");
        goto end;
      }
      else if (t1 == 0 && t2 == 1 && t3 == 0 && t4 == 1) // 2
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind4+1;
        ct1[ntr+1] = ind2+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
        ntr += 2; 
        ptA1[0] = x[ind1]; ptA1[1] = y[ind1]; ptA1[2] = z[ind1];
        ptB1[0] = x[ind2]; ptB1[1] = y[ind2]; ptB1[2] = z[ind2];
        ptC1[0] = x[ind4]; ptC1[1] = y[ind4]; ptC1[2] = z[ind4];
        dir1[0] = (ptB1[1]-ptA1[1])*(ptC1[2]-ptA1[2])-(ptB1[2]-ptA1[2])*(ptC1[1]-ptA1[1]);
        dir1[1] = (ptB1[2]-ptA1[2])*(ptC1[0]-ptA1[0])-(ptB1[0]-ptA1[0])*(ptC1[2]-ptA1[2]);
        dir1[2] = (ptB1[0]-ptA1[0])*(ptC1[1]-ptA1[1])-(ptB1[1]-ptA1[1])*(ptC1[0]-ptA1[0]);
        ndir1 = sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
        ptA3[0] = x[ind3]; ptA3[1] = y[ind3]; ptA3[2] = z[ind3];
        ptB3[0] = x[ind4]; ptB3[1] = y[ind4]; ptB3[2] = z[ind4];
        ptC3[0] = x[ind2]; ptC3[1] = y[ind2]; ptC3[2] = z[ind2];
        dir3[0] = (ptB3[1]-ptA3[1])*(ptC3[2]-ptA3[2])-(ptB3[2]-ptA3[2])*(ptC3[1]-ptA3[1]);
        dir3[1] = (ptB3[2]-ptA3[2])*(ptC3[0]-ptA3[0])-(ptB3[0]-ptA3[0])*(ptC3[2]-ptA3[2]);
        dir3[2] = (ptB3[0]-ptA3[0])*(ptC3[1]-ptA3[1])-(ptB3[1]-ptA3[1])*(ptC3[0]-ptA3[0]);
        ndir3 = sqrt(dir3[0]*dir3[0]+dir3[1]*dir3[1]+dir3[2]*dir3[2]);
        if (ndir1 < 1.e-11 && ndir3 < 1.e-11) printf("A tri degen\n");
        goto end;
      }

      // Angle 2D - 12-14
      ptA1[0] = x[ind1]; ptA1[1] = y[ind1]; ptA1[2] = z[ind1];
      ptB1[0] = x[ind2]; ptB1[1] = y[ind2]; ptB1[2] = z[ind2];
      ptC1[0] = x[ind4]; ptC1[1] = y[ind4]; ptC1[2] = z[ind4];
      dir1[0] = (ptB1[1]-ptA1[1])*(ptC1[2]-ptA1[2])-(ptB1[2]-ptA1[2])*(ptC1[1]-ptA1[1]);
      dir1[1] = (ptB1[2]-ptA1[2])*(ptC1[0]-ptA1[0])-(ptB1[0]-ptA1[0])*(ptC1[2]-ptA1[2]);
      dir1[2] = (ptB1[0]-ptA1[0])*(ptC1[1]-ptA1[1])-(ptB1[1]-ptA1[1])*(ptC1[0]-ptA1[0]);
      ndir1 = sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
        
      // Angle 2D - 23-21
      ptA2[0] = x[ind2]; ptA2[1] = y[ind2]; ptA2[2] = z[ind2];
      ptB2[0] = x[ind3]; ptB2[1] = y[ind3]; ptB2[2] = z[ind3];
      ptC2[0] = x[ind1]; ptC2[1] = y[ind1]; ptC2[2] = z[ind1];
      dir2[0] = (ptB2[1]-ptA2[1])*(ptC2[2]-ptA2[2])-(ptB2[2]-ptA2[2])*(ptC2[1]-ptA2[1]);
      dir2[1] = (ptB2[2]-ptA2[2])*(ptC2[0]-ptA2[0])-(ptB2[0]-ptA2[0])*(ptC2[2]-ptA2[2]);
      dir2[2] = (ptB2[0]-ptA2[0])*(ptC2[1]-ptA2[1])-(ptB2[1]-ptA2[1])*(ptC2[0]-ptA2[0]);
      ndir2 = sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2]);
        
      // Angle 2D - 34-32
      ptA3[0] = x[ind3]; ptA3[1] = y[ind3]; ptA3[2] = z[ind3];
      ptB3[0] = x[ind4]; ptB3[1] = y[ind4]; ptB3[2] = z[ind4];
      ptC3[0] = x[ind2]; ptC3[1] = y[ind2]; ptC3[2] = z[ind2];
      dir3[0] = (ptB3[1]-ptA3[1])*(ptC3[2]-ptA3[2])-(ptB3[2]-ptA3[2])*(ptC3[1]-ptA3[1]);
      dir3[1] = (ptB3[2]-ptA3[2])*(ptC3[0]-ptA3[0])-(ptB3[0]-ptA3[0])*(ptC3[2]-ptA3[2]);
      dir3[2] = (ptB3[0]-ptA3[0])*(ptC3[1]-ptA3[1])-(ptB3[1]-ptA3[1])*(ptC3[0]-ptA3[0]);
      ndir3 = sqrt(dir3[0]*dir3[0]+dir3[1]*dir3[1]+dir3[2]*dir3[2]);
        
      // Angle 2D - 41-43
      ptA4[0] = x[ind4]; ptA4[1] = y[ind4]; ptA4[2] = z[ind4];
      ptB4[0] = x[ind1]; ptB4[1] = y[ind1]; ptB4[2] = z[ind1];
      ptC4[0] = x[ind3]; ptC4[1] = y[ind3]; ptC4[2] = z[ind3];
      dir4[0] = (ptB4[1]-ptA4[1])*(ptC4[2]-ptA4[2])-(ptB4[2]-ptA4[2])*(ptC4[1]-ptA4[1]);
      dir4[1] = (ptB4[2]-ptA4[2])*(ptC4[0]-ptA4[0])-(ptB4[0]-ptA4[0])*(ptC4[2]-ptA4[2]);
      dir4[2] = (ptB4[0]-ptA4[0])*(ptC4[1]-ptA4[1])-(ptB4[1]-ptA4[1])*(ptC4[0]-ptA4[0]);
      ndir4 = sqrt(dir4[0]*dir4[0]+dir4[1]*dir4[1]+dir4[2]*dir4[2]);
        
      if (ndir1 < 1.e-11 && ndir2 < 1.e-11 && ndir3 < 1.e-11 && ndir4 < 1.e-11) printf("elt=" SF_D_ " is degenerated\n", i);

      // Qualite de la coupe1 124-234
      deg1 = false; inverse1 = 0.;
      if (ndir1 < 1.e-11 || ndir3 < 1.e-11) deg1 = true;
      inverse1 = dir1[0]*dir3[0]+dir1[1]*dir3[1]+dir1[2]*dir3[2];
      
      // Qualite de la coupe1 124-234
      deg2 = false; inverse2 = 0.;
      if (ndir2 < 1.e-11 || ndir4 < 1.e-11) deg2 = true;
      inverse2 = dir2[0]*dir4[0]+dir2[1]*dir4[1]+dir2[2]*dir4[2];
      

      if (ndir1 > 1.e-11)
      {
        ndir1 = 1./ndir1;
        dir1[0] = dir1[0]*ndir1; dir1[1] = dir1[1]*ndir1; dir1[2] = dir1[2]*ndir1;
      }
      else if (ndir2 > 1.e-11)
      {
        ndir2 = 1./K_FUNC::E_max(ndir2, 1.e-10);
        dir1[0] = dir2[0]*ndir2; dir1[1] = dir2[1]*ndir2; dir1[2] = dir2[2]*ndir2;
      }
      else if (ndir3 > 1.e-11)
      {
        ndir3 = 1./K_FUNC::E_max(ndir3, 1.e-10);
        dir1[0] = dir3[0]*ndir3; dir1[1] = dir3[1]*ndir3; dir1[2] = dir3[2]*ndir3;
      }
      else
      {
        ndir4 = 1./K_FUNC::E_max(ndir4, 1.e-10);
        dir1[0] = dir4[0]*ndir4; dir1[1] = dir4[1]*ndir4; dir1[2] = dir4[2]*ndir4;
      }
        
      if (ndir2 > 1.e-11)
      {
        ndir2 = 1./ndir2;
        dir2[0] = dir2[0]*ndir2; dir2[1] = dir2[1]*ndir2; dir2[2] = dir2[2]*ndir2;
      }
      else
      {
        dir2[0] = dir1[0]; dir2[1] = dir1[1]; dir2[2] = dir1[2];
      } 

      if (ndir3 > 1.e-11)
      {
        ndir3 = 1./ndir3;
        dir3[0] = dir3[0]*ndir3; dir3[1] = dir3[1]*ndir3; dir3[2] = dir3[2]*ndir3;
      }
      else
      {
        dir3[0] = dir1[0]; dir3[1] = dir1[1]; dir3[2] = dir1[2];
      } 

      if (ndir4 > 1.e-11)
      {
        ndir4 = 1./ndir4;
        dir4[0] = dir4[0]*ndir4; dir4[1] = dir4[1]*ndir4; dir4[2] = dir4[2]*ndir4;
      }
      else 
      {
        dir4[0] = dir1[0]; dir4[1] = dir1[1]; dir4[2] = dir1[2];
      }

      //alpha1 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA1, ptB1, ptA1, ptC1, dir1);
      //alpha2 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA2, ptB2, ptA2, ptC2, dir2);
      //alpha3 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA3, ptB3, ptA3, ptC3, dir3);
      //alpha4 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA4, ptB4, ptA4, ptC4, dir4);

      /*
      printf("elt=" SF_D_ ", ptA=%f %f %f, ptB=%f %f %f, ptC=%f %f %f, ptD=%f %f %f\n",i,
              x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],
              x[ind3],y[ind3],z[ind3],x[ind4],y[ind4],z[ind4]);
      printf("elt=" SF_D_ ", dir1=%f %f %f, dir2=%f %f %f, dir3=%f %f %f, dir4=%f %f %f\n", i, 
              dir1[0], dir1[1], dir1[2], dir2[0], dir2[1], dir2[2], 
              dir3[0], dir3[1], dir3[2], dir4[0], dir4[1], dir4[2]);
      printf("elt=" SF_D_ ", alpha= %f %f %f %f\n", i, alpha1, alpha2, alpha3, alpha4);
      */
      /*
      if (alpha1 > alpha2 || alpha1 > alpha4 || alpha3 > alpha2 || alpha3 > alpha4) // coupe le plus obtu
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind3+1;
        ct1[ntr+1] = ind1+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
        ntr += 2;
      }
      else // decoupage suivant 2
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind4+1;
        ct1[ntr+1] = ind2+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
        ntr += 2;
      }
      */

      if (deg2 == false && inverse2 > 1.e-10)
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind3+1;
        ct1[ntr+1] = ind1+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
        ntr += 2; 
      }
      else if (deg1 == false && inverse1 > 1.e-10)
      {
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind4+1;
        ct1[ntr+1] = ind2+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
        ntr += 2;
      }
      else
      {
        printf("deg or inverse split\n");
        ct1[ntr] = ind1+1; ct2[ntr] = ind2+1; ct3[ntr] = ind4+1;
        ct1[ntr+1] = ind2+1; ct2[ntr+1] = ind3+1; ct3[ntr+1] = ind4+1;
        ntr += 2; 
      }
      end:;
    }
  }

  // Redim 
  /*
  if (ntr != nt)
  {
    ct.reAllocMat(ntr, 3);
  } */

  //E_Float* indic = NULL;
  //if (posi > 0) indic = f->begin(posi);
  
  for (E_Int n = 0; n < 0; n++)
  {
    printf("iteration " SF_D_ "=================\n",n);
    //K_TRANSFORM::flipEdges(ct, f->getSize(), x, y, z,indic);
  }
  PyObject* tpl = K_ARRAY::buildArray3(*f, varString, ct, "TRI", api);
  RELEASESHAREDU(pArray, f, cn);
  return tpl;
}
