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

// Convert array to tetra array

# include "converter.h"
# include "kcore.h"
# include <vector>

using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Convert array to a tetraedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertArray2Tetra(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  res = K_ARRAY::getFromArray3(array, varString, f,
                               ni, nj, nk, cn, eltType);

  if (res == 1)
  {
    PyObject* tpl = convertStruct2Tetra(varString, f, ni, nj, nk);
    RELEASESHAREDS(array, f);
    return tpl; 
  }
  else if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2Tetra: array is invalid.");
    return NULL;
  }

  if (K_STRING::cmp(eltType, "NGON") == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2Tetra: NGON array not supported.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }

  // Check that at least one element type does not belong to the list of
  // elements that does not need to be converted, ie, {NODE, BAR, TRI, TETRA}
  E_Int nc = cn->getNConnect();
  E_Int nfld = f->getNfld();
  E_Int api = f->getApi();
  E_Int npts = f->getSize();

  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  E_Int dim = 3;
  E_Int nelts2 = 0;
  std::vector<E_Bool> convConn(nc, false);
  E_Bool needCoords = false;

  // Offsets of the input ME
  std::vector<E_Int> nepc(nc);
  std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic] = nelts;

    if (K_STRING::cmp(eltTypes[ic], "QUAD") == 0)
    { convConn[ic] = true; nelts2 += 2*nelts; dim = 2; }
    else if (K_STRING::cmp(eltTypes[ic], "PYRA") == 0)
    { convConn[ic] = true; nelts2 += 2*nelts; needCoords = true; }
    else if (K_STRING::cmp(eltTypes[ic], "PENTA") == 0)
    { convConn[ic] = true; nelts2 += 3*nelts; }
    else if (K_STRING::cmp(eltTypes[ic], "HEXA") == 0)
    { convConn[ic] = true; nelts2 += 5*nelts; }
    else { nelts2 += nelts; }
  }
  for (E_Int ic = 0; ic < nc; ic++) cumnepc[ic+1] = cumnepc[ic] + nepc[ic];

  E_Bool foundEltType2Convert = false;
  for (E_Int ic = 0; ic < nc; ic++) foundEltType2Convert |= convConn[ic];
  if (!foundEltType2Convert)
  {
    // Nothing to convert
    RELEASESHAREDU(array, f, cn);
    Py_INCREF(array);
    return array;
  }

  // Coordinates required for PYRA -> TETRA
  E_Float *xp = NULL, *yp = NULL, *zp = NULL;
  if (needCoords)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "convertArray2Tetra: coords must be present in array.");
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }
    posx++; posy++; posz++;
    xp = f->begin(posx);
    yp = f->begin(posy);
    zp = f->begin(posz);
  }

  // Build new connectivity
  char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH]; eltType2[0] = '\0';
  if (dim == 2) strcpy(eltType2, "TRI");
  else strcpy(eltType2, "TETRA");

  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts2,
                                       eltType2, false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  #pragma omp parallel if (nelts2 > __MIN_SIZE_MEAN__)
  {
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      FldArrayI& cm2 = *(cn2->getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int offset = cumnepc[ic];

      if (!convConn[ic])  // simple copy
      {
        E_Int nvpe = cm.getNfld();
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
          for (E_Int j = 1; j <= nvpe; j++) cm2(i+offset,j) = cm(i,j);
      }
      else if (dim == 2)  // QUAD
      {
        E_Int ind, ind1, ind2, ind3, ind4;
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i,1); 
          ind2 = cm(i,2);
          ind3 = cm(i,3);
          ind4 = cm(i,4);
          
          ind = 2*i + offset;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind2;
          cm2(ind,3) = ind3;
          
          ind++;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind4;
        }
      }
      else if (K_STRING::cmp(eltTypes[ic], "PYRA") == 0)
      {
        E_Float l13x, l13y, l13z, l24x, l24y, l24z;
        E_Float sqrDiag1, sqrDiag2;
        E_Int cnt, i1, i2, i3, i4;

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          /* determination des indices de l'element */
          i1 = cm(i,1) - 1;
          i2 = cm(i,2) - 1;
          i3 = cm(i,3) - 1;
          i4 = cm(i,4) - 1;

          l13x = xp[i3] - xp[i1]; l24x = xp[i4] - xp[i2];
          l13y = yp[i3] - yp[i1]; l24y = yp[i4] - yp[i2];
          l13z = zp[i3] - zp[i1]; l24z = zp[i4] - zp[i2];

          /* determination de la plus petite diagonale de la face quad i1i2i3i4 */
          // on retient les 2 sommets de la diag min
          sqrDiag1 = l13x * l13x + l13y * l13y + l13z * l13z;
          sqrDiag2 = l24x * l24x + l24y * l24y + l24z * l24z;

          /* construction des elements tetra */
          if (sqrDiag1 <= sqrDiag2)
          {
            // build tetras: I1I2I3I5,I1I3I4I5
            // t1: I1I2I3I5
            cnt = 2*i + offset;
            cm2(cnt,1) = cm(i,1);
            cm2(cnt,2) = cm(i,2);
            cm2(cnt,3) = cm(i,3);
            cm2(cnt,4) = cm(i,5);

            // t2: I1I3I4I5
            cnt = 2*i+1;
            cm2(cnt,1) = cm(i,1);
            cm2(cnt,2) = cm(i,3);
            cm2(cnt,3) = cm(i,4);
            cm2(cnt,4) = cm(i,5);
          }
          else
          {
            // build tetras: I2I3I4I5, I2I4I1I5
            // t1: I1I2I3I5
            cnt = 2*i + offset;
            cm2(cnt,1) = cm(i,2);
            cm2(cnt,2) = cm(i,3);
            cm2(cnt,3) = cm(i,4);
            cm2(cnt,4) = cm(i,5);

            // t2: I1I5I3I6
            cnt = 2*i+1;
            cm2(cnt,1) = cm(i,2);
            cm2(cnt,2) = cm(i,4);
            cm2(cnt,3) = cm(i,1);
            cm2(cnt,4) = cm(i,5);
          }
        }
      }
      else if (K_STRING::cmp(eltTypes[ic], "PENTA") == 0)
      {
        E_Int cnt, i1, i2, i3, i4, i5, i6, diag;
        E_Int indir[6];

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          /* determination du prisme tourne: imin -> pt 1
            determination du second point min sur la facette quad opposee: diag
            retour du tableau d'indirection I1I2I3I4I5I6 */
          diag = 0;
          buildSortedPrism(i, cm, diag, indir);
        
          i1 = indir[0]; i2 = indir[1];
          i3 = indir[2]; i4 = indir[3];
          i5 = indir[4]; i6 = indir[5];
          if (diag == -1) //config 2-6
          {
            // build tetras: I1I2I3I6,I1I2I6I5,I1I5I6I4 
            // t1: I1I2I3I6
            cnt = 3*i + offset;
            cm2(cnt,1) = cm(i,i1);
            cm2(cnt,2) = cm(i,i2);
            cm2(cnt,3) = cm(i,i3);
            cm2(cnt,4) = cm(i,i6);

            // t2: I1I2I6I5
            cnt = 3*i+1;
            cm2(cnt,1) = cm(i,i1);
            cm2(cnt,2) = cm(i,i2);
            cm2(cnt,3) = cm(i,i6);
            cm2(cnt,4) = cm(i,i5);

            // t3: I1I5I6I4
            cnt = 3*i+2;
            cm2(cnt,1) = cm(i,i1);
            cm2(cnt,2) = cm(i,i5);
            cm2(cnt,3) = cm(i,i6);
            cm2(cnt,4) = cm(i,i4);
          }
          else if (diag == 1)//config 3-5
          {
            // build tetras: I1I2I3I5, I1I5I3I6, I1I5I6I4
            // t1: I1I2I3I5
            cnt = 3*i + offset;
            cm2(cnt,1) = cm(i,i1);
            cm2(cnt,2) = cm(i,i2);
            cm2(cnt,3) = cm(i,i3);
            cm2(cnt,4) = cm(i,i5);

            // t2: I1I5I3I6
            cnt = 3*i+1;
            cm2(cnt,1) = cm(i,i1);
            cm2(cnt,2) = cm(i,i5);
            cm2(cnt,3) = cm(i,i3);
            cm2(cnt,4) = cm(i,i6);

            // t3: I1I5I6I4
            cnt = 3*i+2;
            cm2(cnt,1) = cm(i,i1);
            cm2(cnt,2) = cm(i,i5);
            cm2(cnt,3) = cm(i,i6);
            cm2(cnt,4) = cm(i,i4);
          }
        }
      }
      else if (K_STRING::cmp(eltTypes[ic], "HEXA") == 0)
      {
        E_Int ind, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i,1);
          ind2 = cm(i,2);
          ind3 = cm(i,3);
          ind4 = cm(i,4);
          ind5 = cm(i,5);
          ind6 = cm(i,6);
          ind7 = cm(i,7);
          ind8 = cm(i,8);
          
          if (i%2 == 0) // pair 
          {
            //tetra ABDE
            ind = 5*i + offset;
            cm2(ind,1) = ind1;
            cm2(ind,2) = ind2;
            cm2(ind,3) = ind4;
            cm2(ind,4) = ind5;
            
            //tetra BCDG
            ind++;
            cm2(ind,1) = ind2;
            cm2(ind,2) = ind3;
            cm2(ind,3) = ind4;
            cm2(ind,4) = ind7;
            
            //tetra DEGH
            ind++;
            cm2(ind,1) = ind4;
            cm2(ind,2) = ind5;
            cm2(ind,3) = ind7;
            cm2(ind,4) = ind8;
            
            //tetra BEFG
            ind++;
            cm2(ind,1) = ind2;
            cm2(ind,2) = ind5;
            cm2(ind,3) = ind6;
            cm2(ind,4) = ind7;
            
            //tetra BDEG
            ind++;
            cm2(ind,1) = ind2;
            cm2(ind,2) = ind4;
            cm2(ind,3) = ind5;
            cm2(ind,4) = ind7;
          }
          else // impair 
          {
            //tetra ACDH : 1348
            ind = 5*i + offset;
            cm2(ind,1) = ind1;
            cm2(ind,2) = ind3;
            cm2(ind,3) = ind4;
            cm2(ind,4) = ind8;
            
            //tetra AFBC : 1623
            ind++;
            cm2(ind,1) = ind1;
            cm2(ind,2) = ind2;
            cm2(ind,3) = ind3;
            cm2(ind,4) = ind6;
            //cm2(ind,1) = ind1;
            //cm2(ind,2) = ind6;
            //cm2(ind,3) = ind2;
            //cm2(ind,4) = ind3;

            //tetra HFGC : 8673
            ind++;
            cm2(ind,1) = ind3;
            cm2(ind,2) = ind6;
            cm2(ind,3) = ind7;
            cm2(ind,4) = ind8;
            //cm2(ind,1) = ind8;
            //cm2(ind,2) = ind6;
            //cm2(ind,3) = ind7;
            //cm2(ind,4) = ind3;

            //tetra FHAE : 6815
            ind++;
            cm2(ind,1) = ind6;
            cm2(ind,2) = ind5;
            cm2(ind,3) = ind8;
            cm2(ind,4) = ind1;
            //cm2(ind,1) = ind6;
            //cm2(ind,2) = ind8;
            //cm2(ind,3) = ind1;
            //cm2(ind,4) = ind5;

            //tetra FHAC : 6813
            ind++;
            cm2(ind,1) = ind6;
            cm2(ind,2) = ind3;
            cm2(ind,3) = ind1;
            cm2(ind,4) = ind8;
            //cm2(ind,1) = ind6;
            //cm2(ind,2) = ind8;
            //cm2(ind,3) = ind1;
            //cm2(ind,4) = ind3;
          }
        }
      }
    }

    // Copy fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for nowait
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }

  RELEASESHAREDU(array, f, cn);
  RELEASESHAREDU(tpl, f2, cn2);
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  return tpl;
}

//=============================================================================
/* Determination du prisme ayant le sommet de plus petit indice en bas a gauche
   diag vaut -1 si deuxieme min est I2 ou I6
   diag vaut 1 si deuxieme min est I3 ou I5
 */
//=============================================================================
void K_CONVERTER::buildSortedPrism(E_Int elt, FldArrayI& cn, E_Int& diag,
                                   E_Int* indir)
{
  // Determination de indmin
  E_Int indmin = cn(elt,1);
  E_Int imin = 1;
  E_Int ind;
  for (E_Int i = 2; i <= 6; i++)
  {
    ind = cn(elt,i);
    if (ind < indmin)
    {
      indmin = ind;
      imin = i;
    }
  }
  switch (imin)
  {
    case 1 :
      indir[0] = 1;
      indir[1] = 2;
      indir[2] = 3;
      indir[3] = 4;
      indir[4] = 5;
      indir[5] = 6;
      break;

    case 2 : 
      indir[0] = 2;
      indir[1] = 3;
      indir[2] = 1;
      indir[3] = 5;
      indir[4] = 6;
      indir[5] = 4;
      break;

    case 3 : 
      indir[0] = 3;
      indir[1] = 1;
      indir[2] = 2;
      indir[3] = 6;
      indir[4] = 4;
      indir[5] = 5;
      break;

    case 4 : 
      indir[0] = 4;
      indir[1] = 6;
      indir[2] = 5;
      indir[3] = 1;
      indir[4] = 3;
      indir[5] = 2;
      break;

    case 5 : 
      indir[0] = 5;
      indir[1] = 4;
      indir[2] = 6;
      indir[3] = 2;
      indir[4] = 1;
      indir[5] = 3;
      break;

    case 6 : 
      indir[0] = 6;
      indir[1] = 5;
      indir[2] = 4;
      indir[3] = 3;
      indir[4] = 2;
      indir[5] = 1;
      break;

    default ://erreur de codage
      printf("Error: code error in function buildSortedPrism.\n");
      diag = 0;
      return;
  }
  //determination de l indice min sur la 3eme facette quad
  // soit I2, I6, I3, I5 
  E_Int indI2 = cn(elt,indir[1]);
  E_Int indI3 = cn(elt,indir[2]);
  E_Int indI5 = cn(elt,indir[4]);
  E_Int indI6 = cn(elt,indir[5]);
  
  indmin = indI2;
  diag = -1;

  if (indI6 < indmin)
  {
    indmin = indI6;
    diag = -1;
  } 
  if (indI3 < indmin) 
  {
    indmin = indI3;
    diag = 1;
  }
  if (indI5 < indmin)
  {
    indmin = indI5;
    diag = 1;
  }
}
