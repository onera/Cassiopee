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
# include "metric.h"
# include "Array/Array.h"

//=============================================================================
// Calcul les vecteurs normaux aux triangles
// IN: cn: connectivite elts-noeuds
// IN: xt, yt, zt: coordonnees x, y, z des pts de la grille
// OUT: nsurf: normales aux facettes
//=============================================================================
void K_METRIC::compNormUnstructSurf(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* nsurfx, E_Float* nsurfy, E_Float* nsurfz
)
{
  E_Int offset = 0;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
      
    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      #pragma omp parallel
      {
        E_Int ind1, ind2, ind3, pos;
        E_Float l1x, l1y, l1z;
        E_Float l2x, l2y, l2z;
        E_Float surfx, surfy, surfz;

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i, 1) - 1;
          ind2 = cm(i, 2) - 1;
          ind3 = cm(i, 3) - 1;

          l1x = xt[ind1] - xt[ind2];
          l1y = yt[ind1] - yt[ind2];
          l1z = zt[ind1] - zt[ind2];

          l2x = xt[ind1] - xt[ind3];
          l2y = yt[ind1] - yt[ind3];
          l2z = zt[ind1] - zt[ind3];

          surfx = l1y * l2z - l1z * l2y;
          surfy = l1z * l2x - l1x * l2z;
          surfz = l1x * l2y - l1y * l2x;

          pos = offset + i;
          nsurfx[pos] = K_CONST::ONE_HALF * surfx;
          nsurfy[pos] = K_CONST::ONE_HALF * surfy;
          nsurfz[pos] = K_CONST::ONE_HALF * surfz;
        }
      }
      offset += nelts;
    }
    else
    {
      fprintf(stderr, "Error in K_METRIC::compNormUnstructSurf.\n");
      fprintf(stderr, "Element type can be TRI only, not %s.\n", eltTypes[ic]);
      exit(0);
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

void K_METRIC::compNormUnstructSurft(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* nxt, E_Float* nyt, E_Float* nzt
)
{
  E_Int offset = 0;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    
    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      E_Int nelts = cm.getSize();

      #pragma omp parallel
      {
        E_Int ind1, ind2, ind3, pos;
        E_Float l1x, l1y, l1z;
        E_Float l2x, l2y, l2z;
        E_Float surfx, surfy, surfz;

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i, 1) - 1;
          ind2 = cm(i, 2) - 1;
          ind3 = cm(i, 3) - 1;

          l1x = xt[ind1] - xt[ind2];
          l1y = yt[ind1] - yt[ind2];
          l1z = zt[ind1] - zt[ind2];

          l2x = xt[ind1] - xt[ind3];
          l2y = yt[ind1] - yt[ind3];
          l2z = zt[ind1] - zt[ind3];

          surfx = l1y * l2z - l1z * l2y;
          surfy = l1z * l2x - l1x * l2z;
          surfz = l1x * l2y - l1y * l2x;
          
          pos = offset + i;
          nxt[pos] = K_CONST::ONE_HALF * surfx;
          nyt[pos] = K_CONST::ONE_HALF * surfy;
          nzt[pos] = K_CONST::ONE_HALF * surfz;
        }
      }
      offset += nelts;
    }
    else
    {
      fprintf(stderr, "Error in K_METRIC::compNormUnstructSurft.\n");
      fprintf(stderr, "Element type can be TRI only, not %s.\n", eltTypes[ic]);
      exit(0);
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}
