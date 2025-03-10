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

# include "loc.h"
# include "Connect/connect.h"
# include <vector>
# include <algorithm> 

using namespace K_FLD;

extern "C"
{
  void k6conv2center1_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                       const E_Int& nfld, E_Float* fieldnode, 
                       E_Float* fieldcenter);
  void k6conv2center21_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                        E_Float* fieldnode, 
                        E_Float* fieldcenter);
  void k6conv2center22_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                        E_Float* fieldnode, 
                        E_Float* fieldcenter);
  void k6conv2center12d_(const E_Int& ni, const E_Int& nj, 
                         const E_Int& nfld, E_Float* fieldnode, 
                         E_Float* fieldcenter);
  void k6conv2center212d_(const E_Int& ni, const E_Int& nj, 
                          E_Float* fieldnode, 
                          E_Float* fieldcenter);
  void k6conv2center222d_(const E_Int& ni, const E_Int& nj,
                          E_Float* fieldnode, 
                          E_Float* fieldcenter);
  void k6conv2center11d_(const E_Int& ni, const E_Int& nfld, 
                         E_Float* fieldnode, 
                         E_Float* fieldcenter);
  void k6conv2center211d_(const E_Int& ni, E_Float* fieldnode, 
                          E_Float* fieldcenter);
  void k6conv2center221d_(const E_Int& ni, E_Float* fieldnode, 
                          E_Float* fieldcenter);
}

//=============================================================================
// Convertit un array noeuds en array centres en structure
// Retourne 1 en cas de succes, 0 en cas d'echec.
//=============================================================================
E_Int K_LOC::node2centerStruct_OLD(FldArrayF& FNode, 
                               E_Int ni, E_Int nj, E_Int nk,
                               E_Int cellN, E_Int mod, 
                               FldArrayF& FCenter)
{
  E_Int nv = FNode.getNfld();
  E_Int size, dim;  
  E_Int im, jm, km;

  if (ni != 1 && nj == 1 && nk == 1 )
  {
    size = ni-1; dim = 1; im = ni; 
  }
  else if (ni == 1 && nj != 1 && nk == 1) 
  {
    size = nj-1; dim = 1; im = nj;
  }
  else if (ni == 1 && nj == 1 && nk != 1)
  {
    size = nk-1; dim = 1; im = nk;
  }
  else
  {
    if (ni == 1)
    {
      size = (nj-1)*(nk-1); dim = 2; im = nj; jm = nk;
    }
    else if (nj == 1)
    {
      size = (ni-1)*(nk-1); dim = 2; im = ni; jm = nk;
    }
    else if (nk == 1)
    {
      size = (ni-1)*(nj-1); dim = 2; im = ni; jm = nj;
    }
    else
    {
      size = (ni-1)*(nj-1)*(nk-1); dim = 3; im = ni; jm = nj; km = nk;
    }
  }
  
  // On alloue FCenter seulement s'il n'est pas deja alloue correctement
  if (FCenter.getSize() != size || FCenter.getNfld() != nv)
    FCenter.malloc(size, nv);

  // In of each "block" i, converts field defined in nodes to field 
  // defined in centers
  if (dim == 1)
    k6conv2center11d_(im, nv, FNode.begin(), FCenter.begin());
  
  else if (dim == 2)
    k6conv2center12d_(im, jm, nv, 
                      FNode.begin(), FCenter.begin());
  else
    k6conv2center1_(im, jm, km, nv, 
                    FNode.begin(), FCenter.begin());
  
  // If field contains "cellnaturefield"
  if (cellN != -1)
  {
    switch (mod)
    {
      case 1:
        if (dim == 1)
          k6conv2center211d_(im, FNode.begin(cellN), FCenter.begin(cellN));
        else if (dim == 2)
          k6conv2center212d_(im, jm, FNode.begin(cellN), 
                             FCenter.begin(cellN));
        else
          k6conv2center21_(im, jm, nk, FNode.begin(cellN), 
                           FCenter.begin(cellN));
        break;
        
      case 2:
        if (dim == 1)
          k6conv2center221d_(im, FNode.begin(cellN), FCenter.begin(cellN));
        else if (dim == 2)
          k6conv2center222d_(im, jm, FNode.begin(cellN), 
                             FCenter.begin(cellN));
        else
          k6conv2center22_(im, jm, nk, FNode.begin(cellN), 
                           FCenter.begin(cellN));
        break;
        
      case 3:
        printf("Warning: node2center: this cellN type is not implemented yet.\n");
        return 0;
        break;
        
      default:
        printf("Warning: node2center: unknown cellnaturefield format.\n");
        return 0;
    }
  }
  return 1;

}
//=============================================================================
// Convertit un array noeuds en array centres en non-structure
// Le traitement specifique cellN n'est pas implemente.
// Retourne 1 en cas de succes, 0 en cas d'echec.
// FCenter doit deja etre alloue au nb d'elements
//=============================================================================
E_Int K_LOC::node2centerUnstruct_OLD(FldArrayF& FNode, 
                                 FldArrayI& c,
                                 E_Int cellN, E_Int mod, 
                                 FldArrayF& FCenter)
{
  E_Int ne = c.getSize(); // nombre de centres = nombre d'elements
  E_Int nt = c.getNfld();
  E_Float ntinv = 1./nt;
  E_Int nfld = FNode.getNfld();
  FCenter.setAllValuesAtNull();
  // Les centres sont numerotes comme les elements
  for (E_Int v = 1; v <= nfld; v++)
  {
    E_Float* fnode = FNode.begin(v);
    E_Float* fcen = FCenter.begin(v);
    
    for (E_Int n = 1; n <= nt; n++)
    {
      E_Int* cn = c.begin(n);
      for (E_Int e = 0; e < ne; e++)
      {
        fcen[e] += fnode[cn[e]-1];
      }
    }
    for (E_Int e = 0; e < ne; e++) fcen[e] *= ntinv;
  }
  return 1;
}
//===============================================================================
// Convertit un champ en noeuds en champ en centres en NGON
// Retourne 1 en cas de succes, 0 en cas d'echec.
// FCenter doit deja etre alloue au nb d'elements
// sorted=1: vertices coordinates are sorted for better accuracy in summations 
//===============================================================================
E_Int K_LOC::node2centerNGon_OLD(FldArrayF& FNode, FldArrayI& cNG,
                             FldArrayF& FCenter, E_Int sorted)
{
  E_Int* cnp = cNG.begin();
  E_Int sizeFN = cnp[1];
  E_Int ncells = cnp[sizeFN+2]; 
  E_Int nfld = FNode.getNfld();
  FCenter.setAllValuesAtNull();

  if (sorted == 0)
  {
    std::vector< std::vector<E_Int> > cEV(ncells);
    K_CONNECT::connectNG2EV(cNG, cEV);

    for (E_Int v = 1; v <= nfld; v++)
    {
      E_Float* fnode = FNode.begin(v);
      E_Float* fcen = FCenter.begin(v);
      for (E_Int et = 0; et < ncells; et++)
      {
        std::vector<E_Int>& vertices = cEV[et]; // noeuds associes a l'element et
        E_Int nvert = vertices.size();
        E_Float ntinv = 1./nvert;
        for (E_Int nv = 0; nv < nvert; nv++)
        {
          E_Int indv = vertices[nv]-1;
          fcen[et] += fnode[indv];
        }
        fcen[et] *= ntinv;
      }// loop on elts
    }
  }
  else
  {
    E_Int* ptr = cNG.begin();
    FldArrayI posFace; K_CONNECT::getPosFaces(cNG, posFace);
    FldArrayI posElt; K_CONNECT::getPosElts(cNG, posElt);

# pragma omp parallel for default(shared)
    for (E_Int i = 0; i < ncells; i++)
    {
      E_Int pose = posElt[i];
      E_Int* ptrElt = &ptr[pose];
      E_Int nfaces = ptrElt[0]; 

      std::vector<E_Int> vertices; vertices.reserve(1024);
      for (E_Int n=1; n <= nfaces; n++)
      {
        E_Int ind = ptrElt[n]-1;
        E_Int pos = posFace[ind];
        E_Int* ptrFace = &ptr[pos];
        E_Int nv = ptrFace[0];
      
        for (E_Int p = 1; p <= nv; p++)
        {
          E_Int indv = ptrFace[p]-1;
          vertices.push_back(indv);
        }//loop on vertices
      }
      std::sort(vertices.begin(), vertices.end());
      vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end() );

      std::vector<E_Float> fsort(vertices.size());
      for (E_Int nofld = 1; nofld <= nfld; nofld++)
      {
        E_Float* fnode = FNode.begin(nofld);
        for (size_t nov = 0; nov < vertices.size(); nov++)
        {
          E_Int indv = vertices[nov]; fsort[nov] = fnode[indv];
        }
        std::sort(fsort.begin(), fsort.end());
        E_Float* fcen = FCenter.begin(nofld);
        for (size_t nov = 0; nov < fsort.size(); nov++) {fcen[i] += fsort[nov];}
        E_Float inv = 1./E_Float(fsort.size()); fcen[i] *= inv;   
      }
    } 
  }
  return 1;
}
