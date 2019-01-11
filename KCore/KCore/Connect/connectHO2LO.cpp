/*    
    Copyright 2013-2019 Onera.

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

//=============================================================================
/*
  Change a HO connectivity to a LO connectivity.
  Two modes: coarse or fine.
*/
//=============================================================================
E_Int K_CONNECT::connectHO2LO(const char* eltTypeHO,
                              FldArrayI& cEVHO,
                              char* eltTypeLO,
                              FldArrayI& cEVLO)
{
  
}

// Converti une connectivite HO par elements en connectivite LO coarse
// cEVLO doit etre deja alloue
E_Int K_CONNECT::connectHO2LOCoarse(const char* eltTypeHO,
                                    FldArrayI& cEVHO,
                                    FldArrayI& cEVLO)
{
  // Calcul des strides
  E_Int strideEltHO = 0; // stride pour avancer d'un element
  E_Int strideVertexHO = 0; // stride pour avancer d'un vertex dans un element
  E_Int strideElt = 0; // stride pour avancer d'un element
  E_Int strideVertex = 0; // stride pour avancer d'un vertex dans un element
  E_Int nkeep = 0; // vertex a conserver dans la connectivite HO
  E_Int nskip = 0; // nbre de vertex a dumper dans la connectivite HO

  E_Int nelts = cEVHO.getSize();
  E_Int nvpeHO = cEVHO.getNfld();
  E_Int stride = cEVHO.getStride();

  E_Int nvpe = cEVLO;
  // keep first
  for (E_Int i = 0; i < nelts; i++) 
    for (E_Int n = 1; n <= nvpe; n++) cEVLO(i,n) = cEVHO(i,n);
}

// Retourne toujours une connectivite TRI ou TETRA
E_Int K_CONNECT::connectHO2LOFine(const char* eltTypeHO,
                                  FldArrayI& cEVHO,
                                  char* eltTypeLO,
                                  FldArrayI& cEVLO)
{
  
}
