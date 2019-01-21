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


  for (E_Int i = 0; i < ne; i++)
  {
#ifdef BINARY
    // type d'element
    READI; type = ti;

    // follow
    READI; 

    // ntag
    READI; tagl = ti;

    // Indirection
    READI;

    // ntags
    for (E_Int j = 0; j < tagl; j++) READI;
#else
    // Indirection
    READI;
    
    // type d'element
    READI; type = ti;

    // tag
    READI; tagl = ti;
    for (E_Int j = 0; j < tagl; j++) READI;
#endif
    
    switch (type)
    {
      case 1: // BAR
        nBAR++;
        for (E_Int j = 0; j < 2; j++) READI;
        break;

      case 2: // TRI
        nTRI++;
        for (E_Int j = 0; j < 3; j++) READI;
        break;

      case 3: // QUAD
        nQUAD++;
        for (E_Int j = 0; j < 4; j++) READI;
        break;

      case 4: // TETRA
        nTETRA++;
        for (E_Int j = 0; j < 4; j++) READI;
        break;
        
      case 5: // HEXA
        nHEXA++;
        for (E_Int j = 0; j < 8; j++) READI;
        break;

      case 6: // PENTA
        nPENTA++;
        for (E_Int j = 0; j < 6; j++) READI;
        break;

      case 7: // PYRA
        nPYRA++;
        for (E_Int j = 0; j < 5; j++) READI;
        break;

      case 8: // 3-node BAR (second order)
        nBAR_3++;
        for (E_Int j = 0; j < 3; j++) READI;
        break;
        
      case 9: // 6-node TRI (second order)
        nTRI_6++;
        for (E_Int j = 0; j < 6; j++) READI;
        break;

      case 10: // 9-node QUAD (second order)
        nQUAD_9++;
        for (E_Int j = 0; j < 9; j++) READI;
        break;

      case 11: // 10-node TETRA (second order)
        nTETRA_10++;
        for (E_Int j = 0; j < 10; j++) READI;
        break;

      case 12: // 27-node HEXA (second order)
        nHEXA_27++;
        for (E_Int j = 0; j < 27; j++) READI;
        break;

      case 13: // 18-node PENTA (second order)
        nPENTA_18++;
        for (E_Int j = 0; j < 18; j++) READI;
        break;

      case 14: // 14-node PYRA (second order)
        nPYRA_14++;
        for (E_Int j = 0; j < 14; j++) READI;
        break;

      case 15: // NODE
        nNODE++;
        for (E_Int j = 0; j < 1; j++) READI;
        break;

      case 16: // 8-node QUAD (second order)
        nQUAD_8++;
        for (E_Int j = 0; j < 8; j++) READI;
        break;

      case 17: // 20-node HEXA (second order)
        nHEXA_20++;
        for (E_Int j = 0; j < 20; j++) READI;
        break;

      case 18: // 15-node PENTA (second order)
        nPENTA_15++;
        for (E_Int j = 0; j < 15; j++) READI;
        break;

      case 19: // 13-node PYRA (second order)
        nPYRA_13++;
        for (E_Int j = 0; j < 13; j++) READI;
        break;

      case 20: // 9-node TRI (third order)
        nTRI_9++;
        for (E_Int j = 0; j < 9; j++) READI;
        break;
        
      case 21: // 10-node TRI (third order)
        nTRI_10++;
        for (E_Int j = 0; j < 10; j++) READI;
        break;

      case 22: // 12-node TRI (fourth order)
        nTRI_12++;
        for (E_Int j = 0; j < 12; j++) READI;
        break;

      case 23: // 15-node TRI (fourth order)
        nTRI_15++;
        for (E_Int j = 0; j < 15; j++) READI;
        break;

      case 24: // ??-node TRI (fifth order) -> pas en CGNS
        nTRI++;
        for (E_Int j = 0; j < 15; j++) READI;
        break;

      case 25: // 21-node TRI (fifth order) -> pas en CGNS
        nTRI++;
        for (E_Int j = 0; j < 21; j++) READI;
        break;

      case 26: // 4-node BAR (third order)
        nBAR_4++;
        for (E_Int j = 0; j < 4; j++) READI;
        break;

      case 27: // 5-node BAR (fourth order)
        nBAR_5++;
        for (E_Int j = 0; j < 5; j++) READI;
        break;

      case 28: // 6-node BAR (fifth order)-> pas en CGNS
        nBAR++;
        for (E_Int j = 0; j < 6; j++) READI;
        break;

      case 29: // 20-node TETRA (third order)
        nTETRA_20++;
        for (E_Int j = 0; j < 20; j++) READI;
        break;
        
      case 30: // 35-node TETRA (fourth order)
        nTETRA_35++;
        for (E_Int j = 0; j < 35; j++) READI;
        break;

      case 31: // 56-node TETRA (fifth order)-> pas en CGNS
        nTETRA++;
        for (E_Int j = 0; j < 56; j++) READI;
        break;
        
      case 92: // 64-node HEXA (third order)
        nHEXA_64++;
        for (E_Int j = 0; j < 64; j++) READI;
        break;

      case 93: // 125-node HEXA (fourth order)
        nHEXA_125++;
        for (E_Int j = 0; j < 125; j++) READI;
        break;

      default:
        nDiscard++;
    }
  }
  
