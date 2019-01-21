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
        READI; (*cnBAR)(nBAR,1) = indirNodes[ti-1];
        READI; (*cnBAR)(nBAR,2) = indirNodes[ti-1];
        nBAR++;
        break;

      case 2: // TRI
        READI; (*cnTRI)(nTRI,1) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,2) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,3) = indirNodes[ti-1];
        nTRI++;
        break;

      case 3: // QUAD
        READI; (*cnQUAD)(nQUAD,1) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,2) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,3) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,4) = indirNodes[ti-1];
        nQUAD++;
        break;

      case 4: // TETRA
        READI; (*cnTETRA)(nTETRA,1) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,2) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,3) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,4) = indirNodes[ti-1];
        nTETRA++;
        break;
        
      case 5: // HEXA
        READI; (*cnHEXA)(nHEXA,1) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,2) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,3) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,4) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,5) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,6) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,7) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,8) = indirNodes[ti-1];
        nHEXA++;
        break;

      case 6: // PENTA
        READI; (*cnPENTA)(nPENTA,1) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,2) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,3) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,4) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,5) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,6) = indirNodes[ti-1];
        nPENTA++;
        break;

      case 7: // PYRA
        READI; (*cnPYRA)(nPYRA,1) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,2) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,3) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,4) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,5) = indirNodes[ti-1];
        nPYRA++;
        break;

      case 8: // 3-node BAR (second order)
        for (E_Int j = 1; j <= 3; j++) { READI; (*cnBAR_3)(nBAR_3,j) = indirNodes[ti-1]; }
        nBAR_3++;
        break;
        
      case 9: // 6-node TRI (second order)
        for (E_Int j = 1; j <= 6; j++) { READI; (*cnTRI_6)(nTRI_6,j) = indirNodes[ti-1]; }
        nTRI_6++;
        break;

      case 10: // 9-node QUAD (second order)
        for (E_Int j = 1; j <= 9; j++) { READI; (*cnQUAD_9)(nQUAD_9,j) = indirNodes[ti-1]; }
        nQUAD_9++;
        break;

      case 11: // 10-node TETRA (second order)
        for (E_Int j = 1; j <= 10; j++) { READI; (*cnTETRA_10)(nTETRA_10,j) = indirNodes[ti-1]; }
        nTETRA_10++;
        break;

      case 12: // 27-node HEXA (second order)
        for (E_Int j = 1; j <= 27; j++) { READI; (*cnHEXA_27)(nHEXA_27,j) = indirNodes[ti-1]; }
        nHEXA_27++;
        break;

      case 13: // 18-node PENTA (second order)
        for (E_Int j = 1; j <= 18; j++) { READI; (*cnPENTA_18)(nPENTA_18,j) = indirNodes[ti-1]; }
        nPENTA_18++;
        break;

      case 14: // 14-node PYRA (second order)
        for (E_Int j = 1; j <= 14; j++) { READI; (*cnPYRA_14)(nPYRA_14,j) = indirNodes[ti-1]; }
        nPYRA_14++;
        break;

      case 15: // NODE
        READI; (*indNODE)[nNODE] = ti;
        nNODE++;
        break;

      case 16: // 8-node QUAD (second order)
        for (E_Int j = 1; j <= 8; j++) { READI; (*cnQUAD_8)(nQUAD_8,j) = indirNodes[ti-1]; }
        nQUAD_8++;
        break;

      case 17: // 20-node HEXA (second order)
        for (E_Int j = 1; j <= 20; j++) { READI; (*cnHEXA_20)(nHEXA_20,j) = indirNodes[ti-1]; }
        nHEXA_20++;
        break;

      case 18: // 15-node PENTA (second order)
        for (E_Int j = 1; j <= 15; j++) { READI; (*cnPENTA_15)(nPENTA_15,j) = indirNodes[ti-1]; }
        nPENTA_15++;
        break;

      case 19: // 13-node PYRA (second order)
        for (E_Int j = 1; j <= 13; j++) { READI; (*cnPYRA_13)(nPYRA_13,j) = indirNodes[ti-1]; }
        nPYRA_13++;
        break;

      case 20: // 9-node TRI (third order)
        for (E_Int j = 1; j <= 9; j++) { READI; (*cnTRI_9)(nTRI_9,j) = indirNodes[ti-1]; }
        nTRI_9++;
        break;
        
      case 21: // 10-node TRI (third order)
        for (E_Int j = 1; j <= 10; j++) { READI; (*cnTRI_10)(nTRI_10,j) = indirNodes[ti-1]; }
        nTRI_10++;
        break;

      case 22: // 12-node TRI (fourth order)
        for (E_Int j = 1; j <= 12; j++) { READI; (*cnTRI_12)(nTRI_12,j) = indirNodes[ti-1]; }
        nTRI_12++;
        break;

      case 23: // 15-node TRI (fourth order)
        for (E_Int j = 1; j <= 15; j++) { READI; (*cnTRI_15)(nTRI_15,j) = indirNodes[ti-1]; }
        nTRI_15++;
        break;

      case 24: // 15-node TRI (fifth order)
        for (E_Int j = 1; j <= 3; j++) { READI; (*cnTRI)(nTRI,j) = indirNodes[ti-1]; }
        for (E_Int j = 4; j <= 15; j++) { READI; };
        nTRI++;
        break;

      case 25: // 21-node TRI (fifth order)
        for (E_Int j = 1; j <= 3; j++) { READI; (*cnTRI)(nTRI,j) = indirNodes[ti-1]; }
        for (E_Int j = 4; j <= 21; j++) { READI; }; 
        nTRI++;
        break;

      case 26: // 4-node BAR (third order)
        for (E_Int j = 1; j <= 4; j++) { READI; (*cnBAR_4)(nBAR_4,j) = indirNodes[ti-1]; }
        nBAR_4++;
        break;

      case 27: // 5-node BAR (fourth order)
        for (E_Int j = 1; j <= 5; j++) { READI; (*cnBAR_5)(nBAR_5,j) = indirNodes[ti-1]; }
        nBAR_5++;
        break;

      case 28: // 6-node BAR (fifth order)
        for (E_Int j = 1; j <= 2; j++) { READI; (*cnBAR)(nBAR,j) = indirNodes[ti-1]; }
        for (E_Int j = 3; j <= 6; j++) READI;
        nBAR++;
        break;

      case 29: // 20-node TETRA (third order)
        for (E_Int j = 1; j <= 20; j++) { READI; (*cnTETRA_20)(nTETRA_20,j) = indirNodes[ti-1]; }
        nTETRA_20++;
        break;
        
      case 30: // 35-node TETRA (fourth order)
        for (E_Int j = 1; j <= 35; j++) { READI; (*cnTETRA_35)(nTETRA_35,j) = indirNodes[ti-1]; }
        nTETRA_35++;
        break;

      case 31: // 56-node TETRA (fifth order)
        for (E_Int j = 1; j <= 4; j++) { READI; (*cnTETRA)(nTETRA,j) = indirNodes[ti-1]; }
        for (E_Int j = 5; j <= 56; j++) READI;
        nTETRA++;
        break;
        
      case 92: // 64-node HEXA (third order)
        for (E_Int j = 1; j <= 64; j++) { READI; (*cnHEXA_64)(nHEXA_64,j) = indirNodes[ti-1]; }
        nHEXA_64++;
        break;

      case 93: // 125-node HEXA (fourth order)
        for (E_Int j = 1; j <= 125; j++) { READI; (*cnHEXA_125)(nHEXA_125,j) = indirNodes[ti-1]; }
        nHEXA_125++;
        break;

      default:
        nDiscard++;
    }
  }
