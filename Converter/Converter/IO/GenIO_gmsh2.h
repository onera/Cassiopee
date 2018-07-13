/*    
    Copyright 2013-2018 Onera.

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
        READI; (*cnBAR)(nBAR,1) = indirNodes[ti-1];
        READI; (*cnBAR)(nBAR,2) = indirNodes[ti-1];
        for (E_Int j = 2; j < 3; j++) READI;
        nBAR++;
        break;
        
      case 9: // 6-node TRI (second order)
        READI; (*cnTRI)(nTRI,1) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,2) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,3) = indirNodes[ti-1];
        for (E_Int j = 3; j < 6; j++) READI;
        nTRI++;
        break;

      case 10: // 9-node QUAD (second order)
        READI; (*cnQUAD)(nQUAD,1) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,2) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,3) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,4) = indirNodes[ti-1];
        for (E_Int j = 4; j < 9; j++) READI;
        nQUAD++;
        break;

      case 11: // 10-node TETRA (second order)
        READI; (*cnTETRA)(nQUAD,1) = indirNodes[ti-1];
        READI; (*cnTETRA)(nQUAD,2) = indirNodes[ti-1];
        READI; (*cnTETRA)(nQUAD,3) = indirNodes[ti-1];
        READI; (*cnTETRA)(nQUAD,4) = indirNodes[ti-1];
        nTETRA++;
        for (E_Int j = 4; j < 10; j++) READI;
        break;

      case 12: // 27-node HEXA (second order)
        READI; (*cnHEXA)(nHEXA,1) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,2) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,3) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,4) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,5) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,6) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,7) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,8) = indirNodes[ti-1];
        for (E_Int j = 8; j < 27; j++) READI; 
        nHEXA++;
        break;

      case 13: // 18-node PENTA (second order)
        READI; (*cnPENTA)(nPENTA,1) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,2) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,3) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,4) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,5) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,6) = indirNodes[ti-1];
        for (E_Int j = 6; j < 18; j++) READI; 
        nPENTA++;
        break;

      case 14: // 14-node PYRA (second order)
        READI; (*cnPYRA)(nPYRA,1) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,2) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,3) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,4) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,5) = indirNodes[ti-1];
        for (E_Int j = 5; j < 14; j++) READI; 
        nPYRA++;
        break;

      case 15: // NODE
        READI; (*indNODE)[nNODE] = ti;
        nNODE++;
        break;

      case 16: // 8-node QUAD (second order)
        READI; (*cnQUAD)(nQUAD,1) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,2) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,3) = indirNodes[ti-1];
        READI; (*cnQUAD)(nQUAD,4) = indirNodes[ti-1];
        for (E_Int j = 4; j < 8; j++) READI; 
        nQUAD++;
        break;

      case 17: // 20-node HEXA (second order)
        READI; (*cnHEXA)(nHEXA,1) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,2) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,3) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,4) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,5) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,6) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,7) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,8) = indirNodes[ti-1];
        for (E_Int j = 8; j < 20; j++) READI; 
        nHEXA++;
        break;

      case 18: // 15-node PENTA (second order)
        READI; (*cnPENTA)(nPENTA,1) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,2) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,3) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,4) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,5) = indirNodes[ti-1];
        READI; (*cnPENTA)(nPENTA,6) = indirNodes[ti-1];
        for (E_Int j = 6; j < 15; j++) READI; 
        nPENTA++;
        break;

      case 19: // 13-node PYRA (second order)
        READI; (*cnPYRA)(nPYRA,1) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,2) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,3) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,4) = indirNodes[ti-1];
        READI; (*cnPYRA)(nPYRA,5) = indirNodes[ti-1];
        for (E_Int j = 5; j < 13; j++) READI;
        nPYRA++;
        break;

      case 20: // 9-node TRI (third order)
        READI; (*cnTRI)(nTRI,1) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,2) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,3) = indirNodes[ti-1];
        for (E_Int j = 3; j < 9; j++) READI; 
        nTRI++;
        break;
        
      case 21: // 10-node TRI (third order)
        READI; (*cnTRI)(nTRI,1) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,2) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,3) = indirNodes[ti-1];
        for (E_Int j = 3; j < 10; j++) READI; 
        nTRI++;
        break;

      case 22: // 12-node TRI (fourth order)
        READI; (*cnTRI)(nTRI,1) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,2) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,3) = indirNodes[ti-1];
        for (E_Int j = 3; j < 12; j++) READI; 
        nTRI++;
        break;

      case 23: // 15-node TRI (fourth order)
        READI; (*cnTRI)(nTRI,1) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,2) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,3) = indirNodes[ti-1];
        for (E_Int j = 3; j < 15; j++) READI; 
        nTRI++;
        break;

      case 24: // 15-node TRI (fifth order)
        READI; (*cnTRI)(nTRI,1) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,2) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,3) = indirNodes[ti-1];
        for (E_Int j = 3; j < 15; j++) READI; 
        nTRI++;
        break;

      case 25: // 21-node TRI (fifth order)
        READI; (*cnTRI)(nTRI,1) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,2) = indirNodes[ti-1];
        READI; (*cnTRI)(nTRI,3) = indirNodes[ti-1];
        for (E_Int j = 3; j < 21; j++) READI; 
        nTRI++;
        break;

      case 26: // 4-node BAR (third order)
        READI; (*cnBAR)(nBAR,1) = indirNodes[ti-1];
        READI; (*cnBAR)(nBAR,2) = indirNodes[ti-1];
        for (E_Int j = 2; j < 4; j++) READI; 
        nBAR++;
        break;

      case 27: // 5-node BAR (fourth order)
        READI; (*cnBAR)(nBAR,1) = indirNodes[ti-1];
        READI; (*cnBAR)(nBAR,2) = indirNodes[ti-1];
        for (E_Int j = 2; j < 5; j++) READI; 
        nBAR++;
        break;

      case 28: // 6-node BAR (fifth order)
        READI; (*cnBAR)(nBAR,1) = indirNodes[ti-1];
        READI; (*cnBAR)(nBAR,2) = indirNodes[ti-1];
        for (E_Int j = 2; j < 6; j++) READI; 
        nBAR++;
        break;

      case 29: // 20-node TETRA (third order)
        READI; (*cnTETRA)(nTETRA,1) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,2) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,3) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,4) = indirNodes[ti-1];
        for (E_Int j = 4; j < 20; j++) READI; 
        nTETRA++;
        break;
        
      case 30: // 35-node TETRA (fourth order)
        READI; (*cnTETRA)(nTETRA,1) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,2) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,3) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,4) = indirNodes[ti-1];
        for (E_Int j = 4; j < 35; j++) READI; 
        nTETRA++;
        break;

      case 31: // 56-node TETRA (fifth order)
        READI; (*cnTETRA)(nTETRA,1) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,2) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,3) = indirNodes[ti-1];
        READI; (*cnTETRA)(nTETRA,4) = indirNodes[ti-1];
        for (E_Int j = 4; j < 56; j++) READI; 
        nTETRA++;
        break;
        
      case 92: // 64-node HEXA (third order)
        READI; (*cnHEXA)(nHEXA,1) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,2) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,3) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,4) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,5) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,6) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,7) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,8) = indirNodes[ti-1];
        for (E_Int j = 8; j < 64; j++) READI; 
        nHEXA++;
        break;

      case 93: // 125-node HEXA (fourth order)
        READI; (*cnHEXA)(nHEXA,1) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,2) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,3) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,4) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,5) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,6) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,7) = indirNodes[ti-1];
        READI; (*cnHEXA)(nHEXA,8) = indirNodes[ti-1];
        for (E_Int j = 8; j < 125; j++) READI; 
        nHEXA++;
        break;

      default:
        nDiscard++;
    }
  }
