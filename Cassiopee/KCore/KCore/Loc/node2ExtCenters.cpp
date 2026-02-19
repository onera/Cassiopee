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

# include "loc.h"

using namespace std;
using namespace K_FLD;
//=============================================================================
/* Passage des données de noeuds aux centres étendus 
   Retourne 1 en cas de succès, 0 en cas d'échec*/
//=============================================================================
E_Int K_LOC::node2ExtCenterStruct(E_Int imo, E_Int jmo, E_Int kmo,
                                  FldArrayF& FNode,
                                  E_Int ime, E_Int jme, E_Int kme, 
                                  FldArrayF& FExtCenter)
{
  E_Int i, j, k;
  E_Int pos, pos2, pos2a, pos2b, pos2c, pos2d, pos2e, pos2f, pos2g;
  E_Int nfld = FNode.getNfld();
  E_Int size = ime*jme*kme;

  if (FExtCenter.getSize() != size || FExtCenter.getNfld() != nfld) 
    FExtCenter.malloc(size, nfld);
  FExtCenter.setAllValuesAtNull();
  E_Int imjme = ime*jme;
  E_Int imjmo = imo*jmo;
  for (E_Int v = 1; v <= nfld; v++)
  {
    E_Float* fn = FNode.begin(v);
    E_Float* fe = FExtCenter.begin(v);

    /* corners */
    pos  = 0; pos2 = 0;
    fe[pos] = fn[pos2];

    pos  = ime-1; pos2 = imo-1;
    fe[pos] = fn[pos2];
    
    pos  = (jme-1)*ime; pos2 = (jmo-1)*imo;
    fe[pos] = fn[pos2];

    pos = (kme-1)*imjme; pos2 = (kmo-1)*imjmo;
    fe[pos] = fn[pos2];

    pos  = ime-1 + (jme-1)*ime; pos2 = imo-1 + (jmo-1)*imo;
    fe[pos] = fn[pos2];

    pos = ime-1 + (kme-1)*imjme; pos2 = imo-1 + (kmo-1)*imjmo;
    fe[pos] = fn[pos2];

    pos = (jme-1)*ime + (kme-1)*imjme; pos2 = (jmo-1)*imo + (kmo-1)*imjmo;
    fe[pos] = fn[pos2];
  
    pos = ime-1 + (jme-1)*ime + (kme-1)*imjme; pos2 = imo-1 + (jmo-1)*imo + (kmo-1)*imjmo;
    fe[pos] = fn[pos2];
  
    /* border lines */
    for (i = 2; i < ime; i++)
    {
      pos   = i-1; pos2  = i-1; pos2a = i-2;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);

      pos   = i-1 + (jme-1)*ime;
      pos2  = i-1 + (jmo-1)*imo;
      pos2a = i-2 + (jmo-1)*imo;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);

      pos   = i-1 + (kme-1)*imjme;
      pos2  = i-1 + (kmo-1)*imjmo;
      pos2a = pos2 -1 ;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);
 
      pos   = i-1 + (jme-1)*ime + (kme-1)*imjme;
      pos2  = i-1 + (jmo-1)*imo + (kmo-1)*imjmo;
      pos2a = pos2-1;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);  
    }
  
    for (j = 2; j < jme; j++)
    {
      pos   = (j-1)*ime;
      pos2  = (j-1)*imo; 
      pos2a = pos2 - imo;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);
 
      pos   = ime-1 + (j-1)*ime;
      pos2  = imo-1 + (j-1)*imo; 
      pos2a = pos2 - imo;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);
    
      pos   = (j-1)*ime + (kme-1)*imjme;
      pos2  = (j-1)*imo + (kmo-1)*imjmo;
      pos2a = pos2 - imo;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);
 
      pos   = (ime-1) + (j-1)*ime + (kme-1)*imjme;
      pos2  = (imo-1) + (j-1)*imo + (kmo-1)*imjmo;
      pos2a = pos2 - imo;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);  
    }

    for (k = 2; k < kme; k++)
    {
      pos   = (k-1)*imjme; 
      pos2  = (k-1)*imjmo;
      pos2a = pos2 - imjmo;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);
     
      pos   = (ime-1) + (k-1) * imjme;
      pos2  = (imo-1) + (k-1) * imjmo;
      pos2a = pos2 - imjmo;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);
 
      pos   = (jme-1)*ime + (k-1)*imjme;
      pos2  = (jmo-1)*imo + (k-1)*imjmo;
      pos2a = pos2 - imjmo;;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);
          
      pos   = ime-1 + (jme-1)*ime + (k-1)*imjme;
      pos2  = imo-1 + (jmo-1)*imo + (k-1)*imjmo;
      pos2a = pos2 - imjmo;
      fe[pos] = K_CONST::ONE_HALF*(fn[pos2] + fn[pos2a]);
    }
    
    /* inside plane */
    for (j = 2; j < jme; j++)
      for (i = 2; i < ime; i++)
      {
        pos   = i-1+(j-1)*ime;
        pos2  = i-1+(j-1)*imo;
        pos2a = pos2-1;
        pos2b = pos2-imo;
        pos2c = pos2-1-imo;
        fe[pos] = K_CONST::ONE_FOURTH*(fn[pos2]+fn[pos2a]+fn[pos2b]+fn[pos2c]);
  
        pos   = i-1+(j-1)*ime + (kme-1)*imjme;
        pos2  = i-1+(j-1)*imo + (kmo-1)*imjmo;
        pos2a = pos2-1; pos2b = pos2-imo; pos2c = pos2-1-imo;
        fe[pos] = K_CONST::ONE_FOURTH*(fn[pos2]+fn[pos2a]+fn[pos2b]+fn[pos2c]);  
      }
    
    for (k = 2; k < kme; k++)
      for (i = 2; i < ime; i++)
      {
        pos   = i-1 + (k-1)*imjme;
        pos2  = i-1 + (k-1)*imjmo;
        pos2a = pos2-1;
        pos2b = pos2-imjmo;
        pos2c = pos2-1-imjmo;
        fe[pos] = K_CONST::ONE_FOURTH*(fn[pos2]+fn[pos2a]+fn[pos2b]+fn[pos2c]);
     
        pos   = (i-1) + (jme-1)*ime + (k-1)*imjme;
        pos2  = (i-1) + (jmo-1)*imo + (k-1)*imjmo;
        pos2a = pos2-1;
        pos2b = pos2-imjmo;
        pos2c = pos2-1-imjmo;
        fe[pos] = K_CONST::ONE_FOURTH*(fn[pos2]+fn[pos2a]+fn[pos2b]+fn[pos2c]); 
      }

    for (k = 2; k < kme; k++)
      for (j = 2; j < jme; j++)
      {
        pos   = (j-1)*ime + (k-1)*imjme;
        pos2  = (j-1)*imo + (k-1)*imjmo;
        pos2a = pos2-imjmo;
        pos2b = pos2-imo;
        pos2c = pos2-imo-imjmo;
        fe[pos] = K_CONST::ONE_FOURTH*(fn[pos2]+fn[pos2a]+fn[pos2b]+fn[pos2c]);

        pos   = (ime-1) + (j-1)*ime + (k-1)*imjme;
        pos2  = (imo-1) + (j-1)*imo + (k-1)*imjmo;
        pos2a = pos2-imjmo;
        pos2b = pos2-imo;
        pos2c = pos2-imo-imjmo;
        fe[pos] = K_CONST::ONE_FOURTH*(fn[pos2]+fn[pos2a]+fn[pos2b]+fn[pos2c]);
      }

    /* Points interieurs */
    for (k = 2; k < kme; k++)
      for (j = 2; j < jme; j++)
        for (i = 2; i < ime; i++)
        {
          pos   = i-1 + (j-1)*ime + (k-1)*imjme;
          pos2  = i-1 + (j-1)*imo + (k-1)*imjmo;
          pos2a = pos2-1;
          pos2b = pos2-imo;
          pos2c = pos2-imjmo;
          pos2d = pos2-1-imo;
          pos2e = pos2-1-imjmo;
          pos2f = pos2-1-imo-imjmo;
          pos2g = pos2-imo-imjmo;
          fe[pos] = K_CONST::ONE_EIGHTH*(fn[pos2] + fn[pos2a]+
                                        fn[pos2b] + fn[pos2c]+
                                        fn[pos2d] + fn[pos2e]+
                                        fn[pos2f] + fn[pos2g]);
  
        }
  }
  return 1;
}
