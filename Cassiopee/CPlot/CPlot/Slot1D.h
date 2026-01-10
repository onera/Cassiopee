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

#ifndef _CPLOT_SLOT1D_H_
#define _CPLOT_SLOT1D_H_

#include "Zone1D.h"
#include <vector>

#ifndef MAXSTRINGLENGTH
#define MAXSTRINGLENGTH 128   /* max length of strings */
#endif

/* Define 1D display slot */
class Slot1D
{
public:
  Slot1D(int no, int gridPosI, int gridPosJ, double bgBlend);
  ~Slot1D();
  
  E_Int _no;          // le no du slot
  E_Int _gridPosI;    // position de ce slot dans la grille
  E_Int _gridPosJ;    
  double _bgBlend;  // le % de blending sur le fond (0-1)
  char _title[MAXSTRINGLENGTH]; // titre du slot
  double _r1min;
  double _r1max;
  double _r2min;
  double _r2max;
  std::vector<Zone1D*> _zones; // data a afficher dans ce slot 
  std::vector<E_Int> _var1; // var1 pos
  std::vector<char> _var1NameC1; // char 1 of var1
  std::vector<char> _var1NameC2; // char 2 of var1
  std::vector<E_Int> _var2; // var2 pos
  std::vector<char> _var2NameC1; // char 1 of var2
  std::vector<char> _var2NameC2; // char 2 of var2  
};

#endif
