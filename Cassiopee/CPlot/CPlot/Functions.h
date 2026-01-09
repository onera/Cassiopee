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

#ifndef _CPLOT_FUNCTIONS_H
#define _CPLOT_FUNCTIONS_H

// - Math -
#define ABS(x) (x>0 ? x : -x)
#define MIN(x,y) (x>y ? y : x)
#define MAX(x,y) (x>y ? x : y)

// - Pour les NGONS -
// Nbre de faces
#define NFACES(connect) connect[0]
// Nbre d'elements
#define NELTS(connect) connect[2+connect[1]]
// Pointeur sur le debut des faces
#define PTRFACES(connect) &connect[2]
// Position du debut des faces dans connect
#define POSFACES(connect) 2
// Pointeur sur le debut des elements
#define PTRELTS(connect) &connect[4+connect[1]]
// Position du debut des elements
#define POSELTS(connect) (4+connect[1])
#endif
