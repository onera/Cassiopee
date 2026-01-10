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

#ifndef _CPLOT_ZONE1D_H_
#define _CPLOT_ZONE1D_H_
# include "kcore.h"

/* Define a set of data for 1D display 
   Correspond a un array */
class Zone1D
{
public:
  Zone1D(char* varString, K_FLD::FldArrayF& f, K_FLD::FldArrayI& cn);
  ~Zone1D();

public:
  E_Int _np; // nbre de pts dans f
  E_Int _nv; // nombre de variables dans f
  E_Int _ne; // nbre d'elements dans la BAR
  double* _f; // champ des variables
  E_Int* _cn; // connectivite
  char** _varNames; // nom des variables
};

#endif
