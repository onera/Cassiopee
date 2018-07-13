/*    
    Copyright 2013 Onera.

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

# include "connector.h"
# include "stub.h"

//=============================================================================
/* Calcul et stocke les coefficients d'interpolation pour les centres des 
   cellules
   IN: (Nir,Njr): dimensions (i,j) des blocs interpoles (receveurs)
   IN: coordArrays: pts a interpoler
   IN: interpArrays: domaines d interpolation en noeuds
   IN: interpCellN: cellN en centres pour les domaines d interpolation
   IN: order: order des interpolations
   IN: isEX: indique si on est en depth=2 (isEX = 0) ou en depth=1 (isEX=1)
   OUT: indices des points interpoles et des cellules donneuses et coefficients d'interpolation 
         ( + liste des directions pour les points EX quand isEX=1)
*/
//=============================================================================
PyObject* K_CONNECTOR::setInterpolations(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

