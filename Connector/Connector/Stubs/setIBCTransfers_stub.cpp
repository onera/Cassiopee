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
# include "connector.h"
#include "stub.h"

using namespace std;

//=============================================================================
//Retourne -2: incoherence entre meshtype et le type d interpolation
//         -1: type invalide
//          1: ok
// Entree/Sortie: des variables conservatives ( + ronutildeSA ) 
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar1(
  E_Int bctype,
  E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
  E_Float* xPC, E_Float* yPC, E_Float* zPC,
  E_Float* xPW, E_Float* yPW, E_Float* zPW,
  E_Float* xPI, E_Float* yPI, E_Float* zPI, 
  E_Float* densPtr, E_Float* pressPtr,
  E_Float* vxPtr, E_Float* vyPtr, E_Float* vzPtr,
  E_Float* utauPtr, E_Float* yplusPtr,
  E_Float* d1, E_Float* d2, E_Float* d3, E_Float* d4, E_Float* d5,
  E_Float* tmp, E_Int& size,
  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts, E_Float Pr,
  vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields)
{
  printf("%s\n", STUBMSG);

  return 1;
}
//=============================================================================
//Retourne -2: incoherence entre meshtype et le type d'interpolation
//         -1: type invalide
//          1: ok
// Entree/Sortie:  (ro,u,v,w,t) ( + nutildeSA ) 
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar2(
  E_Int bctype,
  E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
  E_Float* xPC, E_Float* yPC, E_Float* zPC,
  E_Float* xPW, E_Float* yPW, E_Float* zPW,
  E_Float* xPI, E_Float* yPI, E_Float* zPI, 
  E_Float* densPtr, E_Float* pressPtr, 
  E_Float* vxPtr, E_Float* vyPtr, E_Float* vzPtr, 
  E_Float* utauPtr, E_Float* yplusPtr,
  E_Float* d1, E_Float* d2, E_Float* d3, E_Float* d4, E_Float* d5,
  E_Float* tmp, E_Int& size,
  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts, E_Float Pr,
  vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields,
  E_Int nbptslinelets, E_Float* linelets, E_Int* indexlinelets)
{
  printf("%s\n", STUBMSG);
  return 1;
}
//=============================================================================
//Retourne -2 : incoherence entre meshtype et le type d interpolation
//         -1 : type invalide
//          1 : ok
// Entree/Sortie :  (ro,u,v,w,p) ( + ronutildeSA ) 
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar3(
  E_Int bctype,
  E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin,  E_Int& ithread,
  E_Float* xPC, E_Float* yPC, E_Float* zPC,
  E_Float* xPW, E_Float* yPW, E_Float* zPW,
  E_Float* xPI, E_Float* yPI, E_Float* zPI, 
  E_Float* densPtr, E_Float* pressPtr, 
  E_Float* vxPtr, E_Float* vyPtr, E_Float* vzPtr, 
  E_Float* utauPtr, E_Float* yplusPtr,
  E_Float* d1, E_Float* d2, E_Float* d3, E_Float* d4, E_Float* d5,
  E_Float* tmp, E_Int& size,
  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts, E_Float Pr,
  vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields)
{
  printf("%s\n", STUBMSG);
  return 1;
}
//=============================================================================
/* Effectue les transfers IBC */
//=============================================================================
//Stephanie: fonction absente du pytree?
PyObject* K_CONNECTOR::setIBCTransfers(PyObject* self, PyObject* args)
{
  printf("%s\n", STUBMSG);
  return NULL;
}
//=============================================================================
/* Effectue les transfers IBC en in-place */
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers(PyObject* self, PyObject* args)
{
  printf("%s\n", STUBMSG);
  return Py_None;
}
