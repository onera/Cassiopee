/*    
    Copyright 2013-2025 Onera.

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

// convertit un maillage MIX en plusieurs maillages BE

# include "converter.h" 
# include "kcore.h"
# include <string.h>
# include <stdio.h>

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Conversion du numpy MIX en list de numpys BE :
  BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA (None si aucun elt)
 */
//=============================================================================
PyObject* K_CONVERTER::convertMix2BE(PyObject* self, PyObject* args)
{
  PyObject* oMIX;
  if (!PYPARSETUPLE_(args, O_, &oMIX)) return NULL;
 
  // Check numpy
  FldArrayI* cMIX;
  E_Int res = K_NUMPY::getFromNumpyArray(oMIX, cMIX, true);

  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertMix2BE: numpy is invalid.");
    return NULL;
  }
  FldArrayI cBAR; FldArrayI cTRI; FldArrayI cQUAD;
  FldArrayI cTETRA; FldArrayI cPYRA;
  FldArrayI cPENTA; FldArrayI cHEXA;
  K_CONNECT::connectMix2EV(*cMIX,
                           cBAR, cTRI, cQUAD, cTETRA,
                           cPYRA, cPENTA, cHEXA);

  RELEASESHAREDN(oMIX, cMIX);

  PyObject* oBAR;
  if (cBAR.getSize() > 0) oBAR = K_NUMPY::buildNumpyArray(cBAR, false);
  else { Py_INCREF(Py_None); oBAR = Py_None; }
  PyObject* oTRI;
  if (cTRI.getSize() > 0) oTRI = K_NUMPY::buildNumpyArray(cTRI, false);
  else { Py_INCREF(Py_None); oTRI = Py_None; }
  PyObject* oQUAD;
  if (cQUAD.getSize() > 0) oQUAD = K_NUMPY::buildNumpyArray(cQUAD, false);
  else { Py_INCREF(Py_None); oQUAD = Py_None; }
  PyObject* oTETRA;
  if (cTETRA.getSize() > 0) oTETRA = K_NUMPY::buildNumpyArray(cTETRA, false);
  else { Py_INCREF(Py_None); oTETRA = Py_None; }
  PyObject* oPYRA;
  if (cPYRA.getSize() > 0) oPYRA = K_NUMPY::buildNumpyArray(cPYRA, false);
  else { Py_INCREF(Py_None); oPYRA = Py_None; }
  PyObject* oPENTA;
  if (cPENTA.getSize() > 0) oPENTA = K_NUMPY::buildNumpyArray(cPENTA, false);
  else { Py_INCREF(Py_None); oPENTA = Py_None; }
  PyObject* oHEXA;
  if (cHEXA.getSize() > 0) oHEXA = K_NUMPY::buildNumpyArray(cHEXA, false);
  else { Py_INCREF(Py_None); oHEXA = Py_None; }

  PyObject* tpl = Py_BuildValue("[OOOOOOO]", oBAR, oTRI, oQUAD,
                                oTETRA, oPYRA, oPENTA, oHEXA);
  return tpl;
}
