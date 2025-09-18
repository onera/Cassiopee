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

# include "kcore.h"
# include "cplot.h"
# include "Data.h"
#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

using namespace K_FLD;
using namespace std;

//=============================================================================
/* display arrays */
//=============================================================================
PyObject* K_CPLOT::display1D(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Int slot;
  E_Int gridPosI, gridPosJ;
  E_Int gridSizeI, gridSizeJ;
  E_Float bgBlend;
  char* var1; char* var2;
  E_Float r1min, r1max, r2min, r2max;
    
  if (!PYPARSETUPLE_(args, "O_ I_ TII_ TII_ R_ SS_ TRR_ TRR_", 
                        &arrays, &slot, 
                        &gridPosI, &gridPosJ, &gridSizeI, &gridSizeJ,
                        &bgBlend, &var1, &var2,
                        &r1min, &r1max, &r2min, &r2max))
  {
    return NULL;
  }

  // Recuperation du container de donnees
  Data* d = Data::getInstance();

  // Lecture des arrays
  vector<E_Int> res;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = false;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "display1D: invalid list of arrays.");
    E_Int structFSize = structF.size();
    for (E_Int i = 0; i < structFSize; i++)
      RELEASESHAREDS(objs[i], structF[i]);

    E_Int unstrFSize = unstrF.size();
    for (E_Int i = 0; i < unstrFSize; i++)
      RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
    return NULL;
  }
  E_Int unstrFSize = unstrF.size();

  // Grid size
  if (gridSizeI != -1) d->ptrState->gridSizeI = gridSizeI;
  if (gridSizeJ != -1) d->ptrState->gridSizeJ = gridSizeJ;

  // le slot existe-t-il deja?
  E_Int mySlot = slot;
  Slot1D* s;
  vector<Slot1D*>& slots = d->_slots1D;
  size_t ns = slots.size();
  E_Int found = -1;
  for (size_t i = 0; i < ns; i++)
  {
    s = slots[i];
    if (s->_no == mySlot) { found = i; break; } 
  }

  // Switch - Dangerous zone protegee par ptrState->lock
  d->ptrState->syncDisplay();
  if (found != -1) // remplace le slot par celui-ci 
  {
    s = slots[found]; delete s;
    if (unstrFSize != 0) // replace slot
    {
      slots[found] = new Slot1D(mySlot, gridPosI, gridPosJ, bgBlend); 
      s = slots[found];
    }
    else // only delete slot
    {
      vector<Slot1D*>::iterator it = slots.begin()+found;
      slots.erase(it);
    }
  }
  else
  {
    if (unstrFSize != 0)
    {
      s = new Slot1D(mySlot, gridPosI, gridPosJ, bgBlend);
      slots.push_back(s);
    }
  }

  // If 1D plot is linked to view: force r1min, r1max, r2min, r2max
  E_Int linkView = 0;

  // Get data
  E_Int posvar1; E_Int posvar2;
  char c1, c2;
  for (E_Int i = 0; i < unstrFSize; i++)
  {
    Zone1D* z = new Zone1D(unstrVarString[i], *unstrF[i], *cnt[i]);
    s->_zones.push_back(z);
    // Find posvars
    posvar1 = K_ARRAY::isNamePresent(var1, unstrVarString[i]);
    if (posvar1 == -1) posvar1 = 0;
    s->_var1.push_back((int)posvar1);
    d->getCharsFromVarName(var1, c1, c2);
    s->_var1NameC1.push_back(c1);
    s->_var1NameC2.push_back(c2);

    posvar2 = K_ARRAY::isNamePresent(var2, unstrVarString[i]);
    if (posvar2 == -1) posvar2 = 0;
    s->_var2.push_back((int)posvar2);
    d->getCharsFromVarName(var2, c1, c2);
    s->_var2NameC1.push_back(c1);
    s->_var2NameC2.push_back(c2);

    // ranges
    if (linkView == 0)
    {
      s->_r1min = r1min;
      s->_r1max = r1max;
      s->_r2min = r2min;
      s->_r2max = r2max;
    }
    else
    {
      //printf("Init: %f %f %f %f\n", r1min, r1max, r2min, r2max);
      d->link2View(z, posvar1, posvar2, r1min, r1max, r2min, r2max);
      //printf("New: %f %f %f %f\n", r1min, r1max, r2min, r2max);
      s->_r1min = r1min;
      s->_r1max = r1max;
      s->_r2min = r2min;
      s->_r2max = r2max;
    }
  }

  // Here we go
  d->ptrState->render = 1;

  // Free the input arrays
  E_Int structFSize = structF.size();
  for (E_Int i = 0; i < structFSize; i++) 
    RELEASESHAREDS(objs[i], structF[i]);
  for (E_Int i = 0; i < unstrFSize; i++)
    RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);

  // Retourne le hook
  return Py_BuildValue("l", d);
}
