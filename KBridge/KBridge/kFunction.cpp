/*    
    Copyright 2013-2023 Onera.

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

// motion defined by a Cassiopee DesFunction

#include "kbridge.h"

# define E_DOUBLEREAL
//# define E_DOUBLEINT 
# define _DEF_USE_ISO_
# define _DEF_USE_USING_
# define __STD std
# define _ISO_LIST_
# define _ISO_VECTOR_
# define _E_USE_STANDARD_IOSTREAM_
# define _DEF_USE_ISO_

#include <vector>
#include <list>
#include <map>
using namespace std;

# include "Descp/Base/DesApi.h"
# include "Fact/Base/FactBase.h"
# include "Fact/Func/FactMotionRotorP.h"
# include "Fact/Func/FactMotionWingP.h"
# include "Fact/Func/FactMotionTransP.h"
# include "Tbx/Clock/TbxClock.h"

//=============================================================================
PyObject* K_KBRIDGE::evalKDesFunction(PyObject* self, PyObject* args)
{
  E_Float t;
  PyObject* functionName;
  if (!PyArg_ParseTuple(args, "Od", &functionName, &t)) return NULL;

  DesFunction* desFunction;
  if (PyObject_HasAttrString(functionName, "this") == 1) 
  {
    char* s = PyString_AsString(PyObject_GetAttrString(functionName, "this"));
    E_Int ret = DesApi::getInstance()->hackFunction(s, desFunction);
  
    if (ret == 1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "evalKDesFunction: desFunction is not a valid Cartesian Solver DesFunction.");
      return NULL;
    }
  }
  
  TbxMotionBase* mr;
  FldArrayF s0(3);
  FldArrayF omega(3);
  FldArrayF pntAxis(3);
  E_Float tsav;
      
  // Get motion
  TbxString typeFunc = FactBase::getDefS(KEY_FUNCTYPE, 
                                         (const DesBase&) *desFunction);
  if (typeFunc == E_MOTIONROTOR)
  {
    FactMotionRotorP fact;
    mr = (TbxMotionBase*)(fact.createFunc(*desFunction));
  }
  else if (typeFunc == E_MOTIONWING)
  {
    FactMotionWingP fact;
    mr = (TbxMotionBase*)(fact.createFunc(*desFunction));
  }
  else if(typeFunc == E_MOTIONTRANS)
  {
    FactMotionTransP fact;
    mr = (TbxMotionBase*)(fact.createFunc(*desFunction));
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                     "evalKDesFunction: desFunction is not a known motion.");
    return NULL;
  }

  // Compute motion
  tsav = TbxClock::instance()->getTime();
  TbxClock::instance()->setTime(t);
  mr->computeTimeMotion(s0, omega, pntAxis);
  
  // Get matrices
  const FldArrayF& mt  = *(mr->getRotMat()); //matrice de rotation
  const FldArrayF& r0t = *(mr->getR0()); // deplacement
  const FldArrayF& x0t = *(mr->getX0()); //centre de la rotation
  
  TbxClock::instance()->setTime(tsav);

  PyObject* tpl = Py_BuildValue("([ddd],[ddd],[[ddd],[ddd],[ddd]])",
                                r0t[0], r0t[1], r0t[2],
                                x0t[0], x0t[1], x0t[2],
                                mt(0,1), mt(0,2), mt(0,3),
                                mt(1,1), mt(1,2), mt(1,3),
                                mt(2,1), mt(2,2), mt(2,3)); 
  return tpl;
}
