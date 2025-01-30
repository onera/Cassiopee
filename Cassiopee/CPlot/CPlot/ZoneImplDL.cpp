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
#if defined (_WIN32) ||defined(_WIN64)
#  include <windows.h>
#endif
#if defined(__APPLE__)
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

# include "CPlotStateDL.h"
# include "ZoneImplDL.h"

ZoneImplDL::ZoneImplDL() : ZoneImpl(),
			   _DLmesh(0), _DLsolid(0), _DLiso(0),
			   _DLisoField(-1), _DLisoField2(-1), _DLisoField3(-1)
{}
// ========================================================================
void ZoneImplDL::freeGPURes( CPlotState* pt_state, bool freeIso )
{
  CPlotStateDL& state = *static_cast<CPlotStateDL*>(pt_state);
  if (_DLmesh != 0) state.freeDLList.push_back(_DLmesh);
  if (_DLsolid != 0) state.freeDLList.push_back(_DLsolid);
  if (freeIso) {
    if (_DLiso != 0) state.freeDLList.push_back(_DLiso);
    _DLiso = 0;
  }
  _GPUResUse = 1; _DLmesh = 0; _DLsolid = 0;
  state.freeGPURes = 1;
}
// ========================================================================
void
ZoneImplDL::destroyIsoField()
{
  if (_DLiso != 0) { 
    glDeleteLists(_DLiso, 1); 
    _DLisoField = -1; _DLiso = 0; 
  }
}
