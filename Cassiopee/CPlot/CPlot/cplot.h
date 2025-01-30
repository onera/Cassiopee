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

#ifndef _CPLOT_CPLOT_H_
#define _CPLOT_CPLOT_H_

#define KFAILED 0
#define KSUCCESS 1

#include "kcore.h"
#include "Zone.h"
#include <vector>

// Si USEDLMESH=1, on cree des DL pour le maillage
// sinon, on utilise la DLsolid en mode GL_LINE
#define USEDLMESH 0

namespace K_CPLOT
{ 
PyObject* render(PyObject* self, PyObject* args);
PyObject* deletez(PyObject* self, PyObject* args);
PyObject* add(PyObject* self, PyObject* args);
PyObject* replace(PyObject* self, PyObject* args);
PyObject* pressKey(PyObject* self, PyObject* args);
PyObject* displayNew(PyObject* self, PyObject* args);
PyObject* displayAgain(PyObject* self, PyObject* args);
PyObject* setFileName(PyObject* self, PyObject* args);
PyObject* getState(PyObject* self, PyObject* args);
PyObject* getSelectedZone(PyObject* self, PyObject* args);
PyObject* getSelectedZones(PyObject* self, PyObject* args);
PyObject* getSelectedStatus(PyObject* self, PyObject* args);
PyObject* getActiveZones(PyObject* self, PyObject* args);
PyObject* getActiveStatus(PyObject* self, PyObject* args);
PyObject* getActivePoint(PyObject* self, PyObject* args);
PyObject* getActivePointIndex(PyObject* self, PyObject* args);
PyObject* getActivePointF(PyObject* self, PyObject* args);
PyObject* getMouseState(PyObject* self, PyObject* args);
PyObject* getKeyboard(PyObject* self, PyObject* args);
PyObject* resetKeyboard(PyObject* self, PyObject* args);
PyObject* setState(PyObject* self, PyObject* args);
PyObject* setShaderPath(PyObject* self, PyObject* args);
PyObject* setWindowTitle(PyObject* self, PyObject* args);
PyObject* setMode(PyObject* self, PyObject* args);
PyObject* setActivePoint(PyObject* self, PyObject* args);
PyObject* changeVariable(PyObject* self, PyObject* args);
PyObject* changeStyle(PyObject* self, PyObject* args);
PyObject* changeInfoDisplay(PyObject* self, PyObject* args);
PyObject* changeBlanking(PyObject* self, PyObject* args);
PyObject* setDim(PyObject* self, PyObject* args);
PyObject* setSelectedZones(PyObject* self, PyObject* args);
PyObject* unselectAllZones(PyObject* self, PyObject* args);
PyObject* setActiveZones(PyObject* self, PyObject* args);
PyObject* setZoneNames(PyObject* self, PyObject* args);
PyObject* lookFor(PyObject* self, PyObject* args);
PyObject* fitView(PyObject* self, PyObject* args);
PyObject* isDisplayRunning(PyObject* self, PyObject* args);
PyObject* finalizeExport(PyObject* self, PyObject* args);
PyObject* hide(PyObject* self, PyObject* args);
PyObject* show(PyObject* self, PyObject* args);
PyObject* display1D(PyObject* self, PyObject* args);
PyObject* configure(PyObject* self, PyObject* args);
PyObject* panorama(PyObject* self, PyObject* args);
PyObject* panoramaODS(PyObject* self, PyObject* args);
PyObject* blur(PyObject* self, PyObject* args);
}

E_Int getMode(PyObject* modeObject);
E_Int getScalarField(PyObject* scalarFieldMode);
void findMinMax(Zone* zone);
void findFMinMax(Zone* zone);
void globMinMax(Zone** zones, E_Int nz,
                double& xmin, double& xmax,
                double& ymin, double& ymax,
                double& zmin, double& zmax,
                double& epsup, double& epsstrafe, double& dmoy);
void globFMinMax(Zone** zones, E_Int nz, double* minf, double* maxf);
int getStringsFromPyObj(PyObject* obj, std::vector<char*>& strings);
int getStringFromPyObj(PyObject* obj, char*& string);
void insertAfterNz(Zone** zonesp, E_Int& lzonesn, Zone**& zonesn, E_Int nz, Zone* z);
void deleteNz(Zone** zonesp, E_Int& lzonesn, Zone**& zonesn, E_Int nz);
#endif
