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
// read cpacs, export CAD (use tigl)

#include "tigl.h"
#include "tixi.h"
#include "modeler.h"

PyObject* K_MODELER::exportCAD(PyObject* self, PyObject* args)
{
  char* cpacsFileName; char* exportFileName; char* format;
  if (!PYPARSETUPLE_(args, SSS_, &cpacsFileName, &exportFileName, &format)) return NULL;

  TixiDocumentHandle tixiHandle;
  TiglCPACSConfigurationHandle tiglHandle;
  int nWings = 0;
  int nFuselages = 0;
  int i = 0;
  double x, y, z;

  // read in cpacs xml file
  if (tixiOpenDocument(cpacsFileName, &tixiHandle) != SUCCESS)
  {
    return NULL;
  }

  // enable logging of errors and warnings into file
  tiglLogSetFileEnding("txt");
  tiglLogSetTimeInFilenameEnabled(TIGL_FALSE);
  tiglLogToFileEnabled("demolog");

  // now open cpacs file with tigl
  if (tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle) != TIGL_SUCCESS)
  {
    printf("Error: exportCAD: reading in cpacs file with TiGL.\n");
    return NULL;
  }

  // Check cpacs validity
  TiglBoolean isValid;
  tiglIsCPACSConfigurationHandleValid(tiglHandle, &isValid);
  printf("Validation of cpacs: %d\n", isValid);

  // query number of wings and fuselages and their names
  tiglGetWingCount(tiglHandle, &nWings);
  tiglGetFuselageCount(tiglHandle, &nFuselages);
  for (i = 1; i <= nWings; ++i)
  {
    char * wingUID = NULL;
    tiglWingGetUID(tiglHandle, i, &wingUID);
    printf("Wing %d name: %s\n", i, wingUID);
  }

  for (i = 1; i <= nFuselages; ++i)
  {
    char * fuselageUID = NULL;
    tiglFuselageGetUID(tiglHandle, i, &fuselageUID);
    printf("Fuselage %d name: %s\n", i, fuselageUID);
  }

  // query point on the upper wing surface
  if (nWings > 0 && tiglWingGetUpperPoint(tiglHandle, 1, 1, 0.5, 0.5, &x, &y, &z) == TIGL_SUCCESS)
  {
    printf("Point on upper wing surface of wing 1/segment 1: (x,y,z) = (%g, %g, %g)\n", x, y, z);
  }

  // compute intersection of wing with a x-z plane
  if (nWings > 0)
  {
    char * wingUID = NULL;
    char* int_id = NULL;
    int nlines = 0;
    tiglWingGetUID(tiglHandle, 1, &wingUID);
    // do the intersection with a plane at p = (0,0.5,0) and normal vector n = (0,1,0)
    tiglIntersectWithPlane(tiglHandle, wingUID, 0, 0.5, 0, 0, 1, 0, &int_id);

    // get number of intersection wires
    tiglIntersectGetLineCount(tiglHandle, int_id, &nlines);

    // query points on the first line
    if (nlines > 0)
    {
      double zeta = 0.;
      printf("\nIntersection points of a plane with first wing:\n");
      while ( zeta < 1)
      {
        tiglIntersectGetPoint(tiglHandle, int_id, 1, zeta, &x, &y, &z);
        printf("zeta = %g\tp=(%g, %g, %g)\n", zeta, x, y, z);
        zeta += 0.1;
      }
      printf("\n");
    }
  }

  // Export geometry to iges file
  if (strcmp(format, "fmt_iges") == 0)
  {
    if (tiglExportIGES(tiglHandle, exportFileName) != TIGL_SUCCESS)
    {
      printf("Error exporting geometry to IGES!\n");
    }
  }

  // export geometry to step file
  if (strcmp(format, "fmt_step") == 0)
  {
    if (tiglExportSTEP(tiglHandle, exportFileName) != TIGL_SUCCESS)
    {
      printf("Error exporting geometry to STEP!\n");
    }
  }

  // export cpacs
  //tiglSaveCPACSConfiguration("", tiglHandle);

  // close everything
  tiglCloseCPACSConfiguration(tiglHandle);
  tixiCloseDocument(tixiHandle);
  tixiCleanup();

  Py_INCREF(Py_None);
  return Py_None;
}
