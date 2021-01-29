/*    
    Copyright 2013-2021 Onera.

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


# include "intersector.h"
# include "stub.h"

//=============================================================================
/* get COM hook  */
//=============================================================================
void* unpackCOM(PyObject* hook, E_Int *&hook_id, E_Int *&subdiv_type, E_Int *&elt_type, void **&packet)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}
//============================================================================
/* Create a hmesh and returns a hook */
//============================================================================
PyObject* K_INTERSECTOR::createHMesh(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//============================================================================
/* Create a communicator */
//============================================================================
PyObject* K_INTERSECTOR::createCom(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* get hmesh hook  */
//=============================================================================
void* unpackHMesh(PyObject* hook_hmesh, E_Int *&hook_hm_id, E_Int *&subdiv_type, E_Int *&elt_type, std::string *&vString, void **&packet)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* 
   Confomize a hmesh
 */
//=============================================================================
PyObject* K_INTERSECTOR::conformizeHMesh(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//============================================================================
/* Deletes a hmesh */
//============================================================================

PyObject* K_INTERSECTOR::deleteHMesh(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//============================================================================
/* Creates a Geom Sensor */
//============================================================================
PyObject* K_INTERSECTOR::createSensor(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* get sensor hook  */
//=============================================================================
void* unpackSensor(PyObject* hook_sensor, E_Int *&hook_ss_id, E_Int *&sensor_type, E_Int *&smoothing_type, E_Int *&subdiv_type, E_Int *&elt_type, void **&packet_ss)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* delete sensor hook  */
//=============================================================================
PyObject* K_INTERSECTOR::deleteSensor(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* assign data to sensor  */
//=============================================================================
PyObject* K_INTERSECTOR::assignData2Sensor(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* Deletes a COM  */
//=============================================================================
PyObject* K_INTERSECTOR::deleteCOM(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}


//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::splitNonStarCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* todo */
//=============================================================================
PyObject* K_INTERSECTOR::adaptBox(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}


//=======================  Intersector/PolyMeshTools/splitCells.cpp ====================
