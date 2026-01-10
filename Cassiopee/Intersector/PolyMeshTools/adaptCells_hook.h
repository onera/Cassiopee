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

#ifndef ADAPTCELLS_HOOK
#define ADAPTCELLS_HOOK


#define HMESH_HOOK_ID 77
#define SENSOR_HOOK_ID 78

#define HMESH_PACK_SIZE 6


//=============================================================================
/* get hmesh hook  */
//=============================================================================
inline void* unpackHMesh(PyObject* hook_hmesh, int *&hook_hm_id, int *&subdiv_type, int *&elt_type, int *&zid, std::string *&vString, void **&packet)
{
  //std::cout << "unpackHMesh : begin : " << hook_hmesh << std::endl;

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**)PyCObject_AsVoidPtr(hook_hmesh);
#else
  packet = (void**)PyCapsule_GetPointer(hook_hmesh, NULL);
#endif

  //std::cout << "unpackHMesh : after capsule : " << packet << std::endl;

  if (packet == nullptr)
  {
    PyErr_SetString(PyExc_TypeError,
      "unpackHMesh: PyCapsule_GetPointer failure.");
    return nullptr;
  }

  hook_hm_id = (int*)packet[0];        // type of hook

  //std::cout << "unpackHMesh : after type" << std::endl;

  if (*hook_hm_id != HMESH_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError,
      "unpackHMesh: hook id failure.");
    return nullptr;
  }

  //std::cout << "unpackHMesh : before setting vals" << std::endl;
  
  void* hmesh          = packet[1];                // untyped hmesh ptr
  subdiv_type          = (int*)packet[2];        // subdivision type ISO, ISO_HEX, DIR...  
  elt_type             = (int*)packet[3];        // type of elements in hmesh
  vString              = (std::string*)packet[4];  // for buildArray
  zid                  = (int*) packet[5];

  //std::cout << "unpackHMesh : end" << std::endl;

  return hmesh;
}

//=============================================================================
/* get sensor hook  */
//=============================================================================
inline void* unpackSensor(PyObject* hook_sensor, int *&hook_ss_id, int *&sensor_type, int *&smoothing_type, int *&subdiv_type, int *&elt_type, void **&packet_ss)
{

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet_ss = (void**)PyCObject_AsVoidPtr(hook_sensor);
#else
  packet_ss = (void**)PyCapsule_GetPointer(hook_sensor, NULL);
#endif

  hook_ss_id = (int*)packet_ss[0];

  if (*hook_ss_id != SENSOR_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError,
      "unpackSensor: this function requires a identify sensor hook.");
    return nullptr;
  }

  sensor_type    = (int*)packet_ss[1];
  smoothing_type = (int*)packet_ss[3];
  elt_type       = (int*)packet_ss[4];
  subdiv_type    = (int*)packet_ss[5];
  void* sensor   = packet_ss[2];

  return sensor;
}


#endif
