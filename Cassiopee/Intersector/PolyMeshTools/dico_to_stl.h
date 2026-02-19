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

#ifndef PYTHON_TO_STL_UTILS
#define PYTHON_TO_STL_UTILS

#include <map>

//// DICO / MAP utils //////////////////////////
 
inline void convert_dico_to_map___int_int_vecint
(
  PyObject *py_zone_to_rid_to_list,
  std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list)
{
  if (PyDict_Check(py_zone_to_rid_to_list))
  {
    //E_Int nzid = PyDict_Size(py_zone_to_rid_to_list);

    PyObject *py_zid/*key*/, *py_rid_to_list /*value : map rid to ptlist*/;
    Py_ssize_t pos = 0;

    while (PyDict_Next(py_zone_to_rid_to_list, &pos, &py_zid, &py_rid_to_list))
    {
      int zid = (int) PyInt_AsLong(py_zid);

      assert (PyDict_Check(py_rid_to_list) == 1); // it s a map

      PyObject *py_rid/*key*/, *py_ptlist /*value : ptlist*/;
      Py_ssize_t pos1 = 0;

      while (PyDict_Next(py_rid_to_list, &pos1, &py_rid, &py_ptlist))
      {
        int rid = (int) PyInt_AsLong(py_rid);

        assert (PyArray_Check(py_ptlist) == 1) ; // it s a numpy
        
        PyArrayObject* pyarr = reinterpret_cast<PyArrayObject*>(py_ptlist);

        npy_intp* dims = PyArray_SHAPE(pyarr);

        E_Int ptl_sz = dims[0];
        
        //long* dataPtr = static_cast<long*>(PyArray_DATA(pyarr));
        E_Int* dataPtr = (E_Int*)PyArray_DATA(pyarr);

        std::vector<E_Int> ptl(ptl_sz);
        for (size_t u=0; u < (size_t)ptl_sz; ++u) ptl[u] = dataPtr[u];

        //std::cout << "max in C is : " << *std::max_element(ALL(ptl)) << std::endl;

        zone_to_rid_to_list[zid][rid]=ptl;

      }
    }
  }
}

inline void convert_dico_to_map__int_pairint
(
  PyObject *py_rid_to_zones,
  std::map<int, std::pair<int,int>>& rid_to_zones)
{
  if (PyDict_Check(py_rid_to_zones))
  {
    // E_Int nzid = PyDict_Size(py_rid_to_zones);

    PyObject *py_rid/*key*/, *py_pair /*value : map zid to ptlist*/;
    Py_ssize_t pos = 0;

    while (PyDict_Next(py_rid_to_zones, &pos, &py_rid, &py_pair))
    {
      int rid = (int) PyInt_AsLong(py_rid);

      assert (PyTuple_Check(py_pair) == 1); // is it a tuple ?

      PyTupleObject* pytup = reinterpret_cast<PyTupleObject*>(py_pair);    

      // -----

      std::pair<int,int> pair_zid;

      PyObject * z1 PyTuple_GET_ITEM(pytup, 0);
      PyObject * z2 PyTuple_GET_ITEM(pytup, 1);

      pair_zid.first  = (double) PyFloat_AsDouble(z1);
      pair_zid.second = (double) PyFloat_AsDouble(z2);

      rid_to_zones[rid] = pair_zid;
    }
  }
}

struct transf_t {
  double t[6];
  bool operator==(const transf_t& r) const
  {
    for (size_t k=0; k < 6; ++k)
      if (t[k] != r.t[k]) return false;
    return true;
  }
  bool operator<(const transf_t& r) const
  {
    if (*this == r) return false;
    for (size_t k=0; k < 6; ++k)
      if (t[k] <r.t[k]) return true;
    return false;
  }
};

inline int convert_dico_to_map___transfo_to_vecint
(
  PyObject *py_transfo_to_list,
  std::map<transf_t, std::vector<int>>& transfo_to_list
)
{
  if (PyDict_Check(py_transfo_to_list) == 0) return 1;
  
  //E_Int nzid = PyDict_Size(transfo_to_list);

  PyObject *py_transfo/*key*/, *py_vecint /*value : vector<int>*/;
  Py_ssize_t pos = 0;

  transf_t t;

  while (PyDict_Next(py_transfo_to_list, &pos, &py_transfo, &py_vecint))
  {
    // key
    assert (PyTuple_Check(py_transfo) == 1) ; // it s a tuple (Xx, Yc, Zc, R)
    PyTupleObject* pytup = reinterpret_cast<PyTupleObject*>(py_transfo);

    for (size_t i=0; i < 6; ++i)
    {
      PyObject * p PyTuple_GET_ITEM(pytup, i);
      t.t[i] = (double) PyFloat_AsDouble(p);
      //std::cout << "transfo " << i << " : " << t.t[i] << std::endl;
    }

    // val
    assert (PyArray_Check(py_vecint) == 1) ; // it s a numpy
    PyArrayObject* pyarr = reinterpret_cast<PyArrayObject*>(py_vecint);

    npy_intp* dims = PyArray_SHAPE(pyarr);

    E_Int sz = dims[0];
    
    //long* dataPtr = static_cast<long*>(PyArray_DATA(pyarr));
    E_Int* dataPtr = (E_Int*)PyArray_DATA(pyarr);

    transfo_to_list[t].resize(sz);
    for (size_t u=0; u < (size_t)sz; ++u) transfo_to_list[t][u] = dataPtr[u];
  }
  return 0;
}

////////////////////////////////////////////////
#endif
