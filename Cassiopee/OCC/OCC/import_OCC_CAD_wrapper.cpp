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
//Author : Sâm Landier (sam.landier@onera.fr)

#include "import_OCC_CAD_wrapper.h"
#include "CADviaOCC.h"

E_Int K_OCC::import_OCC_CAD_wrapper::import_cad
(
  const char* fname, const char* format,
  std::vector<K_FLD::FloatArray> & crds,
  std::vector<K_FLD::IntArray>& connectMs,
  E_Float h, E_Float chordal_err, E_Float gr,
  bool aniso, bool do_join)
{
#ifdef DEBUG_CAD_READER
  std::cout << "import_OCC_CAD_wrapper::import_cad..." << std::endl;
  std::cout << "import_cad..." << std::endl;
#endif
  
  // CAD --> OCC Shape with associated homemade graph to link flat storage ids between faces and edges.
  CADviaOCC reader;
  E_Int err = reader.import_cad(fname, format, h, chordal_err, gr);
  if (err) return err;
  
#ifdef DEBUG_CAD_READER
  std::cout << "import_cad done." << std::endl;
  std::cout << "mesh_edges..." << std::endl;
#endif
    
  // Mesh the edges.
  std::vector<K_FLD::IntArray> connectEs;
  K_FLD::FloatArray coords;
  err = reader.mesh_edges(coords, connectEs);
  if (err) return err;
  
#ifdef DEBUG_CAD_READER
  std::cout << "mesh_edges done." << std::endl;
  std::cout << "build_loops..." << std::endl;
#endif
  
  // Prepare loops.
  std::vector<K_FLD::IntArray> connectBs;
  err = reader.build_loops(coords, connectEs, connectBs);
  if (err) return err;
  
#ifdef DEBUG_CAD_READER
  std::cout << "build_loops done." << std::endl;
  std::cout << "mesh_faces..." << std::endl;
#endif
    
  // Mesh the surfaces.
  connectMs.clear();
  err = reader.mesh_faces(coords, connectBs, crds, connectMs, aniso, do_join);
  
#ifdef DEBUG_CAD_READER 
  std::cout << "import_OCC_CAD_wrapper::import_cad done." << std::endl;
#endif

  return err;
}

K_OCC::import_OCC_CAD_wrapper::import_OCC_CAD_wrapper() {
}

K_OCC::import_OCC_CAD_wrapper::import_OCC_CAD_wrapper(const import_OCC_CAD_wrapper& orig) {
}

K_OCC::import_OCC_CAD_wrapper::~import_OCC_CAD_wrapper() {
}

