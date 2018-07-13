/*    
    Copyright 2013-2018 Onera.

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

#ifndef __CADVIAOCC_H__
#define __CADVIAOCC_H__

#include "TopoDS_Face.hxx"
#include "TopExp.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "Geom_Surface.hxx"
#include "Def/DefTypes.h"
#include <vector>
#include <map>
#include "Fld/DynArray.h"
#include "OCCSurface.h"

namespace K_OCC
{

///
class CADviaOCC
{
public:
	CADviaOCC();
	~CADviaOCC();
    
    ///
    E_Int import_cad(const char* fname, const char*format, E_Float h=0., E_Float chordal_err=0.);
    ///
    E_Int compute_h_sizing(K_FLD::FloatArray& coords, std::vector<E_Int>& Ns);
    ///
    E_Int update_with_chordal_sizing(std::vector<E_Int>& Ns);
    ///
    E_Int mesh_edges(K_FLD::FloatArray& coords, std::vector<K_FLD::IntArray>& connectEs);
    ///
    E_Int build_loops (K_FLD::FloatArray& coords, const std::vector<K_FLD::IntArray>& connectEs, std::vector<K_FLD::IntArray>& connectBs);
    ///
    E_Int mesh_faces(K_FLD::FloatArray& coords, const std::vector<K_FLD::IntArray>& connectEs, std::vector<K_FLD::FloatArray>& crds, std::vector<K_FLD::IntArray>& connectMs);
        
private:
  
    E_Int __build_graph(const TopoDS_Shape& occ_shape, std::vector<OCCSurface*>& vFG);
    void __traverse_face_edges(const TopoDS_Face& F, TopExp_Explorer& edge_expl, std::vector<E_Int>& edges);
    E_Int __split_surface_of_revolution(const TopoDS_Face&, std::vector<OCCSurface*>& vFG, E_Int nb_solid=-1);
    
    E_Int __h_sizing(const TopoDS_Edge& E, E_Float& L);
    E_Int __chord_sizing(const TopoDS_Edge& E, E_Float chordal_err, E_Int& nb_points);
    
    E_Int __clean(const K_FLD::ArrayAccessor<K_FLD::FloatArray>& crdA, K_FLD::IntArray& connectB, E_Float& tol2);
    
    E_Int __mesh_edge(const TopoDS_Edge& E, E_Int& N, K_FLD::FloatArray& coords, K_FLD::IntArray& connectEs);
    E_Int __remove_degenerated(K_FLD::IntArray& connectE);
    
    void __computeOrient(const K_FLD::FloatArray crd2D, const K_FLD::IntArray& cnt, E_Int&o);
    
    E_Float _chordal_err, _h, _merge_tol, _Lmin, _Lmax, _Lmean;
    bool _hrelative;
    
    TopoDS_Shape _occ_shape;
    TopTools_IndexedMapOfShape _surfs, _edges;
    
	std::vector<OCCSurface*> _faces;
    
    // to avoid reallocs
    std::vector<E_Int> _end_nodes;
    std::map<E_Int, E_Int> _enodes;
    
    K_FLD::FloatArray _coord2D;
};

}

#endif
