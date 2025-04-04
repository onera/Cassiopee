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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef __DELAUNAY_MESHER_H__
#define __DELAUNAY_MESHER_H__

#include "Nuga/include/Kernel.h"
#include "Nuga/include/delaunay_preds.h"
#include "Nuga/include/MeshData.h"
#include "Nuga/include/MesherMode.h"
#include "Nuga/include/Metric.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/Refiner.h"
#include "Nuga/include/KdTree.h"
#include "Nuga/include/macros.h"
#include "Nuga/include/Edge.h"
#ifdef E_TIME
#include "Nuga/include/chrono.h"
#endif
#include "Nuga/include/DelaunayMath.h"
#include "Nuga/include/IdTool.h"
#include "Nuga/include/random.h"
#include <list>
#include <deque>
#include <set>

#if defined(DEBUG_METRIC)
#include "iodata.h"
#include <sstream>
#endif
#include <iostream>

#ifdef DEBUG_MESHER
#include "Nuga/include/medit.hxx"
#endif

//#define FORCING_EDGE_ERROR 999
#define NODE_IN_EDGE_ERROR 3
#define COINCIDENT_NODE_ERROR 4

namespace DELAUNAY
{
  struct edge_error_t
  {
    edge_error_t():Ni(IDX_NONE),Nj(IDX_NONE){};
    edge_error_t(const edge_error_t& e)
    {
      Ni = e.Ni;
      Nj = e.Nj;
      nodes = e.nodes;
    }
    E_Int Ni,Nj;
    std::set<E_Int> nodes;
  };

  template <typename T, typename MetricType>
  class Mesher
  {

  public:
    typedef NUGA::size_type                                   size_type;
    typedef NUGA::int_set_type                                int_set_type;
    typedef NUGA::int_vector_type                             int_vector_type;
    typedef NUGA::bool_vector_type                            bool_vector_type;
    typedef NUGA::int_pair_type                               int_pair_type;
    typedef NUGA::int_pair_vector_type                        int_pair_vector_type;
    typedef NUGA::non_oriented_edge_set_type                  non_oriented_edge_set_type;
    typedef K_MESH::Triangle                                  element_type;
    typedef K_MESH::Edge                                      edge_type;
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray>           coord_access_type;
    typedef K_SEARCH::KdTree<>                                tree_type;
    typedef typename DELAUNAY::Kernel<T>                      kernel_type;


  public:
    Mesher();
    Mesher(MetricType& metric);
    ~Mesher(void);

    void set(MetricType& metric);
    void clear();
    void seed_random(long int s);

    E_Int run(MeshData& data);

  private:

    E_Int initialize();

    E_Int triangulate();

    E_Int restoreBoundaries(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                            K_FLD::IntArray& neighbors, int_vector_type& ancestors);

    E_Int setColors(size_type Nbox, MeshData& data);

    E_Int finalize(MeshData& data, E_Int N0);

    E_Int refine();

    E_Int clean_data(MeshData& data, const bool_vector_type& mask);

  private:

    E_Int __compute_bounding_box(const int_vector_type& cloud, E_Float& minX, E_Float& minY, E_Float& maxX, E_Float& maxY);

    E_Int __getPipe(size_type N0, size_type N1, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
      const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, int_set_type& pipe,
      int_pair_vector_type& X_edges);
    
    E_Int __get_xedge_on_shell
    (size_type N0, size_type N1, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
     const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, E_Int& S, E_Int& n, E_Float &tolerance);

    E_Int __forceEdge(size_type N0, size_type N1, /*int_set_type& pipe,*/ int_pair_vector_type& Xedges, 
      const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, K_FLD::IntArray& neighbors,
      int_vector_type& ancestors);    

    E_Int __store_edge_error(E_Int Ni, E_Int Nj, const K_FLD::IntArray& connect, const int_pair_vector_type& Xedges);
  
  protected:

    MeshData*                   _data;

    MetricType*                 _metric;

    NUGA::MeshTool*             _tool;

    coord_access_type*          _posAcc;
    tree_type*                  _tree;

    kernel_type*                _kernel;

    int_set_type                _box_nodes;

    E_Int                       _err;

    E_Int                       _N0;

    NUGA::random                _random;

  public:
    MesherMode        mode;
    std::vector<edge_error_t>  _edge_errors;
    
    
#if defined(DEBUG_MESHER)
    public:
      bool dbg_flag;
#endif
  };

  //
  template <typename T, typename MetricType>
  Mesher<T, MetricType>::Mesher()
    :_data(nullptr), _metric(nullptr), _tool(nullptr), _posAcc(nullptr), _tree(nullptr), _kernel(nullptr)
#ifdef DEBUG_MESHER
    , dbg_flag(false)
#endif
  {}


  //
  template <typename T, typename MetricType>
  Mesher<T, MetricType>::Mesher(MetricType& metric)
    :_data(nullptr), _metric(&metric), _tool(nullptr), _posAcc(nullptr), _tree(nullptr), _kernel(nullptr)
#ifdef DEBUG_MESHER
    , dbg_flag(false)
#endif
  {}

  template <typename T, typename MetricType>
  void Mesher<T, MetricType>::clear()
  {
    _box_nodes.clear();
    _edge_errors.clear();
  }

  template <typename T, typename MetricType>
  void Mesher<T, MetricType>::set(MetricType& metric)
  {
    _metric = &metric;// _metric is not freed as it is not owned
  }

  template <typename T, typename MetricType>
  void Mesher<T, MetricType>::seed_random(long int s)
  {
    _random.srand(s);
  }

  template <typename T, typename MetricType>
  Mesher<T, MetricType>::~Mesher(void)
  {
    // _metric is not freed as it is not owned

    if (_tool)
    {
      delete _tool;
      _tool = nullptr;
    }

    if (_posAcc)
    {
      delete _posAcc;
      _posAcc = nullptr;
    }

    if (_tree)
    {
      delete _tree;
      _tree = nullptr;
    }

    if (_kernel)
    {
      delete _kernel;
      _kernel = nullptr;
    }
  }

  //
  template <typename T, typename MetricType>
  E_Int Mesher<T, MetricType>::run(MeshData& data)
  {
    _err = 0;
    _data = &data;

#ifdef E_TIME
    NUGA::chrono c;
    c.start();
#endif

#ifdef E_TIME
    std::cout << "INIT" << std::endl;
#endif

    // Check data etc...
    _err = initialize();

#ifdef E_TIME
    std::cout << c.elapsed() << std::endl;
    c.start();
    std::cout << "TRIANGULATE" << std::endl;
#endif

    // Triangulate.
    _err = triangulate();

    if (_err)
    {
      if (!mode.silent_errors) std::cout << "Warning: mesher: error triangulating." << std::endl;
#ifdef DEBUG_MESHER
      medith::write("err_tria.mesh", *_data->pos, *_data->connectB, "BAR");
#endif
      return _err;
    }

#ifdef E_TIME
    std::cout << c.elapsed() << std::endl;
    c.start();
    std::cout << "RESTORE BOUNDARIES" << std::endl;
#endif

    //DELAUNAY::iodata::read("out.mdat", data);
    // Restore missing hard edges.
    _err = restoreBoundaries(*data.pos, data.connectM, data.neighbors, data.ancestors);
    if (_err)
    {
      if (!mode.silent_errors) std::cout << "Warning: mesher: error restoring boundaries." << std::endl;
      return _err;
    }
    

#ifdef E_TIME
    std::cout << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_MESHER
    if (dbg_flag)
      medith::write("triangulationC", *_data->pos, _data->connectM, "TRI", &_data->mask);
#endif

#ifdef E_TIME
    std::cout << "SET COLORS" << std::endl;
    c.start();
#endif

    // Find out the subdomains.
    _err = setColors(data.pos->cols()-1, data); //fixme : Nbox
    if (_err)
    {
      if (!mode.silent_errors) std::cout << "Warning: mesher: error setting colors." << std::endl;
      return _err;
    }

#ifdef DEBUG_MESHER
    //if (dbg_flag)
      medith::write("triangulationC_color",*_data->pos, _data->connectM, "TRI", &_data->mask, &data.colors);
#endif
    
    if (mode.mesh_mode == MesherMode::REFINE_MODE)
    {
#ifdef E_TIME
      c.start();
      std::cout << "REFINE" << std::endl;
#endif

      // Refine
      _err = refine();
      if (_err) return _err;
    }

#ifdef E_TIME
    std::cout << c.elapsed() << std::endl;
    c.start();
    std::cout << "FINALIZE MESH" << std::endl;
#endif

    _err = finalize(data, _N0);

#ifdef E_TIME
    std::cout << " " << c.elapsed() << std::endl;
    c.start();
#endif

    if (_err && !mode.silent_errors) std::cout << "Warning: mesher: error finalizing." << std::endl;
    
    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int Mesher<T, MetricType>::initialize()
  {
    // Fast returns.
    if (_err) return _err;

    const K_FLD::IntArray& connectB = *_data->connectB;
    size_type COLS(connectB.cols());
    
    // Check that IDs in connect are referring to column in pos.

    // Check that in 2D, connectB and pos have 2 rows

    //etc...

    // Keep the minimum index in pos that can be changed.
    _N0 = _data->pos->cols();

    // Store the hard edges in a non oriented set.
    E_Int Ni, Nj;
    K_FLD::IntArray::const_iterator pS = connectB.begin();
    std::set<E_Int> hNodes;
    for (size_type i = 0; i < COLS; ++i, pS = pS+2)
    {
      Ni = *pS;
      Nj = *(pS+1);
      _data->hardEdges.insert(K_MESH::NO_Edge(Ni, Nj));
      hNodes.insert(Ni);
      hNodes.insert(Nj);
    }

    // Reset hard nodes to be consistent with hard edges
    hNodes.insert(ALL(_data->hardNodes)); // Append with the input hard nodes.
    _data->hardNodes.clear();
    E_Int idmax=-1;
    for (int_set_type::const_iterator it = hNodes.begin(); it != hNodes.end(); ++it){
      _data->hardNodes.push_back(*it);
      idmax = std::max(idmax, *it);
    }
    
    if (mode.ignore_unforceable_edges) mode.ignore_coincident_nodes=true;
    // 
    if (mode.ignore_coincident_nodes == true)
      K_CONNECT::IdTool::init_inc(_data->hnids, idmax+1);

    // Build the initial mesh.
    E_Float minX, minY, maxX, maxY, L;
    __compute_bounding_box(_data->hardNodes, minX, minY, maxX, maxY);

    L = std::max(maxY-minY, maxX-minX);

    double factor = 0.1;

    minX -= factor*L;
    minY -= factor*L;
    maxX += factor*L;
    maxY += factor*L;

    E_Float c1[2] = {minX,minY};
    E_Float c2[2] = {maxX,minY};
    E_Float c3[2] = {maxX,maxY};
    E_Float c4[2] = {minX,maxY};

    size_type C1,C2,C3,C4;
    assert(_data->pos->rows() == 2);
    _data->pos->pushBack(c1, c1+2); C1 = _data->pos->cols()-1;
    _data->pos->pushBack(c2, c2+2); C2 = _data->pos->cols()-1;
    _data->pos->pushBack(c3, c3+2); C3 = _data->pos->cols()-1;
    _data->pos->pushBack(c4, c4+2); C4 = _data->pos->cols()-1;

    _box_nodes.insert(C1);
    _box_nodes.insert(C2);
    _box_nodes.insert(C3);
    _box_nodes.insert(C4);

    E_Int T1[3] = {C1,C2,C4};
    E_Int T2[3] = {C2,C3,C4};

    _data->connectM.pushBack (T1, T1+3);
    _data->connectM.pushBack (T2, T2+3);

    _data->ancestors.resize(_data->pos->cols(), IDX_NONE);
    _data->ancestors[C1] = _data->ancestors[C2] = _data->ancestors[C4] = 0;
    _data->ancestors[C3] = 1;

    size_type def = IDX_NONE;
    _data->neighbors.resize(3, 2, &def);
    _data->neighbors(0,0) = 1;
    _data->neighbors(1,1) = 0;

    std::vector<E_Int> indices;
    indices.push_back(C1);
    indices.push_back(C2);
    indices.push_back(C3);
    indices.push_back(C4);

    // KdTree initialisation
    if (_posAcc == nullptr)
      _posAcc = new K_FLD::ArrayAccessor<K_FLD::FloatArray>(*_data->pos);
    else
      _posAcc->set(*_data->pos);

    if (_tree == nullptr)
      _tree = new tree_type(*_posAcc, indices); //add the box nodes only.
    else
      _tree->build(&indices);

    if (_tool == nullptr)
      _tool = new NUGA::MeshTool(*_tree);
    else
      _tool->set(*_tree);

    if (_kernel == nullptr)
      _kernel = new kernel_type(*_data, *_tool);
    else
      _kernel->set(*_data, *_tool);

    //fixme
    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int Mesher<T, MetricType>::triangulate()
  {
    // Fast returns.
    if (_err) return _err;

    if (!mode.do_not_shuffle)
      std::shuffle (ALL(_data->hardNodes), _random.gen);

    size_type nb_nodes(_data->hardNodes.size()), Ni;

    std::vector<E_Int> newIds;
    int unconstrained = 0;
    for (size_type i = 0; (i < nb_nodes) && (_err == 0); ++i)
    {
      Ni = _data->hardNodes[i];

      _err = _kernel->insertNode(Ni, (*_metric)[Ni], unconstrained);
      if (_err == 0) _tree->insert(Ni);
      if (_err == 2 && mode.ignore_coincident_nodes)
      {
        _data->unsync_nodes = true;
        _data->hnids[Ni]=_kernel->_Nmatch;
        _err = 0; //to carry on
      }

#ifdef DEBUG_MESHER
      if (dbg_flag)
      {
        std::ostringstream o;
        o << "t_" << i << ".mesh";
        //K_CONVERTER::DynArrayIO::write(o.str().c_str(), _data->pos, _data->connectM, "TRI", &_data->mask);
        medith::write(o.str().c_str(), *_data->pos, _data->connectM, "TRI", &_data->mask);
      }
#endif
    }

    if (_err) return _err;

#ifdef DEBUG_MESHER
    if (dbg_flag)
      medith::write("triangulation0.mesh", *_data->pos, _data->connectM, "TRI", &_data->mask);
#endif

    clean_data(*_data, _data->mask); // Clean by removing invalidated elements and update the data structure. 


#ifdef E_TIME
    std::cout << "cav time        : "           << _kernel->cavity_time << std::endl;
    std::cout << "    init cav time   : "       << _kernel->init_cavity_time << std::endl;
    std::cout << "            base time   : "   << _kernel->_base_time << std::endl;
    std::cout << "            append time : "   << _kernel->_append_time << std::endl;
    std::cout << "    fix cav time    : "       << _kernel->fix_cavity_time << std::endl;
    std::cout << "    sort bound time : "       << _kernel->sorting_bound_time << std::endl;
    std::cout << "remesh time     : "           << _kernel->remesh_time << std::endl;
    std::cout << "inval time      : "           << _kernel->inval_time << std::endl;
    std::cout << std::endl;
#endif

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int Mesher<T, MetricType>::restoreBoundaries
    (const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
    K_FLD::IntArray& neighbors, int_vector_type& ancestors)
  {

    std::vector<K_MESH::NO_Edge> missing_edges;
    int_set_type pipe;
    int_pair_vector_type Xedges;

    size_type cols(connect.cols());
    K_FLD::IntArray::const_iterator pS;

    // Get all the triangulation edges.
    non_oriented_edge_set_type all_edges;
    for (size_type j = 0; j < cols; ++j)
    {
      pS = connect.col(j);
      for (size_type i = 0; i < element_type::NB_NODES; ++i)
        all_edges.insert(K_MESH::NO_Edge(*(pS+i), *(pS + (i+1)%element_type::NB_NODES)));
    }

    // Get the missing hard edges.
    std::set_difference (ALL(_data->hardEdges), ALL(all_edges), std::back_inserter(missing_edges));

    size_type sz = (size_type)missing_edges.size();
    
    E_Int Ni, Nj;
    
    for (size_type i = 0; (i < sz) && !_err; ++i)
    {
      Ni = missing_edges[i].node(0);
      Nj = missing_edges[i].node(1);
      
#ifdef DEBUG_MESHER
      if (dbg_flag)
        medith::write("triangulationi.mesh", pos, connect, "TRI", &(_data->mask));
#endif
      _err = __getPipe(Ni, Nj, pos, connect, neighbors, ancestors, pipe, Xedges);
      
      if (_err)
      {
        if (mode.ignore_unforceable_edges)
        {
          //store this edge and the faulty entities and carry on with other missing edges
          __store_edge_error(Ni,Nj, connect, Xedges);
          _err = 0;
          continue;
        }
        else if (!mode.silent_errors)
          std::cout << "error getting pipe : " << _err << std::endl;
      }

      if (!_err)
        _err = __forceEdge(Ni, Nj, /*pipe,*/ Xedges, pos, connect, neighbors, ancestors);
      
      if (_err)
      {
        if (mode.ignore_unforceable_edges)
        {
          //store this edge and the faulty entities and carry on with other missing edges
          __store_edge_error(Ni, Nj, connect, Xedges);
          _err = 0;
          continue;
        }
        else if (!mode.silent_errors)
          std::cout << "error forcing edge" << std::endl;
      }
    }

    return _err;
  }
  
  ///
  template <typename T, typename MetricType>
  E_Int Mesher<T, MetricType>::__store_edge_error(E_Int Ni, E_Int Nj, const K_FLD::IntArray& connect, const int_pair_vector_type& Xedges)
  {
    edge_error_t e;
    e.Ni = Ni;
    e.Nj = Nj;
    
    for (size_t i=0; i < Xedges.size(); ++i)
    {
      E_Int K = Xedges[i].first;
      E_Int n = Xedges[i].second;
      const E_Int* pK = connect.col(K);
      E_Int N0 = *(pK+(n+1)%3);
      E_Int N1 = *(pK+(n+2)%3);
      
      e.nodes.insert(N0);
      e.nodes.insert(N1);    
    }
    
    _edge_errors.push_back(e);
    
    return 0;
  }

  ///
  template <typename T, typename MetricType>
  E_Int Mesher<T, MetricType>::setColors(size_type Nbox, MeshData& data)
  {
    size_type cols(data.connectM.cols()), Sseed(IDX_NONE), S, Sn, Ni, Nj;
    K_FLD::IntArray::const_iterator pS;

    data.colors.clear();
    data.colors.resize(cols, IDX_NONE);

    // Get an element connected to the box.
    for (size_type i = 0; (i < cols) && (Sseed == IDX_NONE); ++i)
    {
      pS = data.connectM.col(i);
      for (size_type j = 0; (j < element_type::NB_NODES) && (Sseed == IDX_NONE); ++j)
      {
        if (*(pS+j) == Nbox) Sseed = i;
      }
    }

    if (Sseed == IDX_NONE) return 1; // Error

    std::vector<E_Int> cpool;
    size_type color = 0;
    size_type colored = 0;
    while (colored != cols)
    {
      cpool.push_back(Sseed);
      data.colors[Sseed] = color;
      ++colored;

      while (!cpool.empty())
      {
        S = cpool.back();
        cpool.pop_back();
        pS = data.connectM.col(S);

        for (size_type i = 0; i < element_type::NB_NODES; ++i)
        {
          Sn = data.neighbors(i, S);

          if (Sn == IDX_NONE) continue;

          if (data.colors[Sn] != IDX_NONE) continue;

          Ni = *(pS + (i+1) % element_type::NB_NODES);
          Nj = *(pS + (i+2) % element_type::NB_NODES);

          if (_data->hardEdges.find(K_MESH::NO_Edge(Ni, Nj)) != _data->hardEdges.end())
            continue;

          data.colors[Sn] = color;
          ++colored;
          cpool.push_back(Sn);
        }
      }

      ++color;
      if (colored != cols)
      {
        bool found = false;
        size_type i = 0;
        for (; (i < cols) && !found; ++i)
          found = (data.colors[i] == IDX_NONE);
        Sseed = i-1;
      }
    }

    data.mono_connex = (color == 2); // box-elts color + one color

    if ((color > 2) && (mode.remove_holes == true)) // detect eventual interior
    {
      std::vector<K_FLD::IntArray> connects;
      connects.resize(color);
      E_Int c;
      K_FLD::IntArray bound;
      // Split by color.
      for (E_Int Si = 0; Si < data.connectM.cols(); ++Si)
      {
        c = data.colors[Si];
        pS = data.connectM.col(Si);
        if (c >= 1) connects[c].pushBack(pS, pS+3);
      }

      // Store the hard edges in an oriented set.
      E_Int Ni, Nj;
      K_FLD::IntArray::const_iterator pS = data.connectB->begin();
      NUGA::oriented_edge_set_type hard_edges;
      for (size_type i = 0; i < data.connectB->cols(); ++i, pS = pS+2)
      {
        Ni = *pS;
        Nj = *(pS+1);
        hard_edges.insert(K_MESH::Edge(Ni, Nj));
      }

      // Reset the color to zero for interior parts.
      E_Int nbc = (E_Int)connects.size(), c1, Si, nbound;
      bool good_color;
      for (c1 = 1; c1 < nbc; ++c1)
      {
        //DynArrayIO::write("part.mesh", data.pos, connects[c]);
        _tool->getBoundary(connects[c1], bound);
        //medith::write("bound.mesh", *data.pos, bound, "BAR");
        good_color = true;
        nbound = bound.cols();

        for (Si = 0; (Si < nbound) && good_color; ++Si)
          good_color &= (hard_edges.find(K_MESH::Edge(bound.col(Si))) != hard_edges.end());

        if (!good_color)
        {
          for (Si = 0; Si < data.connectM.cols(); ++Si)
          {
            if (data.colors[Si] == c1) data.colors[Si] = 0;
          }
        }
      }
    }

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int Mesher<T, MetricType>::refine()
  {
    
#ifdef DEBUG_METRIC
    {
      std::ostringstream o;
      o << "ellipse_triangulation.mesh";
      _metric->draw_ellipse_field(o.str().c_str(), *_data->pos, _data->connectM, &_data->mask);
      
      std::vector<bool> tmask = _data->mask;
      E_Int cols = _data->connectM.cols();
      tmask.resize(cols);
      for (size_type i = 0; i < cols; ++i) // mask box elements in top of invalidated ones.
      {
        if (tmask[i])
          tmask[i] = (_data->colors[i] != 0);
      }
      
      medith::write("triangulation.mesh", *_data->pos, _data->connectM, "TRI", &tmask);
    }
#endif
    
    K_FLD::FloatArray& pos = *_data->pos;

    _kernel->setConstraint(_data->hardEdges);

    std::vector<size_type> refine_nodes;

    coord_access_type posAcc(pos);
    tree_type filter_tree(posAcc);

    size_type Ni, nb_refine_nodes;

    Refiner<MetricType> saturator(*_metric, mode.growth_ratio, mode.nb_smooth_iter, mode.symmetrize);

#ifdef E_TIME
    NUGA::chrono c;
#endif

    float constrained = 0.;
    E_Int iter = 0;
    bool carry_on = false;
    
    do
    {
      carry_on = false;
#ifdef E_TIME
      c.start();
#endif

#ifdef DEBUG_METRIC
    {
      std::ostringstream o;
      o << "ellipse_refine_iter_" << iter << ".mesh";
      _metric->draw_ellipse_field(o.str().c_str(), *_data->pos, _data->connectM, &_data->mask);
    }
#endif

      // genere des points dans refine_nodes
      saturator.computeRefinePoints(iter, *_data, _box_nodes, _data->hardEdges, refine_nodes, _N0);

#ifdef E_TIME
      std::cout << "__compute_refine_points : " << c.elapsed() << std::endl;
      std::cout << "nb proposed refined points : " << refine_nodes.size() << std::endl;
      c.start();
#endif
#ifdef DEBUG_METRIC
    {
      std::ostringstream o;
      o << "ellipse_refine_points_iter_" << iter << ".mesh";
      _metric->draw_ellipse_field(o.str().c_str(), *_data->pos, refine_nodes);
    }
    //if (iter == 5) saturator._debug = true;
#endif

      // must be here!
      std::shuffle(ALL(refine_nodes), _random.gen);

      // filtre refine_nodes
      saturator.filterRefinePoints(*_data, _box_nodes, refine_nodes, filter_tree);
#ifdef E_TIME
      std::cout << "__filter_refine_points : " << c.elapsed() << std::endl;
      std::cout << "nb points to insert : " << refine_nodes.size() << std::endl;
      c.start();
#endif

      // insere les noeuds resultants
      nb_refine_nodes = refine_nodes.size();
      carry_on = (nb_refine_nodes != 0);

      // must not be here!
      //std::shuffle(ALL(refine_nodes), _random.gen);

      _data->ancestors.resize(_data->pos->cols(), IDX_NONE);

      for (size_type i = 0; (i < nb_refine_nodes) && !_err; ++i)
      {
        Ni = refine_nodes[i];
        _err = _kernel->insertNode(Ni, (*_metric)[Ni], constrained);
        _tree->insert(Ni);
      }

      if (_err) return _err;

#ifdef E_TIME
      std::cout << "insertion : " << c.elapsed() << std::endl;
      c.start();
#endif

      this->clean_data(*_data, _data->mask); // Clean the current mesh by removing invalidated elements.
      
#ifdef DEBUG_METRIC
    {
      std::ostringstream o;
      o << "ellipse_mesh_iter_" << iter << ".mesh";
      _metric->draw_ellipse_field(o.str().c_str(), *_data->pos, _data->connectM);
    }
#endif

#ifdef E_TIME
      std::cout << "compacting : " << c.elapsed() << std::endl;
      std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
#endif
      ++iter;

      //printf("iterating %d\n", iter); fflush(stdout);  
      if (iter > 12) 
      {
        printf("Warning: too much iterations: stopped\n"); fflush(stdout);
        carry_on = false; // hard break
      }
    
    }
    while (carry_on);

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::__compute_bounding_box
    (const int_vector_type& cloud, E_Float& minX, E_Float& minY, E_Float& maxX, E_Float& maxY)
  {
    minX = minY = NUGA::FLOAT_MAX;
    maxX = maxY = - NUGA::FLOAT_MAX;

    K_FLD::FloatArray::const_iterator pK;
    size_type nb_nodes = (size_type)cloud.size();

    for (size_type i = 0; i < nb_nodes; ++i)
    {
      pK = _data->pos->col(cloud[i]);
      minX = std::min(minX, *pK);
      maxX = std::max(maxX, *pK);
      minY = std::min(minY, *(pK+1));
      maxY = std::max(maxY, *(pK+1));
    }

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::__getPipe
    (size_type N0, size_type N1, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
    const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, int_set_type& pipe, int_pair_vector_type& X_edges)
  {
    int_vector_type Ancs;
    size_type ni, S(IDX_NONE), n(IDX_NONE), count(0), nopp, nstart, Ni, Nj;
    K_FLD::IntArray::const_iterator pS;
    E_Bool done(false), intersect;
    E_Float tolerance(_tool->getTolerance())/*fixmetol*/, tol2(tolerance*tolerance)/*fixmetol*//*, u00, u01, u10, u11, bestu=-1.;*//*used in case of poor ancestors to pick the right one*/;
    size_type dim = pos.rows();

    pipe.clear();
    X_edges.clear();
    
    E_Int err = __get_xedge_on_shell(N0, N1, pos, connect, neighbors, ancestors, S, n, tolerance);
    if (err || S==IDX_NONE) return err;
    
    X_edges.push_back(int_pair_type(S,n));
    pipe.insert(S);

    while (!done)
    {
      count  = 0;
      nopp   = element_type::getOppLocalNodeId(S, n, connect, neighbors);
      nstart = (nopp+1)%element_type::NB_NODES;
        
      S = neighbors(n, S);
      pS = connect.col(S);
        
      for (size_type i = 0; (i < edge_type::NB_NODES) && !done; ++i)
      {
        ni  = (nstart+i) % element_type::NB_NODES;
        Ni = *(pS + (ni+1) % element_type::NB_NODES);
        Nj = *(pS + (ni+2) % element_type::NB_NODES);

        done = (Nj == N1); // End of pipe reached.

        intersect  = !done && (element_type::surface(pos.col(N0), pos.col(Ni), pos.col(N1), dim) > tol2);
        intersect &= (element_type::surface(pos.col(N0), pos.col(N1), pos.col(Nj), dim) > tol2);
        // intersect = !done && edge_type::intersect<2>(pos.col(N0), pos.col(N1), pos.col(Ni), pos.col(Nj), tolerance, true/*absolute*/,u00, u01, u10, u11, overlap);
        if (intersect)
        {
          ++count;
          n = ni;
          X_edges.push_back(int_pair_type(S,n));
        }
      }

      if (done && (count != 0))  // Error : we must have Nj = N1 at i = 0 when reaching the end of pipe.
      {err=2; break;}
      if (!done && (count != 1)) // Error : one node is on the edge N0N1.
      {
        err=NODE_IN_EDGE_ERROR;
        break;
      }

      pipe.insert(S);
    }
    
    return err;
  }
  
  inline E_Int sign2D(const E_Float*P0, const E_Float*P1, const E_Float*P)
  {
    E_Float perpdot = (P0[1]-P1[1])*(P[0]-P0[0]) + (P1[0]-P0[0])*(P[1]-P0[1]);
    //return (perpdot < -EPSILON) ? -1 : (perpdot > EPSILON) ? 1 : 0;
    return (perpdot < 0.) ? -1 : (perpdot > 0.) ? 1 : 0;
  }
  
  inline E_Float proj(const E_Float*P0, const E_Float*P1, const E_Float*P)
  {
    return (P1[0]-P0[0])*(P[0]-P0[0]) + (P1[1]-P0[1])*(P[1]-P0[1]);
  }
  
  ///
  template <typename T, typename MetricType>
  E_Int
  Mesher<T, MetricType>::__get_xedge_on_shell
  (size_type N0, size_type N1, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
   const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, E_Int& Sx, E_Int& nx, E_Float &tolerance)
  {
    int_vector_type Ancs;
    size_type Si, ni, Ni, Nj;
    K_FLD::IntArray::const_iterator pS;
    //E_Bool done(false), intersect, overlap;
    E_Float dum1,dum2,dum3;
    E_Bool dum4;
            
    if (ancestors[N0] == IDX_NONE) return 4;
    
    const E_Float* P0=pos.col(N0);
    const E_Float* P1=pos.col(N1);

    pS = connect.col(ancestors[N0]);
    _tool->getAncestors(N0, connect, ancestors, neighbors, std::back_inserter(Ancs));
    
    // Get the ancestor(s) intersecting N0N1. Should be only one, numerical error otherwise.
    size_type AncNb = (size_type)Ancs.size();
    E_Int SIGNi(-2), SIGNip1(-2);
    E_Float ux0=-NUGA::FLOAT_MAX, ux;
    for (size_type i = 0; i < AncNb; ++i)
    {
      Si = Ancs[i];
      pS = connect.col(Si);
      ni = element_type::getLocalNodeId(pS, N0);

      Ni = *(pS + (ni+1)%element_type::NB_NODES);
      Nj = *(pS + (ni+2)%element_type::NB_NODES);

      if ((Ni == N1) || (Nj == N1))
      {
        Sx=IDX_NONE;
        return 0; // N0N1 is already in, so return OK
      }
      
      if (neighbors(ni, Si) == IDX_NONE) continue;
        
      SIGNi = sign2D(P0,P1, pos.col(Ni));
      SIGNip1 = sign2D(P0,P1, pos.col(Nj));
      
      if (SIGNi*SIGNip1 <= 0)
      {
        //ux = proj(P0,P1, pos.col(Nj));
        //bool intersect = 
        edge_type::intersect<2>(P0, P1, pos.col(Ni), pos.col(Nj), tolerance, true/*absolute*/,ux, dum1, dum2, dum3, dum4);
        
        if (SIGNi*SIGNip1 < 0 && ux > 0. && ux < 1.)
        {
          Sx=Si;
          nx=ni;
          return 0;
        }
      
        if (ux0 < ux && ux < 1.1 && ux > -1.1)
        {
          Sx=Si;
          nx=ni;
          ux0=ux;
        } 
      }
    }
    return 0;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::__forceEdge
    (size_type N00, size_type N11, /*int_set_type& pipe,*/ int_pair_vector_type& Xedges,
    const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, K_FLD::IntArray& neighbors,
    int_vector_type& ancestors)
  {
    size_type r, S, b, Sn, bn, Ni, Nj, Nk, Nl, S1, S2, dim(pos.rows()), b1, b2, Xnb;
    size_t nb_try(0);
    K_FLD::IntArray::iterator pS, pSn;
    E_Bool isConvex;
    E_Float tolerance(_tool->getTolerance()), s1, s2, tol2(tolerance*tolerance);
    const size_type& NB_NODES = element_type::NB_NODES;

    while (!Xedges.empty()  && (nb_try < 10*Xedges.size()))
    {
      // Choose randomly an intersecting edge.

      Xnb = (size_type)Xedges.size();
      r   = _random.rand() % Xnb;
      //std::cout << "Mesher : rand : " << r << std::endl;

      int_pair_type& E = Xedges[r];

      S = E.first;
      b = E.second;

      pS = connect.col(S);
      Ni = *(pS + b);
      Nj = *(pS + (b+1) % NB_NODES);
      Nl = *(pS + (b+2) % NB_NODES);

      Sn = neighbors(b, S);

      pSn = connect.col(Sn);
      bn = element_type::getOppLocalNodeId(S, b, connect, neighbors);

      Nk = *(pSn + bn);

      // Convexity test : skip if the pair of element (S, Sn) is not convex (i.e. edge swapping is not doable).
      // Test done by visibility criterion.

      s1 = element_type::surface(pos.col(Ni), pos.col(Nj), pos.col(Nk), dim);
      s2 = element_type::surface(pos.col(Ni), pos.col(Nk), pos.col(Nl), dim);
      isConvex  = (s1 > tol2) && (s2 > tol2);

      if (!isConvex)
      {
        ++nb_try;
        continue;
      }
      nb_try = 0;

      // Neighbors to modify
      S1 = neighbors((b+1) % NB_NODES, S);
      S2 = neighbors((bn+1) % NB_NODES, Sn);
      b1 = element_type::getOppLocalNodeId(S, (b+1) % NB_NODES, connect, neighbors);
      b2 = element_type::getOppLocalNodeId(Sn, (bn+1) % NB_NODES, connect, neighbors);

      // Update elements S and Sn (connect)
      *(pS  + (b+2)  % NB_NODES) = Nk;
      *(pSn + (bn+2) % NB_NODES) = Ni;

      // Update the neighboring (neighbors)
      neighbors((b+1)  % NB_NODES, S)  = Sn;
      neighbors((bn+1) % NB_NODES, Sn) = S;
      if ((S1 != IDX_NONE) && (b1 != IDX_NONE))
        neighbors(b1, S1)              = Sn;
      neighbors(bn, Sn)                = S1;
      if ((S2 != IDX_NONE) && (b2 != IDX_NONE))
        neighbors(b2, S2)              = S;
      neighbors(b, S)                  = S2;

      // Update the ancestors.
      ancestors[Ni] = ancestors[Nj] = ancestors[Nk] = S;
      ancestors[Nl] = Sn;

      // Remove the edge.
      E = Xedges.back(); // compacting : put the last in the current for the next pass
      Xedges.pop_back(); // and remove the last.
      if (--Xnb == 0) return 0;

      // Update (if they exist) the intersecting edge references:
      // (Sn, (bn+1)) -> (S,b)
      // (S, (b+1))   -> (Sn, bn)
      for (size_type i = 0; i < Xnb; ++i)
      {
        int_pair_type& E = Xedges[i];
        if ((E.first == Sn) && (E.second == (bn+1) % NB_NODES))
        {
          E.first = S;
          E.second = b;
        }
        else if ((E.first == S) && (E.second == (b+1) % NB_NODES))
        {
          E.first = Sn;
          E.second = bn;
        }
      }

      // If the three nodes of the new S (or Sn) are in the same side of N0N1, remove S(or Sn).
      // Otherwise insert swapE in Xedges. 
      size_type count = 0;
      E_Float N0N1[2], V[2], v;
      const E_Float* P0 = pos.col(N00);
      NUGA::diff<2> (P0, pos.col(N11), N0N1);

      for (size_type i = 0; i < NB_NODES; ++i)
      {
        NUGA::diff<2> (P0, pos.col(*(pS+i)), V);
        NUGA::crossProduct<2> (N0N1, V, &v);

        if (v > 0.)
          ++count;
        else if (v < 0.)
          --count;
      }

      if (::abs(count) > 1)
      {
        //pipe.erase(S);
        continue;
      }

      count = 0;
      for (size_type i = 0; i < NB_NODES; ++i)
      {
        NUGA::diff<2> (P0, pos.col(*(pSn+i)), V);
        NUGA::crossProduct<2> (N0N1, V, &v);
        if (v > 0.)
          ++count;
        else if (v < 0.)
          --count;
      }

      if (::abs(count) > 1)
      {/*pipe.erase(Sn);*/}
      else
        Xedges.push_back(int_pair_type(S, (b+1) % NB_NODES));// Swapped edge.
    }// While loop

    if (Xedges.empty())
      return 0;
    else
    {
      return 1;//fixme : handle errors in a better way (nb_try...)
    }
  }

  ///
  template <typename T, typename MetricType>
  E_Int
  Mesher<T, MetricType>::clean_data(MeshData& data, const bool_vector_type& mask)
  {
    std::vector<size_type> newIds;

    // Remove invalidated elements
    K_FLD::IntArray::compact (_data->connectM, mask, newIds);
    // Update neighbors
    K_FLD::IntArray::compact (_data->neighbors, newIds);
    K_FLD::IntArray::changeIndices (_data->neighbors, newIds);
    
    // Update colors
    size_type cols(data.connectM.cols()), nb_nodes(data.ancestors.size());
    if (!_data->colors.empty())
    {
      size_type sz = (size_type)newIds.size(), newId;
      for (size_type i = 0; i < sz; ++i)
      {
        newId = newIds[i];
        if (newId != IDX_NONE)
          _data->colors[newId] = _data->colors[i];
      }
      _data->colors.resize(cols);
    }

    // Update mask
    _data->mask.clear();
    _data->mask.resize(_data->connectM.cols(), true);

    // Update ancestors
    for (size_type i = 0; i < nb_nodes; ++i)
    {
      if (_data->ancestors[i] != IDX_NONE)
        _data->ancestors[i] = newIds[_data->ancestors[i]];
    }
    
#ifdef E_DEBUG
    //std::cout << _data->connectM;
    //std::cout << _data->neighbors;
#endif
    
    // synchronize hard entitities regarding the removed nodes
    if (_data->unsync_nodes) // if and only if ignore_coincident_nodes = true
      _data->sync_hards();

    return (cols - data.connectM.cols());
  }

  ///
  template <typename T, typename MetricType>
  E_Int Mesher<T, MetricType>::finalize(MeshData& data, E_Int N0)
  {
    size_type cols(data.connectM.cols());

    data.mask.resize(cols); // ???

    for (size_type i = 0; i < cols; ++i) // mask box elements in top of invalidated ones.
    {
      if (data.mask[i])
        data.mask[i] = (data.colors[i] != 0);
    }

    clean_data(data, data.mask);

    if (cols && data.connectM.cols() == 0) // all the elements are masked
      return 77; // error finalize

    std::vector<E_Int> new_IDs;
    _tool->compact_to_mesh(*data.pos, data.connectM, new_IDs, &N0); // Remove unused nodes above N0.

    // update ancestors.
    for (size_t i = 0; i < data.ancestors.size(); ++i)
    {
      if (new_IDs[i] != IDX_NONE) data.ancestors[new_IDs[i]] = data.ancestors[i];
    }
    data.ancestors.resize(data.pos->cols());

    return 0;
  }
}
#endif
