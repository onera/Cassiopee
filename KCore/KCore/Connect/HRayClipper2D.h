/*
    Copyright 2013-2019 Onera.

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

#ifndef BND_GEOM_HRAYCLIPPER2D_H
#define BND_GEOM_HRAYCLIPPER2D_H

#include "Connect/connect.h"
#include "Def/DefTypes.h"
#include "Fld/DynArray.h"
#include "MeshElement//Polygon.h"
#include <vector>

#define NB_NODES_MAX 8 //ASSUMING that the rays clipping will never create PGs with more than 8 nodes (true for TRI and QUADS).

#ifdef DEBUG_CLIP
#include "IO/io.h"
//#include "chrono.h"
void draw_PG(const char* fname, K_FLD::FloatArray& crd);
#endif

template <typename Connectivity_t>
class HRayClipper2D {
  public:

    // polygon loop
    struct pg_t {
      E_Int _n[NB_NODES_MAX], count;
      pg_t():count(0){}
    };

    ///
    HRayClipper2D(const E_Float* xs,
                  const E_Float* ys,
                  const Connectivity_t& cnt,        //ngon or conhyb
                        bool     conhyb,
                        E_Int    index_start,       // nodes indices minimum in cnt (0 or 1)
                        E_Int    facets_start,      // position of the first polygons in cnt (2 for Cassiope NGON, 0 for conhyb...))
                        E_Int    nb_tot_facets,     // nb of polygons in cnt
                  const E_Int*   to_process,        // list of faces to process (0-based)
                        E_Int    sz_to_process,     // size of the above list
                        E_Float  TOL,
                  const E_Float* rays,
                        E_Int    nrays,
                        bool     force_inclusion=false)
    : _xs(xs), _ys(ys),
      _ngon(cnt.begin()),
      _conhyb(conhyb),
      _index_start(index_start),
      _to_process(to_process),
      _sz_to_process(sz_to_process),
      _tolerance(TOL),
      _nrays(nrays),
      _clip_data(NULL),
      _allocated_clip_sz(0)
      {
        __init(force_inclusion, rays, nrays, facets_start, nb_tot_facets);
        assert (conhyb || (_index_start==1)); // because we use getPosFacets that has this assumption
      }

    ~HRayClipper2D()
    {
      if (_clip_data != NULL)
        delete [] _clip_data;
    }

    /// Main function
    void run(const E_Float* w, E_Int nfield, E_Float* Wb, E_Float* Sb);

  private:
    /// Computes rays vector and node colors
    void __init(bool force_inclusion, const E_Float * rays, E_Int nrays, E_Int facets_start=0, E_Int nb_facets=0);
    /// Resizes appropriately the pg_t array for clipped bits.
    void __init_clip_data(E_Int n);
    /// Refines the input polygon
    void __intersect_hray(const E_Int* pK, E_Int nK, const E_Float* xs, const E_Float* ys, const std::vector<E_Float>& colors,
                         std::vector<E_Float>& hrays, E_Int nrays, E_Float TOL, K_FLD::FloatArray &refine_crd,
                         std::vector<E_Float>& refine_colors, E_Int& nbands, E_Int & C0);
    /// Clips the polygon by the rays set
    void __clip(const K_FLD::FloatArray& refinedPG, const std::vector<E_Float>& colors, E_Int C0, E_Int nbands);
    // Computes a polygon surface
    E_Float __surface(const E_Float* xs, const E_Float* ys, const E_Int* nodes, E_Int nb_nodes);
    /// Sets the colors for each non-already-colored nodes.
    void __set_colors(const E_Int* pK, E_Int nK, const E_Float* ys, std::vector<E_Float>& colors, const std::vector<E_Float>& hrays, E_Int nrays, E_Float TOL);
    /// Inner recursive function to pick by dichotomy the right band, i.e the right color
    E_Float __recurse_col(const std::vector<E_Float>& hrays, E_Float y0, E_Float TOL, E_Int first, E_Int last);
    /// if rays are assumed to surround the full input mesh, this changes the end rays in order to avoid extra clipping
    void force_inclusion_in_rays();

  //
  private:
    const E_Float *_xs, *_ys;
    const E_Int* _ngon;
    const E_Int* _to_process;
    E_Int _sz_to_process;
    E_Float _tolerance;
    E_Int _nrays;
    bool _conhyb;
    E_Int _index_start;
    std::vector<E_Float> _colors, _hrays;
    std::vector<E_Int> _facets;

    pg_t* _clip_data;
    E_Int _allocated_clip_sz;
};

///
template <typename Connectivity_t>
void HRayClipper2D<Connectivity_t>::run(const E_Float* w, E_Int nfield, E_Float* Wb, E_Float* Sb)
{
  //
  E_Int nK, C0, xbands/*intersected bands*/;
  const E_Int* pK;
  K_FLD::FloatArray refine_crd;
  std::vector<E_Float> refine_colors;
  E_Float s0, si;

  E_Int NBANDS=_nrays-1;

  //init
  for (size_t i=0; i < NBANDS; ++i)
    Sb[i]=0.;
  for (size_t i=0; i < NBANDS*nfield; ++i)
    Wb[i]=0.;

  //
  for (size_t i=0; i < _sz_to_process; ++i)
  {
    //
    //std::cout << i << std::endl;
    const E_Int& Ki = _to_process[i];
    //std::cout << _facets[Ki] << std::endl;
    
    if (!_conhyb)
    {
      nK = _ngon[_facets[Ki]];
      pK = &_ngon[_facets[Ki]+1];
    }
    else
    {
      nK = (_ngon[_facets[Ki]+3]==-1 ? 3:4);
      pK = &_ngon[_facets[Ki]];
    }

    __intersect_hray(pK, nK, _xs, _ys, _colors, _hrays, _nrays, _tolerance,
                     refine_crd, refine_colors, xbands, C0);

#ifdef DEBUG_CLIP
    assert (xbands <= NBANDS+2);
    draw_PG("/home/slandier/tmp/hrayclipper/split.mesh", refine_crd);
#endif

    if ((C0 > _nrays) || (C0 < 0 && xbands==1)) // completly outside the rays
    {
#ifdef DEBUG_CLIP
      std::cout << "discarded PG: " << Ki << std::endl;
#endif
       continue;
    }
    
#ifdef DEBUG_CLIP
    if ( (C0 < 0)               || /* some polygon nodes are under the first ray*/
         (C0+xbands > NBANDS)   ) /* some polygon nodes are above the last ray*/
      std::cout << "WARNING : some parts of polygon no " << Ki << " are out of the rays" << std::endl;
#endif
    
     __clip(refine_crd, refine_colors, C0, xbands);

#ifdef DEBUG_CLIP
    /*for (size_t b=0; b<xbands; ++b)
    {
      for (size_t j=0; j< _clip_data[b].count; ++j)
        std::cout << _clip_data[b]._n[j] << " ";
      std::cout << std::endl;
    }*/
#endif

    s0 = __surface(_xs, _ys, pK, nK);

    if (xbands == 1) // Face fully inside band C0
    {
      if (C0 < 0) // under all rays
        continue;
      Sb[C0] += s0;
      for (size_t f=0; f < nfield;++f)
        Wb[C0+f*NBANDS] += s0*w[i+f*_sz_to_process];
    }
    else // Face reparted to xband from C0
    {
      for (E_Int n=0; n < xbands; ++n)
      {
        // Band C0+n
        if (((C0+n) < 0) || ((C0+n) >= NBANDS)) //discard pieces outside the rays.
          continue;
        pg_t& PGi = _clip_data[n];
        si = K_MESH::Polygon::surface<K_FLD::FloatArray,2>(refine_crd, &PGi._n[0], PGi.count, 0);
        //assert(si < s0);
        Sb[C0+n] += si;
        for (size_t f=0; f < nfield;++f)
          Wb[C0+n+f*NBANDS] += (si*w[i+f*_sz_to_process]);
      }
    }
  }
}

///
template <typename Connectivity_t>
void HRayClipper2D<Connectivity_t>::__init(bool force_inclusion, const E_Float * rays, E_Int nrays, E_Int facets_start, E_Int nb_tot_facets){


  // Init ngon facets
  if (!_conhyb)
    K_CONNECT::getPosFacets(_ngon, facets_start, nb_tot_facets, _facets);
  else
  {
    _facets.clear();
    _facets.resize(nb_tot_facets, 1);
    for (E_Int i = 0; i < nb_tot_facets; i++)
    {
      _facets[i] = facets_start+2;
      facets_start += 6;
    }
  }
  
  _hrays.resize(_nrays, 0.);
  for (size_t i=0; i < _nrays; ++i)
    _hrays[i]=rays[i];
  
  if (force_inclusion) // to capture all the pieces by ensuring that botom ansd top rays are surrounding the full mesh
    force_inclusion_in_rays();   

  // set node colors
  E_Int maxId=-1, nb_nodes, id;
  // Get first the max node id to resize _colors
  if (!_conhyb)
  {
    for (E_Int i = 0; i < _sz_to_process; ++i)
    {
      const E_Int& Ki = _to_process[i];
      nb_nodes = _ngon[_facets[Ki]];
      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        id = _ngon[_facets[Ki]+1+j]-_index_start;
        maxId = (maxId < id) ? id : maxId;
      }
    }
  }
  else
  {
    for (E_Int i = 0; i < _sz_to_process; ++i)
    {
      const E_Int& Ki = _to_process[i];
      nb_nodes = (_ngon[_facets[Ki]+3]==-1) ? 3:4;
      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        id = _ngon[_facets[Ki]+j]-_index_start;
        maxId = (maxId < id) ? id : maxId;
      }
    }
  }

  _colors.resize(maxId+1, -1);

  E_Int nK;
  const E_Int* pK;
  //
  if (!_conhyb)
  {
    for (size_t i=0; i < _sz_to_process; ++i)
    {
      const E_Int& Ki = _to_process[i];
      nK = _ngon[_facets[Ki]];
      pK = &_ngon[_facets[Ki]+1];

      // Set vertex colors if required
      __set_colors(pK, nK, _ys, _colors, _hrays, _nrays, _tolerance);
    }
  }
  else
  {
    for (size_t i=0; i < _sz_to_process; ++i)
    {
      const E_Int& Ki = _to_process[i];
      nK = (_ngon[_facets[Ki]+3]==-1 ? 3:4);
      pK = &_ngon[_facets[Ki]];

      // Set vertex colors if required
      __set_colors(pK, nK, _ys, _colors, _hrays, _nrays, _tolerance);
    }
  }
}

template <typename Connectivity_t>
void HRayClipper2D<Connectivity_t>::__init_clip_data(E_Int xbands)
{
  if (_allocated_clip_sz < xbands)
  {
    delete [] _clip_data;
    _clip_data = new pg_t[xbands];
    _allocated_clip_sz = xbands;
    //++nb_resize;
  }
  else
  {
    for (size_t n=0; n < xbands; ++n)
      _clip_data[n].count=0;
  }
}

template <typename Connectivity_t>
E_Float HRayClipper2D<Connectivity_t>::__recurse_col
(const std::vector<E_Float>& hrays, E_Float y0, E_Float TOL, E_Int first, E_Int last)
{
  if (last==(first+1))
  {
    if ((y0 - hrays[first]) <= -TOL) return first-1.; //under all rays
    if ((y0 - hrays[first]) <   TOL) return first;
    if ((hrays[last] - y0)  <= -TOL) return last+1.; //above all rays
    if ((hrays[last] - y0)  <   TOL) return last;
    return first+0.5;
  }

  E_Int middle = (last+first) >> 1;

  if (y0 < hrays[middle])
    return __recurse_col(hrays, y0, TOL, first, middle);
  else
    return __recurse_col(hrays, y0, TOL, middle, last);
}

template <typename Connectivity_t>
void HRayClipper2D<Connectivity_t>::__set_colors(const E_Int* pK, E_Int nK, const E_Float* ys, std::vector<E_Float>& colors, const std::vector<E_Float>& hrays, E_Int nrays, E_Float TOL)
{
  for (size_t i=0; i < nK; ++i)
  {
    E_Int N = *(pK+i)-_index_start;
    E_Float y0 = ys[N];
    if (colors[N] < 0.) colors[N] = __recurse_col(hrays, y0, TOL, 0, nrays);
  }
}
/*
template <typename Connectivity_t>
void HRayClipper2D<Connectivity_t>::__set_colors
(const E_Int* pK, E_Int nK, const E_Float* ys, std::vector<E_Float>& colors, const E_Float* hrays, E_Int nrays, E_Float TOL)
{
  E_Int istart=-1, Ni, i0, i1, h0;
  E_Float y0, c0, y1;

  for (size_t i=0; i < nK; ++i)
  {
    Ni = *(pK+i) -_index_start;
    y0 = ys[Ni];
    if (colors[Ni] < 0.)
      continue;
    istart=i;
    break;
  }

  if (istart==-1)
  {
    y0 = ys[*(pK)-_index_start];
    colors[*(pK)-_index_start] = __recurse_col(hrays, y0, TOL, 0, nrays);
  }

  for (size_t i=0; i < nK; ++i)
  {
    i0 = (istart+i)%nK;
    i1 = (i0+1)%nK;
    c0 = colors[*(pK+i0)-_index_start];
    E_Float &c1 = colors[*(pK+i1)-_index_start];
    if (c1 >-1.)
      continue;
    y0 = ys[*(pK+i0)-_index_start];
    y1 = ys[*(pK+i1)-_index_start];
    if (y1 > y0)
    {
      h0 = (E_Int)c0;
      while ((hrays[h0] < y1) && (h0<nrays))++h0;
      c1=(hrays[h0] < y1) ? h0+0.5 : h0;
    }
    else if (y1 < y0)
    {
      h0 = (E_Int)c0;
      while ((hrays[h0] > y1) && (h0>=0))--h0;
      c1=(hrays[h0] > y1) ? h0+0.5 : h0;
    }
    else
      c1=c0;
  }
}*/

#define FABS(a,b) ( (a-b) < 0. ? b-a : a-b )

template <typename Connectivity_t>
void HRayClipper2D<Connectivity_t>::__intersect_hray
(const E_Int* pK, E_Int nK, const E_Float* xs, const E_Float* ys, const std::vector<E_Float>& colors,
 std::vector<E_Float>& hrays, E_Int nrays, E_Float TOL, K_FLD::FloatArray &refine_crd,
 std::vector<E_Float>& refine_colors, E_Int& nbands, E_Int & C0)
{
  //
  nbands = 1; // number of bands covering the polygon
  C0 = nrays+2;//lowest band color for this polygon
  E_Int C1 = -1; // greatest band color for this polygon ==> nband = C1 - C0

  E_Int Ni, Nip1, Nstart, nr, Cb/*bottom ray color*/, Ct/*top ray color*/, Cbi;
  E_Float dcol, u[2], Pt[2], uyinv;
  bool reversed;

  refine_crd.clear();
  refine_colors.clear();

  // Loop on edges : create intersection points
  for (size_t i=0; i < nK; ++i)
  {
    Nstart = Ni = *(pK+i)-_index_start; //Nstart == first edge index == inital Ni (in case of reverse Ni becomes Nip1)
    Nip1 = *(pK+(i+1)%nK)-_index_start;

    // work with Ni as lowest color
    reversed = (colors[Nip1] < colors[Ni]);
    if (reversed)
      std::swap(Ni, Nip1);

    Cbi = E_Int(colors[Ni]);                                         //bottom ray for Ni
    Cb  = std::min(Cbi, E_Int(colors[Nip1]));                        //botom ray for this edge
    Ct  = std::max(E_Int(colors[Ni]+0.5), E_Int(colors[Nip1]+0.5));  //top ray for this edge

    C0 = std::min(Cb, C0);
    C1 = std::max(Ct, C1);

    if ((Cb == Ct) && (Cbi==Cb))// lying on the same ray (will be recovered by clipping automatically)
      continue;

    Pt[0]=xs[Nstart];
    Pt[1]=ys[Nstart];

    refine_crd.pushBack(Pt,Pt+2);
    refine_colors.push_back(colors[Nstart]);

    if (Ct-Cb <= 1) // belong to the same band
      continue;

    // must be intersected with nr rays

    ++Cb;
    nr = Ct-Cb;
    assert (nr > 0);

    //director
    u[0]=xs[Nip1]-xs[Ni];
    u[1]=ys[Nip1]-ys[Ni];

    uyinv=(1./u[1]);

    if (!reversed)
    {
      for (E_Int r=0; r < nr; ++r)
      {
        //t = uyinv * (hrays[Cb+r] - yNi);
        Pt[1] = hrays[Cb+r];
        Pt[0] = xs[Ni] + uyinv * (Pt[1] - ys[Ni]) * u[0];

        refine_crd.pushBack(Pt, Pt+2);
        refine_colors.push_back(Cb+r);
      }
    }
    else
    {
      for (E_Int r=nr-1; r >=0; --r)
      {
        //t = uyinv * (hrays[Cb+r] - yNi);
        Pt[1] = hrays[Cb+r];
        Pt[0] = xs[Ni] + uyinv * (Pt[1] - ys[Ni]) * u[0];

        refine_crd.pushBack(Pt, Pt+2);
        refine_colors.push_back(Cb+r);
      }
    }
  }
  
  nbands = std::max(C1-C0, nbands);
  
}

#define zMIN(a,b) ((a<b) ? a: b)

///
template <typename Connectivity_t>
void HRayClipper2D<Connectivity_t>::__clip
(const K_FLD::FloatArray& refinedPG, const std::vector<E_Float>& colors, E_Int C0, E_Int nbands)
{
  E_Int sz(refinedPG.cols());

  __init_clip_data(nbands);

  //
  for (size_t p=0; p< sz; ++p)
  {
    const E_Int& Ni = p;
    const E_Int& Nip1 = (p+1)%sz;
    E_Int col = zMIN(colors[Ni],colors[Nip1]) - C0;

#ifdef DEBUG_CLIP
    assert ((col > -1) && (col < nbands));
#endif

    E_Int* n = _clip_data[col]._n;
    E_Int& c = _clip_data[col].count;
    n[c]=Ni; n[c+1]=Nip1;
    c +=2;
  }
}

///
template <typename Connectivity_t>
E_Float HRayClipper2D<Connectivity_t>::__surface(const E_Float* xs, const E_Float* ys, const E_Int* nodes, E_Int nb_nodes)
{
  E_Float s(0.);
  E_Int N2, N3;
  // N1 = origin = (xs[0], ys[0])
  for (size_t n=0; n < nb_nodes; ++n)
  {
    N2 = *(nodes+n)-_index_start;
    N3 = *(nodes+(n+1)%nb_nodes)-_index_start;
    s += 0.5 * ((xs[0]-xs[N3])*(ys[N2]-ys[N3]) - (ys[0]-ys[N3])*(xs[N2]-xs[N3]));
  }

  return s;
}

///
template <typename Connectivity_t>
void HRayClipper2D<Connectivity_t>::force_inclusion_in_rays()
{
  E_Int id;
  E_Float miny(K_CONST::E_MAX_FLOAT), maxy(-K_CONST::E_MAX_FLOAT),y;
  
  if (!_conhyb)
  {
    for (E_Int i = 0; i < _sz_to_process; ++i)
    {
      const E_Int& Ki = _to_process[i];
      const E_Int& nb_nodes = _ngon[_facets[Ki]];
      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        id = _ngon[_facets[Ki]+1+j]-_index_start;
        y=_ys[id];
        miny = (y<miny) ? y : miny;
        maxy = (maxy < y) ? y : maxy;
      }
    }
  }
  else
  {
    for (E_Int i = 0; i < _sz_to_process; ++i)
    {
      const E_Int& Ki = _to_process[i];     
      const E_Int& nb_nodes = (_ngon[_facets[Ki]+3]==-1) ? 3:4;
      
      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        id = _ngon[_facets[Ki]+j]-_index_start;
        y=_ys[id];
        miny = (y<miny) ? y : miny;
        maxy = (maxy < y) ? y : maxy;
        
      }
    }
  }
  
  _hrays[0] = (miny < _hrays[0]) ? miny : _hrays[0];
  _hrays[_nrays-1] = (_hrays[_nrays-1] < maxy) ? maxy : _hrays[_nrays-1];
  
}

#ifdef DEBUG_CLIP
void draw_PG(const char* fname, K_FLD::FloatArray& crd)
{
  K_FLD::IntArray edges;//(2, PG.size());
  E_Int E[2];
  size_t sz = crd.cols();
  for (size_t j=0; j < sz; ++j)
  {
    E[0]=j;
    E[1]=(j+1)%sz;
    edges.pushBack(E, E+2);
  }
  MIO::write(fname, crd, edges, "BAR");
}
#endif

#endif  /* BND_GEOM_HRAYCLIPPER2D_H */
