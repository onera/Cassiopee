/*
 
 
 
              NUGA 
 
 
 
 */

#include "Fld/DynArray.h"
#include "Fld/ngon_t.hxx"

using ngon_type = ngon_t<K_FLD::IntArray>;

// DIMS : 1(1D), 15("1.5D" : 3D-lineic), 2(2D), 25("2.5D" : 3D-surfacic), 3(3D)

// : default impl : fixed stride (Basic element in any DIM)
template <short DIM, bool fixed_stride>
struct connect_trait
{
  using cnt_t = K_FLD::IntArray;
  static int ncells(const cnt_t& c){return c.cols();}
};

// SURF == 2.5D
template <>
struct connect_trait<25, false>
{
  using cnt_t = ngon_unit;
  static int ncells(const cnt_t& c){return c.size();}
};
// VOL == 3D
template <>
struct connect_trait<3, false>
{
  using cnt_t = ngon_type;
  static int ncells(const cnt_t& c){return c.PHs.size();}
};

template <short DIM, bool FIXSTRIDE>
struct mesh_t
{
  using cnt_t = typename connect_trait<DIM, FIXSTRIDE>::cnt_t;
  K_FLD::FloatArray  crd;
  cnt_t              cnt;
  
  int ncells() const {return connect_trait<DIM, FIXSTRIDE>::ncells(cnt);}
  
  void append(const mesh_t<DIM, FIXSTRIDE>& m)
  {
    cnt_t mcnt = m.cnt;//fixme copy just because of the non-constness of the shift : shift should alos ahev an option to start from an inex to make it at the end
    mcnt.shift(crd.cols());
    crd.pushBack(m.crd);
    cnt.append(mcnt);
  }
  
};
