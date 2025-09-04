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

#ifndef __KCORE_ARRAY_WRITER_H__
#define __KCORE_ARRAY_WRITER_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"
#include "Nuga/include/ArrayAccessor.h"

namespace K_FLD
{
  /// Base used by FldArrays
  template <typename  ArrayType>
  class ArrayWriter
  {
  public: /** Typedefs */
    typedef           E_Int                 size_type;
    typedef           ArrayType             array_type;
    typedef  typename ArrayType::value_type value_type;

  public:
    /// Constuctor
    explicit ArrayWriter(array_type& arr, size_type posx, size_type posy, size_type posz = -1, E_Int shift = 0):_arr(arr)
    {
      _stride = (posz == -1) ? 2 : 3;
      _posX = new size_type[_stride];
      _posX[0] = posx-NUMFIELD0;  _posX[1] = posy-NUMFIELD0;
      if (_stride == 3)_posX[2] = posz-NUMFIELD0;
    }
    
    explicit ArrayWriter(array_type& arr, E_Int shift = 0):_arr(&arr), _stride(arr.getNfld())
    {
      _posX = new size_type[_stride];
      for (size_t i = 0; i < _stride; ++i)_posX[i]=i;
    }
    
    /// Destructor
    ~ArrayWriter(){delete [] _posX;}


    /// Returns the total number of points.
    inline E_Int size() const {return _arr.getSize();}

    /// Returns the number of used rows
    inline E_Int stride() const {return _stride;}
    
    /// Returns the array in read-only
    inline ArrayType& array() const {return *_arr;}

    /// Returns the i-th field of the j-th entry.
    inline value_type& getVal(E_Int j, E_Int i) {return *(_arr.begin(_posX[i]+NUMFIELD0) + j);}

    /// Returns the j-th entry.
    inline void getEntry(const E_Int& j, value_type* entry) const
    {for (E_Int k = 0; k < _stride; ++k)entry[k] =*(_arr.begin(_posX[k]+NUMFIELD0) + j);}
    
  private:
    ArrayWriter();

  protected:
    /// Coordinates array.
    ArrayType&                  _arr;
    /// number of used fields (equal to _posX size)
    size_type                   _stride;
    /// fields positions
    size_type*                  _posX;
  };

  // Specialisation for FloatArray and IntArray
  template <typename T>
  class ArrayWriter<K_FLD::DynArray<T> >
  {

  public: /** Typedefs */
    typedef typename  K_FLD::DynArray<T>  array_type;
    typedef           E_Int               size_type;

  public:
    /// Constuctor
    explicit ArrayWriter(array_type& arr, E_Int dummyshift = 0):_arr(arr), _stride(arr.rows()){}
    
    ArrayWriter (const ArrayAccessor<K_FLD::DynArray<T> >& accessor, array_type& arr):_arr(arr), _stride(accessor.stride()){}

    /// Destructor
    ~ArrayWriter(){};

    /// Returns the total number of entries (points/elements).
    inline E_Int size() const {return _arr.cols();}
    /// Returns the number of used rows
    inline E_Int stride() const {return _stride;}
    
    /// Returns the array in read-only
    inline K_FLD::DynArray<T>& array() {return *_arr;}
    
    inline void reserve(E_Int stride, E_Int n) { _arr.reserve(stride, n);}
    
    inline void clear() { _arr.clear();}
    
    inline void push_back(T* ptr, E_Int n) {_arr.pushBack(ptr, ptr+n);}
    ///
    //inline E_Int shift() const { return 0;}

    /// Returns the (i,j) value. as a matrix pount of view "i" is the row, "j" the column.
    inline T& getVal(E_Int j, E_Int i) const {return *(_arr.begin() + _stride*j + i);}

    /// Returns the j-th entry (j-th column).
    inline void getEntry(const E_Int& j, T* pE) const
    {
      const typename array_type::value_type* p = (_arr.begin() + _stride*j);
      for (E_Int k = 0; k < _stride; ++k)
        pE[k] = p[k];
    }

  private:
    ArrayWriter();

  protected:
    /// Coordinates array.
    array_type&  _arr;
    /// number of used fields (equal to _posX size)
    size_type    _stride;
  };
}

#endif
