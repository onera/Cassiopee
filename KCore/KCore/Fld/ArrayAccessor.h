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
#ifndef __KCORE_ARRAY_ACCESSOR_H__
#define __KCORE_ARRAY_ACCESSOR_H__

#include <iostream>
#include "Def/DefTypes.h"
#include "Fld/DynArray.h"

namespace K_FLD
{
  /// Base used by FldArrays
  template <typename  ArrayType>
  class ArrayAccessor
  {
  public: /** Typedefs */
    typedef           E_Int                 size_type;
    typedef           ArrayType             array_type;
    typedef  typename ArrayType::value_type value_type;

  public:
    /// Constuctors
    explicit ArrayAccessor(const array_type& arr, size_type posx, size_type posy, size_type posz = -1, E_Int shift = 0):_arr(arr), _shift(shift)
    {
      _stride = (posz == -1) ? 2 : 3;
      _posX = new size_type[_stride];
      _posX[0] = posx-NUMFIELD0;  _posX[1] = posy-NUMFIELD0;
      if (_stride == 3)_posX[2] = posz-NUMFIELD0;
    }
    
    explicit ArrayAccessor(const array_type& arr, E_Int shift = 0):_arr(arr), _stride(arr.getNfld()), _shift(shift)
    {
      _posX = new size_type[_stride];
      for (E_Int i = 0; i < _stride; ++i)_posX[i]=i;
    }
    
    ArrayAccessor (const ArrayAccessor& aA):_arr(aA.array()), _stride(aA.stride()), _shift(aA.shift()){
      _posX = new size_type[_stride];
      for (E_Int i = 0; i < _stride; ++i)_posX[i]=aA._posX[i];
    }

    /// Destructor
    ~ArrayAccessor(){delete [] _posX;}


    /// Returns the total number of points.
    inline E_Int size() const {return _arr.getSize();}

    /// Returns the number of used rows
    inline E_Int stride() const {return _stride;}
    /// Returns the number of used rows for the i-th element (NGON specific)
    inline E_Int stride(E_Int i) const {return _stride;}
    
    /// Returns the stride length
    inline E_Int row_stride() const { return _arr._sizeMax;}
    
    /// Returns the array in read-only
    inline const ArrayType& array() const {return _arr;}

    ///
    inline E_Int shift() const { return _shift;}
    
    ///
    inline E_Int posX(E_Int i) const { return _posX[i]+NUMFIELD0;}

    /// Returns the i-th field of the j-th entry.
    inline value_type getVal(const E_Int& j, const E_Int& i) const {return *((_arr._data+((_posX[i])*_arr._sizeMax)) + j) + _shift;}

    /// Returns the j-th entry.
    inline void getEntry(const E_Int& j, value_type* entry) const
    {for (E_Int k = 0; k < _stride; ++k)entry[k] =*((_arr._data+((_posX[k])*_arr._sizeMax)) + j) + _shift;}
    /// Returns the j-th entry's pointer to the first field.
    inline const value_type* getEntry(const E_Int& j) const { return _arr._data+j+(_posX[0])*_arr._sizeMax;}
    
    /// Returns the j-th entry : MUST BE A FIXED STRIDE (i.e E2N element)
    template <typename ELT>
    inline void getEntry(const E_Int& j, ELT& PHj) const
    {
      E_Int* p = PHj.nodes();
      for (E_Int k = 0; k < _stride; ++k)p[k] =*((_arr._data+((_posX[k])*_arr._sizeMax)) + j) + _shift;
    }

    /// Checks whether the index is out of range or not.
    inline bool isOutOfRange(const E_Int& j) const {return (j >= _arr.getSize());}

    /// square of the distance between two nodes.
    inline E_Float dist2(E_Int n, E_Int m) const
    {
      E_Float d2 = 0.;
      for (E_Int k = 0; k < _stride; ++k)
      {
        const E_Float* X = (_arr._data+((_posX[k])*_arr._sizeMax));
        d2 += (*(X+n) - *(X+m))*(*(X+n) - *(X+m));
      }
      return d2;
    }

    /// square of the distance between a point and a node.
    inline E_Float dist2(const E_Float* pt, E_Int j) const
    {
      E_Float d2 = 0.;
      for (E_Int k = 0; k < _stride; ++k)
      {
        const E_Float* X = (_arr._data+((_posX[k])*_arr._sizeMax))+j;
        d2 += (*X - pt[k])*(*X - pt[k]);
      }
      return d2;
    }
    
  private:
    ArrayAccessor();

  protected:
    /// Coordinates array.
    const ArrayType&           _arr;
    /// number of used fields (equal to _posX size)
    size_type                  _stride;
    /// fields positions
    size_type*                  _posX;
    ///
    E_Int                       _shift; //node inidices in FldArrays starts at 1
  };

  // Specialisation for FloatArray and IntArray
  template <typename T>
  class ArrayAccessor<K_FLD::DynArray<T> >
  {

  public: /** Typedefs */
    typedef typename  K_FLD::DynArray<T>  array_type;
    typedef           E_Int               size_type;

  public:
    /// Constuctor
    explicit ArrayAccessor(const array_type& arr, E_Int dummyshift = 0):_arr(arr), _stride(arr.rows())
    {}

    /// Destructor
    ~ArrayAccessor(){};

    /// Returns the total number of entries (points/elements).
    inline E_Int size() const {return _arr.cols();}
    /// Returns the number of used rows
    inline E_Int stride() const {return _stride;}
    /// Returns the number of used rows for the i-th element (NGON specific)
    inline E_Int stride(E_Int i) const {return _stride;}
    
    /// Returns the array in read-only
    inline const K_FLD::DynArray<T>& array() const {return _arr;}
        
    ///
    inline E_Int shift() { return 0;}

    /// Returns the (i,j) value. as a matrix pount of view "i" is the row, "j" the column.
    inline T getVal(const E_Int& j, const E_Int& i) const {return *(_arr.begin() + _arr._rowsMax*j + i);}

    /// Returns the j-th entry (j-th column).
    inline void getEntry(const E_Int& j, T* pE) const
    {
      const typename array_type::value_type* p = (_arr.begin() + _arr._rowsMax*j);
      for (E_Int k = 0; k < _stride; ++k)
        pE[k] = p[k];
    }
    
    /// Returns the j-th entry's pointer to the first field.
    inline const T* getEntry(const E_Int& j) const { return (_arr.begin() + _arr._rowsMax*j);}

    /// Checks whether the index is out of range or not.
    inline bool isOutOfRange(const E_Int& j) const {return (j >= _arr.cols());}

    /// square of the distance between two nodes.
    inline E_Float dist2(E_Int n, E_Int m) const
    {
      const E_Float *pn(_arr.begin() + _arr._rowsMax*n), *pm(_arr.begin() + _arr._rowsMax*m);
      if (_stride == 3)
        return (pn[0] - pm[0])*(pn[0] - pm[0]) + (pn[1] - pm[1])*(pn[1] - pm[1]) + (pn[2] - pm[2])*(pn[2] - pm[2]);
      else
        return (pn[0] - pm[0])*(pn[0] - pm[0]) + (pn[1] - pm[1])*(pn[1] - pm[1]);
    }

    /// square of the distance between a point and a node.
    inline E_Float dist2(const E_Float* pt, E_Int j) const
    {
      const E_Float *xj(_arr.begin() + _arr._rowsMax*j);
      if (_stride == 3)
        return (pt[0] - xj[0])*(pt[0] - xj[0]) + (pt[1] - xj[1])*(pt[1] - xj[1]) + (pt[2] - xj[2])*(pt[2] - xj[2]);
      else
        return (pt[0] - xj[0])*(pt[0] - xj[0]) + (pt[1] - xj[1])*(pt[1] - xj[1]);
    }
    

  private:
    ArrayAccessor();

  protected:
    /// Coordinates array.
    const array_type& _arr;
    /// number of used fields (equal to _posX size)
    size_type         _stride;
  };
  
  inline std::ostream & operator<<(std::ostream& out, const FldArrayF& arr){

    out << "####################################" << std::endl;

    // Print out the matrix.
    for (E_Int i = 0; i < arr.getNfld(); ++i){
      for (E_Int j = 0; j < arr.getSize(); ++j)
        out << ((arr(j,i+1) == E_IDX_NONE) ? -1 : arr(j,i+1)) << " ";
      out << std::endl;
    }

    out << std::endl;
    out << "####################################" << std::endl;
    return out;
  }
  
  inline std::ostream & operator<<(std::ostream& out, const FldArrayI& arr){

    out << "####################################" << std::endl;

    // Print out the matrix.
    for (E_Int i = 0; i < arr.getNfld(); ++i){
      for (E_Int j = 0; j < arr.getSize(); ++j)
        out << ((arr(j,i+1) == E_IDX_NONE) ? -1 : arr(j,i+1)) << " ";
      out << std::endl;
    }

    out << std::endl;
    out << "####################################" << std::endl;
    return out;
  }

}

#endif
