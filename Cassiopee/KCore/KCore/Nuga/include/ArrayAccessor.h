/*    
    Copyright 2013-2024 Onera.

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

#ifndef __NUGA_ACCESSOR_H__
#define __NUGA_ACCESSOR_H__

#include <iostream>
#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"

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

    struct pt_t
    {
      pt_t() {};
      pt_t(value_type* p) {e[0]=p[0]; e[1]=p[1]; e[2]=p[2];}
      pt_t& operator=(pt_t const& r) {e[0]=r.e[0]; e[1]=r.e[1]; e[2]=r.e[2]; return *this;}
      pt_t(pt_t const & r){*this = r;}
      value_type e[3];
      operator value_type*() { return e;}
      //operator typename ArrayType::iterator() { return e;}
      //value_type* operator+(E_Int i) {return e+i;}
    };


  public:
    /// Constuctors

#ifndef NUGALIB
    explicit ArrayAccessor(const array_type& arr, size_type posx, size_type posy, size_type posz = -1, E_Int shift = 0):_arr(&arr), _shift(shift)
    {
      _stride = (posz == -1) ? 2 : 3;
      _posX = new size_type[_stride];
      _posX[0] = posx-NUMFIELD0;  _posX[1] = posy-NUMFIELD0;
      if (_stride == 3)_posX[2] = posz-NUMFIELD0;
    }
#endif
    
    explicit ArrayAccessor(const array_type& arr, E_Int shift = 0):_arr(&arr), _stride(arr.getNfld()), _shift(shift)
    {
      _posX = new size_type[_stride];
      for (E_Int i = 0; i < _stride; ++i)_posX[i]=i;
    }
    
    ArrayAccessor (const ArrayAccessor& aA):_arr(&aA.array()), _stride(aA.stride()), _shift(aA.shift()){
      _posX = new size_type[_stride];
      for (E_Int i = 0; i < _stride; ++i)_posX[i]=aA._posX[i];
    }

    /// Destructor
    ~ArrayAccessor(){delete [] _posX;}


    /// Returns the total number of points.
    inline E_Int size() const {return _arr->getSize();}

    /// Returns the number of used rows
    inline E_Int stride() const {return _stride;}
    /// Returns the number of used rows for the i-th element (NGON specific)
    inline E_Int stride(E_Int i) const {return _stride;}
    
    /// Returns the stride length
    inline E_Int row_stride() const { return _arr->_sizeMax;}
    
    /// Returns the array in read-only
    inline const ArrayType& array() const {return *_arr;}

    ///
    inline E_Int shift() const { return _shift;}
    
    ///
#ifndef NUGALIB
    inline E_Int posX(E_Int i) const { return _posX[i]+NUMFIELD0;}
#endif

    /// Returns the j-th entry's pointer to the i-th field.
    inline value_type* getValPtr(const E_Int& j, const E_Int& i) const
    {
        if (_arr->_compact) // array1
        {
          return _arr->_data + _posX[i]*_arr->_sizeMax + j;
        }
        else // array2/3
        {
          return _arr->_rake[i] + j;
        }
    }
    
    /// Returns the j-th entry to the i-th field.
    inline value_type getVal(const E_Int& j, const E_Int& i) const
    {
        return *getValPtr(j,i) + _shift;
    }

    /// Returns the j-th entry.
    inline void getEntry(const E_Int& j, value_type* entry) const
    {
      for (E_Int k = 0; k < _stride; ++k) entry[k] = getVal(j,k);
    }
    
    /// Returns the j-th entry's pointer to the first field.
    inline const value_type* getEntry(const E_Int& j) const {return getValPtr(j,0);}
    
    /// Returns the j-th entry : MUST BE A FIXED STRIDE (i.e E2N element)
    template <typename ELT>
    inline void getEntry(const E_Int& j, ELT& PHj) const
    {
      E_Int* p = PHj.nodes();
      for (E_Int k = 0; k < _stride; ++k) p[k] = getVal(j,k);
    }

    pt_t col(E_Int j) const { 
      value_type e[3];
      getEntry(j, e);
      return pt_t(e);
    }

    /// Checks whether the index is out of range or not.
    inline E_Bool isOutOfRange(const E_Int& j) const {return (j >= _arr->getSize());}

    /// square of the distance between two nodes.
    inline E_Float dist2(E_Int n, E_Int m) const
    {
      E_Float d2 = 0.;
      for (E_Int k = 0; k < _stride; ++k)
      {
        const E_Float* X = getValPtr(0,k);
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
        const E_Float* X = getValPtr(j,k);
        d2 += (*X - pt[k])*(*X - pt[k]);
      }
      return d2;
    }
    
  private:
    ArrayAccessor();

  protected:
    /// Coordinates array.
    const ArrayType*           _arr;
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
    typedef           T const *           pt_t;//entry type
    
  public:
    /// Constuctors
    explicit ArrayAccessor():_arr(nullptr), _stride(0){}

    explicit ArrayAccessor(const array_type& arr, E_Int dummyshift = 0):_arr(&arr), _stride(arr.rows())
    {}

    /// Destructor
    ~ArrayAccessor(){};

    void set(const array_type& arr){_arr=&arr; _stride = arr.rows();}

    /// Returns the total number of entries (points/elements).
    inline E_Int size() const {return _arr->cols();}
    /// Returns the number of used rows
    inline E_Int stride() const {return _stride;}
    /// Returns the number of used rows for the i-th element (NGON specific)
    inline E_Int stride(E_Int i) const {return _stride;}
    
    /// Returns the array in read-only
    inline const K_FLD::DynArray<T>& array() const {return *_arr;}
        
    ///
    inline E_Int shift() { return 0;}

    /// Returns the (i,j) value. as a matrix pount of view "i" is the row, "j" the column.
    inline T getVal(const E_Int& j, const E_Int& i) const {return *(_arr->begin() + _arr->_rowsMax*j + i);}

    /// Returns the j-th entry (j-th column).
    inline void getEntry(const E_Int& j, T* pE) const
    {
      const typename array_type::value_type* p = (_arr->begin() + _arr->_rowsMax*j);
      for (E_Int k = 0; k < _stride; ++k)
        pE[k] = p[k];
    }

    inline const T* col(E_Int j) const { return _arr->col(j);}
    
    /// Returns the j-th entry's pointer to the first field.
    inline const T* getEntry(const E_Int& j) const { return (_arr->begin() + _arr->_rowsMax*j);}

    /// Checks whether the index is out of range or not.
    inline E_Bool isOutOfRange(const E_Int& j) const {return (j >= _arr->cols());}

    /// square of the distance between two nodes.
    inline E_Float dist2(E_Int n, E_Int m) const
    {
      const E_Float *pn(_arr->begin() + _arr->_rowsMax*n), *pm(_arr->begin() + _arr->_rowsMax*m);
      if (_stride == 3)
        return (pn[0] - pm[0])*(pn[0] - pm[0]) + (pn[1] - pm[1])*(pn[1] - pm[1]) + (pn[2] - pm[2])*(pn[2] - pm[2]);
      else
        return (pn[0] - pm[0])*(pn[0] - pm[0]) + (pn[1] - pm[1])*(pn[1] - pm[1]);
    }

    /// square of the distance between a point and a node.
    inline E_Float dist2(const E_Float* pt, E_Int j) const
    {
      const E_Float *xj(_arr->begin() + _arr->_rowsMax*j);
      if (_stride == 3)
        return (pt[0] - xj[0])*(pt[0] - xj[0]) + (pt[1] - xj[1])*(pt[1] - xj[1]) + (pt[2] - xj[2])*(pt[2] - xj[2]);
      else
        return (pt[0] - xj[0])*(pt[0] - xj[0]) + (pt[1] - xj[1])*(pt[1] - xj[1]);
    }

  protected:
    /// Coordinates array.
    const array_type* _arr;
    /// number of used fields (equal to _posX size)
    size_type         _stride;
  };
  
#ifndef NUGALIB
  inline std::ostream & operator<<(std::ostream& out, const FldArrayF& arr){

    out << "####################################" << std::endl;

    // Print out the matrix.
    for (E_Int i = 0; i < arr.getNfld(); ++i){
      for (E_Int j = 0; j < arr.getSize(); ++j)
        out << ((arr(j,i+1) == IDX_NONE) ? -1 : arr(j,i+1)) << " ";
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
        out << ((arr(j,i+1) == IDX_NONE) ? -1 : arr(j,i+1)) << " ";
      out << std::endl;
    }

    out << std::endl;
    out << "####################################" << std::endl;
    return out;
  }
#endif
}

#endif
