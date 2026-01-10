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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef _KCORE_DYNARRAY_H_
#define _KCORE_DYNARRAY_H_
#include <iterator>
#include "Nuga/include/defs.h"
#include "Fld/FldArray.h"
#include "Nuga/include/allocator.hxx"

#include <assert.h>
#include <stdio.h>
#include <ostream>
#include <set>
#include <vector>
#include <algorithm>

namespace K_FLD 
{

// ============================================================================
// @Name DynArray<T>
// @Memo Container of type T (E_Float, E_Int, E_Bool...)
/* @Text

Design
> Container of type T
> Basic encapsulation of a (single dimension) dynamic C array
> (attribute _data)

> Caution : indices for elements and number of fields go from 0 to size-1.
  raw data is stored column-wise.
>
*/
// ============================================================================

  template <typename T> class ArrayAccessor; //forward dec for frienship (gcc)
  
  template <typename T>
  class DynArray
  {
    /// To output the DynArray data.
    template<typename U> friend std::ostream& operator<<(std::ostream&, const DynArray<U>&);

  public:
    template<typename U> friend class ArrayAccessor;
    
  
  public: /** Typedefs */

    typedef       DynArray       self_type;
    typedef       T              value_type;
    typedef       value_type*    iterator;
    typedef const value_type*    const_iterator;
    typedef       E_Int          size_type;
    typedef       const_iterator pt_t;

  public: /** Constructors and Destructor */

    /// Default Constructor.
    explicit DynArray();
    ///Constructor with a specified allocation option
    DynArray(bool calloc); // false => CPP, true => C
    /// Constructor with number of rows (number of variables) and columns.
    DynArray(size_type rows, size_type cols, bool calloc = false);
    /// Constructor with number of rows, columns and default value.
    DynArray(size_type rows, size_type cols, const value_type& val, bool calloc = false);
    /// Constructor with input data : WARNING becomes owner
    DynArray(T*& compactdata, size_type rows, size_type cols, bool calloc = false);

    /// Copy constructor (cast if U is not T).
    template <typename U> explicit DynArray(const DynArray<U>& rhs);
    ///
    DynArray(const DynArray<T>& rhs);//fixme : why required ? i.e why the above doesn't work (not called when returning an object created inside a function scope)
    DynArray(DynArray<T>&& rhs); //move version

    /// Destructor.
    ~DynArray(){__destroy();}

    void set_alloc(const self_type& arr) {
      if (_data != nullptr) return; //has allocated data so cannot change style
      _calloc = arr._calloc;
    }

    /// Constructor by type conversion (for FldArrays).
    inline explicit DynArray(const FldArray<T>& i, value_type shift = value_type(0.));
    inline explicit DynArray(const FldArray<T>& i, E_Int posx, E_Int posy, E_Int posz=-1);

    /// Converts the DynArray to a FldArray.
    inline void convert(FldArray<T>& i, value_type shift = value_type(0.)) const ;
    inline void convert(DynArray<T>& out, T shift) const {out = *this;/*fixme : dum to convert to itself*/}

    // Extracts a given field
    inline void extract_field(size_type i, std::vector<T>& f);
    // Sets a given field on an existing row
    inline void set_field(size_type i, const std::vector<T>& f);

  public: /** Attributes */

    /// Gets the number of columns.
    inline size_type cols() const {return _cols;}
    inline size_type getSize() const {return _cols;} //to match FldArray interface
    /// Gets number of rows.
    inline size_type rows() const {return _rows;}
    /// Gets the array size.
    inline size_type size() const {return _cols*_rows;}

  public: /** Dynamic memory management*/

    /// Clears the array. Do not free the memory.
    void clear();
    
    /// Releases the memory
    void release();

    /// Relays the mem to caller, detach from it
    void relay_mem(T*& data, E_Int& rows, E_Int& cols, bool& calloc);
    
    /// Gets the capacity
    inline size_type capacity() { return _allocated_sz;}

    /// Expands the capacity.
    int reserve(size_type rows, size_type cols);

    /// Resizes the array preserving the relevant data and set the missing ones to the input val.
    void resize(size_type rows, size_type cols, const value_type* val = 0);
    void resize(size_type rows, size_type cols, const value_type& val);

    /// Appends an array at the end of the array. array must have the same nb of rows.
    int pushBack(const self_type& a);

    /// pushBack a colum defined by values between begin and end.
    template <typename Iterator> int pushBack(Iterator begin, Iterator end);
    
    /// pushBack a std::vector for a one-row DynArray
    int pushBack(const std::vector<T>& a);
    
    ///
    int pushBack(const self_type& a, const E_Int* fields, E_Int sz);

    /// pushBack a constant col colum .
    int pushBack(value_type v, E_Int rows);
    
    /// Extract a selection of entries given by ids from arr.
    void append_selection (const self_type& arr, const std::vector<E_Int>& ids);

  public: /** Iterators */

    /// Returns an iterator at the beginning of data.
    inline iterator begin() {return _data;}

    /// Returns an iterator at the beginning of data (const version).
    inline const_iterator begin() const {return _data;}

    /// Returns an iterator at the end of data.
    inline iterator end() {return _data + (_cols*_rowsMax);}

    /// Returns an iterator at the end of data (const version).
    inline const_iterator end() const {return _data + (_cols*_rowsMax);}

    /// Returns an iterator pointing at the beginning of the i-th column.
    inline iterator col(size_type i)
    {assert(i < _cols); return _data + (i*_rowsMax);}
    // synonym
    inline const_iterator get(size_type i) const { return col(i); }

    /// Returns an iterator pointing at the beginning of the i-th column (const version).
    inline const_iterator col(size_type i) const
    {assert(i < _cols);return _data + (i*_rowsMax);}
    // synonym
    inline iterator get(size_type i) { return col(i); }
    
    /// Returns an iterator pointing at the beginning of the i-th column.
    inline iterator begin(size_type i)
    {assert(i < _cols); return _data + (i*_rowsMax);}

    /// Returns an iterator pointing at the beginning of the i-th column (const version).
    inline const_iterator begin(size_type i) const
    {assert(i < _cols);return _data + (i*_rowsMax);} 
    
    inline size_type stride(int i) { return _rows; }

  public: /** Operators and operations*/

    /// Assignment operator
    self_type& operator=(const self_type& rhs);
    /// Move assignment operator
    self_type& operator=(self_type&& rhs);
    
     /// Assignment operator
    template <typename FldArrayType> self_type& operator=(const FldArrayType& rhs) { *this = DynArray<T>(rhs); return *this;}

    /// Returns a reference to the entry at i-th row and j-th column.
    inline value_type& operator()(size_type i, size_type j){assert((i<_rows)&&(j<_cols));return *(_data + j*_rowsMax + i);}

    /// Returns a reference to the entry at i-th row and j-th column (const version).
    inline const value_type& operator()(size_type i, size_type j) const {assert((i<_rows)&&(j<_cols));return *(_data + j*_rowsMax + i);}

    /// Returns the i-th element of the underneath data array. Equal to view array when the array is a column/vector.
    inline const value_type& operator[](size_type i) const { assert (i < _rows*_cols); return _data[i];}

    /// Returns the i-th element of the underneath data array (write mode).
    /**Equal to view array when the array is a column/vector.*/
    inline value_type& operator[](size_type i) { assert (i < _rows*_cols); return _data[i];}

    /// Adds two DynArrays. Do the job on the common dimensions.
    self_type operator+(const self_type& rhs) const;

    /// Multiplies two DynArrays. Do the job on the common dimensions.
    const self_type operator*(const self_type& rhs) const;
    
    /// 
    self_type& operator*=(T factor) ;
    
    /// Check equality
    bool operator==(const self_type& rhs) const;

    /// Transposes the array.
    self_type& transpose();

    /// Returns the transposed array.
    void transpose(self_type&) const;

    /// Inverse
    static void inverse2(self_type&);
    static void inverse3(self_type&);

    /// Compacts the array using the input flag vector.
    /** Any column with a flag to false(true) is removed(kept). */
    template <typename Vector, typename Vector2>
    static size_type compact (self_type&, const Vector& flag, Vector2& new_Ids);

    template <typename Vector>
    static size_type compact (self_type& a, const Vector& new_Ids);

    template <typename Vector>
    static void changeIndices (self_type& a, const Vector& new_Ids);

    template <typename Vector>
    void uniqueVals(Vector& vals) const;

    inline void uniqueVals(std::set<E_Int>& vals) const;

    void shift(value_type val, size_type istart=0);

  private:

    /// Allocate memory (_data).
    iterator __create(size_type size);

    /// Deallocate memory (_data).
    void __destroy();

    /// Copies the entries from begin to end at location 'there'.
    template <class InputIterator>
    inline static void __copy (InputIterator begin, InputIterator end, iterator there);

    inline static void __copy(value_type v, E_Int rows, iterator there);//fill a column with a given val
    
    ///
    inline void __mappingcopy(const self_type& src_arr, iterator tgt_data, E_Int tgt_stride, E_Int nb_fields, E_Int nb_seq);
    inline void __mappingcopyandfill(const self_type& arr, iterator there, E_Int tgt_stride, E_Int nb_fields, E_Int nb_seq, const T* val);

    /// Assigns v for all the entries between begin and end.
    void __assign (iterator begin, iterator end, const value_type& v);

    /// Reads in an FldArray and sets the data for the DynArray.
    template <typename FldArrayType>
    void __importFldArray (const FldArrayType& i, value_type shift = value_type(0.));
    ///
    template <typename FldArrayType>
    void __importFldArray (const FldArrayType& i, const E_Float** pos, E_Int npos);

    /// Writes out an FldArray from a DynArray.
    template <typename FldArrayType>
    void __exportFldArray (FldArrayType& i, value_type shift = value_type(0.)) const ;

  private:
    /// Total allocated size
    size_type _allocated_sz;
    /// Number of used rows.
    size_type  _rows;
    /// Number of used columns.
    size_type  _cols;
    /// Number of allocated rows.
    size_type  _rowsMax;
    /// Number of allocated columns.
    size_type  _colsMax;
    /// The data storage
    iterator   _data;
    /// allocation style (0 : CPP, 1 : C)
    bool       _calloc;

  }; // End class DynArray

  //----------------------------------------------------------------------------

  typedef DynArray<E_Float>  FloatArray;
  typedef DynArray<E_Int>    IntArray;

  //---------------------------------------------------------------------------

  // Implementation : DynArray

  /// Default constructor.
  template <typename T>
  DynArray<T>::DynArray()
    :_allocated_sz(0), _rows(0), _cols(0), _rowsMax (0), _colsMax (0), _data(0), _calloc(false){}

  /// Constructor with as specified allocation option (C or CPP).
  template <typename T>
  DynArray<T>::DynArray(bool calloc) // false => CPP, true => C
    : _allocated_sz(0), _rows(0), _cols(0), _rowsMax(0), _colsMax(0), _data(0), _calloc(calloc) {}

  /// Constructor with a specified size.
  template <typename T>
  DynArray<T>::DynArray(size_type rows, size_type cols, bool calloc)
    :_allocated_sz(0), _rows(rows), _cols(cols), _rowsMax(rows), _colsMax(cols), _data(0), _calloc(calloc)
  {  
    if (_cols * _rows == 0) // No data so make the attributes consistent.
      _allocated_sz = _cols = _rows = _colsMax = _rowsMax = 0;
    else
      _data = __create(rows*cols);        
  }

  /// Constructor with a specified size.
  template <typename T>
  DynArray<T>::DynArray(size_type rows, size_type cols, const value_type& val, bool calloc)
    :_allocated_sz(0), _rows(0), _cols(0), _rowsMax (0), _colsMax (0), _data(0), _calloc(calloc)
  {
    resize (rows, cols, val);
  }

  /// Constructor with input data : WARNING becomes owner => set input point to null
  template <typename T>
  DynArray<T>::DynArray(T*& compactdata, size_type rows, size_type cols, bool calloc)
  {
    _data = compactdata;
    compactdata = nullptr;

    _cols = _colsMax = cols;
    _rows = _rowsMax = rows;
    _calloc = calloc;
    _allocated_sz = cols * rows;
  }

  /// Copy constructor.
  template <typename T>
  template <typename U>
  DynArray<T>::DynArray(const DynArray<U>& rhs)
    :_allocated_sz(0), _rows(rhs._rows), _cols(rhs._cols), _rowsMax (rhs._rows),_colsMax (rhs._cols), _data(0), _calloc(rhs._calloc)
  {
    if (_cols * _rows == 0) // No data so make the attributes consistent.
      _allocated_sz = _cols = _rows = _colsMax = _rowsMax = 0;
    else // Assign the data.
    {
      _data = __create(_rows * _cols);
      __mappingcopy(rhs, _data, _rowsMax, _rows, _cols);
    }
  }

  ///
  template <typename T>
  DynArray<T>::DynArray(const DynArray<T>& rhs)
    :_allocated_sz(0), _rows(rhs._rows), _cols(rhs._cols), _rowsMax(rhs._rows),_colsMax(rhs._cols), _data(0), _calloc(rhs._calloc)
  {
    if (_cols * _rows == 0) // No data so make the attributes conssistent.
      _allocated_sz = _cols = _rows = _colsMax = _rowsMax = 0;
    else // Assign the data.
    {
      _data = __create(_rows * _cols);
      __mappingcopy(rhs, _data, _rowsMax, _rows, _cols);
    }
  }

  template <typename T>
  DynArray<T>::DynArray(DynArray<T>&& rhs)
    :_allocated_sz(rhs._allocated_sz), _rows(rhs._rows), _cols(rhs._cols), _rowsMax(rhs._rowsMax), _colsMax(rhs._colsMax), _data(rhs._data), _calloc(rhs._calloc)
  {
    rhs._data = nullptr;
    rhs._rows = rhs._rowsMax = rhs._cols = rhs._colsMax = rhs._allocated_sz = 0;
  }

  /** Conversion specialization */
  ///FldArray --> DynArray
  template <typename T> inline
    DynArray<T>::DynArray(const FldArray<T>& a, T shift)
    :_allocated_sz(0), _rows(0), _cols(0), _rowsMax (0), _colsMax (0),_data(0), _calloc(false)
  {__importFldArray(a, shift);}

  ///DynArray --> FldArray
  template <typename T> inline
  void
  DynArray<T>::convert(FldArray<T>& out, T shift) const
  {__exportFldArray(out, shift);}

  template <typename T> inline
  void
  DynArray<T>::extract_field(size_type i, std::vector<T>& f)
  {
    f.resize(_cols, {});
    for (size_type j = 0; j < _cols; ++j)
      f[j] = *(_data + j * _rowsMax + i);
  }
  
  template <typename T> inline
  void
  DynArray<T>::set_field(size_type i, const std::vector<T>& f)
  {
    assert(f.size() == (size_t)_cols);
    assert(i < _rows);

    for (size_type j = 0; j < _cols; ++j)
      *(_data + j * _rowsMax + i) = f[j];
  }

  ///Clear the array. Do not deallocate.
  template <typename T>
  void DynArray<T>::clear(){_rows = _cols = 0;}
  
  /// Deallocate memory
  template <typename T>
  void DynArray<T>::release(){__destroy(); _allocated_sz = _cols = _rows = _colsMax = _rowsMax = 0;}

  /// Pass memory responsability to caller, and detach from it
  template <typename T>
  void DynArray<T>::relay_mem(T*& data, E_Int& rows, E_Int& cols, bool& calloc)
  {
    data = _data;
    rows = _rows;
    cols = _cols;
    calloc = _calloc;

    _data = nullptr;
    _allocated_sz = _cols = _rows = _colsMax = _rowsMax = 0;
  }

  /** Expands the memory with cols and rows.
  *   Memory is reallocated (preserve data). 
  */
  template <typename T>
  int
  DynArray<T>::reserve(size_type rows, size_type cols)
  {
    if ((rows == _rowsMax) && (cols == _colsMax)) // Allocated memory is OK.
      return 0;

    if (cols < 0)
    {
      std::cout << "DynArray error: ask to reserve with a negative value => int32 limit ?" << std::endl;
      return 1;
    }

    size_type r = std::max(rows, _rowsMax);
    size_type c = std::max(cols, _colsMax);

    iterator newdata(__create(r*c)), start(0), there(0);

    if (newdata == nullptr) return 1;

    // Assign the values column by colum.
    for (size_type i = 0; i < _cols; ++i)
    {
      start = _data + i*_rowsMax;
      there = newdata + i*r;
      __copy (start, start + _rows, there);
    }

    __destroy();
    _data = newdata;
    _rowsMax = r;
    _colsMax = c;

    return 0;
  }

  template <typename T>
  void
  DynArray<T>::resize(size_type rows, size_type cols, const value_type* val)
  {
    size_type size(rows*cols);
    size_type cols0(_cols);
      
    if (cols == _cols && rows <= _rows)
    {
      _rows = rows;
      return;
    }
    if ((rows != _rows) || (cols > _colsMax)) //need to reallocate
    {
      iterator newdata(__create(size));
       
      // reassign the values column by colum.
      if (_data != 0)
      {
        __mappingcopyandfill(*this, newdata, rows, _rows, _cols, val);
        __destroy();
      }
        
      _rowsMax = _rows = rows;
      _colsMax = _cols = cols;
      _data = newdata;
    }
    if (val && cols > cols0)
      __assign (_data + cols0*rows, _data + size, *val); // Fill the data's tail with val
      
    _rows = rows;
    _cols = cols;
  }
  
  template <typename T>
  void
  DynArray<T>::resize(size_type rows, size_type cols, const value_type& val)
  {
    size_type size(rows*cols);
    size_type cols0(_cols);
      
    if ((rows != _rows) || (cols > _colsMax)) //need to reallocate
    {
      iterator newdata(__create(size));
      
      if (_data != 0)
      {
        // Assign the values column by colum.
        __mappingcopyandfill(*this, newdata, rows, _rows, _cols, &val);
        __destroy();
      }
      
      _rowsMax = _rows = rows;
      _colsMax = _cols = cols;
      _data = newdata;
    }
    if (cols > cols0)
      __assign (_data + cols0*rows, _data + size, val); // Fill the data's tail with val
    
    _rows = rows;
    _cols = cols;
  }

  ///
  template <typename T>
  int
  DynArray<T>::pushBack(const self_type& a)
  {
    if (_rows == 0)// Empty array
    {
      *this = a;
      return 0;
    }
    else if (_rows != a._rows) return 1;

    if (a._cols > (_colsMax - _cols))
    {
      // If the spare room is too tight.
      int err = reserve(_rows, 2*(_cols+a._cols));  // Double the columns capacity.
      if (err) return 1;
    }
    
    E_Int last = _cols;
    _cols += a._cols;
    __mappingcopy(a, _data+last*_rowsMax, _rowsMax, a._rows, a._cols);

    return 0;
  }
  
  ///
  template <typename T>
  int
  DynArray<T>::pushBack(const self_type& a, const E_Int* fields, E_Int sz)
  {
    if (sz==0) return 0;
    
    assert (_rows==0 || sz == _rows);
    
    if (_rows == 0) // Empty array
      _rows=sz;
    
    if (a._cols > (_colsMax - _cols)) 
    {    
      // If the spare room is too tight.
      int err = reserve(_rows, 2*(_cols+a._cols));  // Double the columns capacity.
      if (err) return 1;
    }

    for (size_type j = _cols; j < _cols + a._cols; ++j){
      for (size_type i = 0; i < sz; ++i)
        *(_data + j*_rowsMax +i) = a(fields[i], j-_cols);
    }
          
    _cols += a._cols;

    return 0;
  }
  
  ///
  template <typename T>
  int
  DynArray<T>::pushBack(const std::vector<T>& a)
  {
    if (a.empty())  return 0;
    if (_rows == 0) _rows=1;
    if (_rows > 1)  return 1; //shape inconsistency

    if (E_Int(a.size()) > (_colsMax - _cols)) // If the spare room is too tight.
    {
      int err = reserve(_rows, 2*(_cols+a.size()));  // Double the columns capacity.
      if (err) return 1;
    }

    assert(_data != nullptr);

    iterator there = _data + (_cols*_rowsMax);//col(_cols);
    __copy(a.begin(), a.end(), there);
    
    _cols += a.size();

    return 0;
  }

  ///
  template <typename T>
  template <typename Iterator>
  int
  DynArray<T>::pushBack(Iterator begin, Iterator end)
  {
    if (_rows == 0)// Empty array
      _rows = (end - begin);
    else if (_rows != (end - begin))
      return 1;

    if (_cols == _colsMax)
    {
      int err = reserve (_rows, 2*(_cols+1));
      if (err) return 1;
    }

    iterator there(col(_cols++));
    __copy(begin, end, there);

    return 0;
  }

  ///
  template <typename T>
  int
  DynArray<T>::pushBack(value_type val, E_Int rows) 
  {      
    if (_rows == 0)// Empty array
      _rows = rows;
    else if (_rows != rows)
      return 1; // shape inconsistency

    if (_cols == _colsMax)
    {
      int err = reserve(_rows, 2 * (_cols + 1));
      if (err) return 1;
    }

    iterator there(col(_cols++));
    __copy(val, rows, there);

    return 0;
  }
  
  ///
  template <typename T>
  void
  DynArray<T>::append_selection(const self_type& arr, const std::vector<E_Int>& ids)
  {
    assert ((_rows==0) || (_rows==arr._rows));
    
    for (size_t i = 0; i < ids.size(); ++i)
      this->pushBack(arr.col(ids[i]), arr.col(ids[i])+arr._rows);
  }

  /// Assignment operator.
  template <typename T>
  DynArray<T>&
  DynArray<T>::operator=(const self_type& rhs)
  {
    if (this == &rhs) //avoid self assignment.
      return *this;
      
    E_Int required_sz = rhs._cols*rhs._rows;
    if (required_sz == 0)
    {
      this->clear();
      return *this;
    }
      
    if (_allocated_sz < required_sz)
    {
      __destroy(); // Delete old memory.
      _calloc = rhs._calloc; //use same policy as rhs
      _data = __create(required_sz); // Reallocate.
      _rowsMax = _rows = rhs._rows;
      _colsMax = _cols = rhs._cols;
      _allocated_sz = required_sz;
    }
    else
    {
      // WARNING : compacted assignment
      _rowsMax = _rows = rhs._rows;
      assert (_rowsMax > 0);
      _colsMax = _allocated_sz / _rowsMax;
      _rows = rhs._rows;
      _cols = rhs._cols;
      assert (_cols <= _colsMax);
    }
      
    // Assign the rhs data.
    __mappingcopy(rhs, _data, _rowsMax, _rows, _cols);

    return *this;
  }

  /// Move Assignment operator.
  template <typename T>
  DynArray<T>&
  DynArray<T>::operator=(self_type&& rhs) 
  {
    if (this == &rhs) //avoid self move.
      return *this;

    __destroy(); // delete _data

    _data         = rhs._data;
    _calloc       = rhs._calloc;
    _rows         = rhs._rows;
    _rowsMax      = rhs._rowsMax;
    _cols         = rhs._cols;
    _colsMax      = rhs._colsMax;
    _allocated_sz = rhs._allocated_sz;

    rhs._data = nullptr;
    rhs._rows = rhs._rowsMax = rhs._cols = rhs._colsMax = rhs._allocated_sz = 0;

    return *this;
  }

  /// Adds two DynArrays. Do the job on the common dimensions.
  template <typename T>
  DynArray<T>
  DynArray<T>::operator+(const self_type& rhs) const 
  {
    size_type r = std::min(_rows, rhs._rows);
    size_type c = std::min(_cols, rhs._cols);
    DynArray<T> result(r,c);

    for (size_type j = 0; j < c; ++j)
      for (size_type i = 0; i < r; ++i)
        result(i,j) = (*this)(i,j)+ rhs(i,j);

    return result;
  }

  /// Multiplies two DynArrays. Do the job on the common dimensions.
  template <typename T>
  const DynArray<T>
  DynArray<T>::operator*(const self_type& rhs) const 
  {
    size_type n = std::min(_rows, rhs._cols);
    size_type m = std::min(_cols, rhs._rows);

    DynArray<T> result(n,n, value_type(0.));//fixme
      
    for (size_type j = 0; j < n; ++j)
      for (size_type i = 0; i < n; ++i)
        for (size_type k = 0; k < m; ++k)
          result(i,j) += (*this)(i,k) * rhs(k,j);

    return result;
  }
  
  /// Multiplies two DynArrays. Do the job on the common dimensions.
  template <typename T>
  DynArray<T>& DynArray<T>::operator*=(T factor) 
  {    
    for (size_type i = 0; i < _rows; ++i)
      for (size_type j = 0; j < _cols; ++j)
          (*this)(i,j) *= factor;

    return *this;
  }
  
  /// Checks equality
  template <typename T>
  bool
  DynArray<T>::operator==(const self_type& rhs) const
  {
    size_t r=rows();
    size_t c=cols();
    
    if (rhs.cols() != c)
      return false;
    if (rhs.rows() != r)
      return false;
    
    for (size_t j=0; j<c; ++j)
    {
      for (size_t i=0; i<r; ++i)
      {
        if ((*this)(i,j) != rhs(i,j))
          return false;
      }
    } 
    return true;  
  }

  /// Transposes the input array.
  template <typename T>
  DynArray<T>&
  DynArray<T>::transpose()
  {
    DynArray<T> result(_cols, _rows);

    for (size_type j = 0; j < _rows; ++j)
      for (size_type i = 0; i < _cols; ++i)
        result(i,j) = (*this)(j,i);

    *(this) = result;
    return *this;
  }

  template <typename T>
  void
  DynArray<T>::transpose(self_type& result) const 
  {
    result.resize(_cols, _rows);

    for (size_type j = 0; j < _rows; ++j)
      for (size_type i = 0; i < _cols; ++i)
        result(i,j) = (*this)(j,i);
  }

  template <typename T>
  template <typename Vector1, typename Vector2>
  typename DynArray<T>::size_type
  DynArray<T>::compact (self_type& a, const Vector1& keep, Vector2& new_Ids)
  {
    size_type cols(a._cols);
    // Fast returns
    if (cols == 0) return 0;

    assert (keep.size() == (size_t)cols);

    bool carry_on(false);
    size_type i1(0), i2(cols-1);

    new_Ids.clear();
    new_Ids.resize(cols, IDX_NONE);

    do 
    {
      while ((i1 <= i2) && keep[i1]){new_Ids[i1] = i1; ++i1;}  // Get the first empty column.
      while ((i1 <= i2) && !keep[i2]){--i2;} // Get the first column to move from the tail.

      carry_on = (i1 < i2);

      if (carry_on)
      { // Copy column i2 in column i1.
        new_Ids[i2] = i1;
        const_iterator  it2 = a.col(i2--);
        iterator        it1 = a.col(i1++);
        __copy(it2, it2 + a._rows, it1);
      }
    }
    while (carry_on);

    a.resize (a._rows, i1);

    return (cols - i1);
  }

  template <typename T>
  template <typename Vector>
  typename DynArray<T>::size_type
  DynArray<T>::compact (self_type& a, const Vector& new_Ids)
  {
    size_type       cols(a.cols()), newId, count(0);
    const_iterator  it2;
    iterator        it1;

    assert (new_Ids.size() == (size_t)cols);

    for (size_type i = 0; i < cols; ++i)
    {
      newId = new_Ids[i];
      if (newId == IDX_NONE)
        ++count;
      else if (newId != i)
      {
        it2 = a.col(i);
        it1 = a.col(newId);
        __copy (it2, it2 + a._rows, it1);
      }
    }
    a.resize(a._rows, cols-count);

    return (cols - a.cols());
  }

  ///
  template <>
  template <typename Vector>
  void
  DynArray<E_Int>::changeIndices (self_type& a, const Vector& new_Ids)
  {
    for (size_type i=0; i<a._cols; ++i)
    {
      for (size_type j=0; j<a._rows; ++j)
      {
        E_Int& v = *(a._data + i*a._rowsMax + j);
        if (v != IDX_NONE)
          v = new_Ids[v];
      }
    }
  }

  ///
  template <typename T>
  template <typename Vector>
  void
  DynArray<T>::uniqueVals(Vector& vals) const
  {
    vals.clear();
    std::set<T> pool;
    T v;
    
    for (size_type i=0; i<_cols; ++i)
    {
      for (size_type j=0; j<_rows; ++j)
      {
        v = *(_data + i*_rowsMax + j);
        if (pool.insert(v).second)
          vals.push_back(v);
      }
    }
  }
  
  ///
  template <>
  template <typename Vector>
  void
  DynArray<E_Int>::uniqueVals(Vector& vals) const
  {
    vals.clear();
    std::set<E_Int> pool;
    E_Int v;
    
    pool.insert(IDX_NONE);//to avoid to extract none indices
    
    for (size_type i=0; i<_cols; ++i)
    {
      for (size_type j=0; j<_rows; ++j)
      {
        v = *(_data + i*_rowsMax + j);
        if (pool.insert(v).second)
          vals.push_back(v);
      }
    }
  }

  ///
  template <> inline
  void
  DynArray<E_Int>::uniqueVals(std::set<E_Int>& vals) const
  {
    vals.clear();
    E_Int v;

    for (size_type i = 0; i<_cols; ++i)
    {
      for (size_type j = 0; j<_rows; ++j)
      {
        v = *(_data + i * _rowsMax + j);
        if (v != IDX_NONE) vals.insert(v);
      }
    }
  }

  template <typename T>
  void
  DynArray<T>::shift(T val, size_type istart)
  {
    size_type beg=istart*_rowsMax;
    T* p = _data + beg;
    size_type end = _rowsMax*_cols;
    for (size_type i = beg; i < end; ++i)
    {
      if (*p != IDX_NONE)
        *(p++) += val;
    }
  }

  /*****************  private ***********************************/

  template <typename T>
  T*
    DynArray<T>::__create(size_type size)
    {
      if (size == 0)
        return nullptr;

      T* p = nullptr;
      if (!_calloc)
        p = NUGA::allocator<false>::allocate<T>(size);
      else
        p = NUGA::allocator<true>::allocate<T>(size);

      if (p != nullptr) _allocated_sz = size;
      assert(p != nullptr);
      return p;
  }

  template <typename T>
  void
    DynArray<T>::__destroy(){
      
    if (!_calloc)
      NUGA::allocator<false>::deallocate(_data);
    else
      NUGA::allocator<true>::deallocate(_data);
  }

  template <typename T>
  template <typename InputIterator>
  void
  DynArray<T>::__copy(InputIterator begin, InputIterator end, iterator there){

    assert(there != nullptr);

    for (InputIterator it = begin; it != end; ++it)
      *(there++) = *(it);

    //std::copy (begin, end, there);
  }

  template <typename T>
  void
    DynArray<T>::__copy(value_type val, E_Int rows, iterator there) {

    assert(there != nullptr);
    
    for (E_Int i=0; i < rows; ++i)
      *(there++) = val;
  }
  
  template <typename T>
  inline void
  DynArray<T>::__mappingcopy(const self_type& arr, iterator tgt_data, E_Int tgt_stride, E_Int nb_fields, E_Int nb_entries)
  {  
    if (nb_fields*nb_entries == 0)
      return;

    assert (nb_fields >=0);
    assert (nb_entries >=0);
    assert (nb_fields <= arr._rows);               // the source is a subpart of arr
    assert (nb_entries <= arr._cols);              // the source is a subpart of arr
    assert (arr._rows*arr._cols <= _allocated_sz); //required space is enough
    assert (tgt_stride >= arr._rows);                //shape is adapted
   
    if (tgt_stride == arr._rowsMax && nb_fields == arr._rows) //same strides and plain copy-> flat copy instead
    {
      __copy(arr.begin(), arr.end(), tgt_data);
      return;
    }

    T *start, *there;
    // Assign the values column by colum.
    for (size_type i = 0; i < nb_entries; ++i)
    {
      start = arr._data + i*arr._rowsMax;
      there = tgt_data + i*tgt_stride;
      __copy(start, start + nb_fields, there);
    }
  }

  template <typename T>
  inline void
  DynArray<T>::__mappingcopyandfill(const self_type& arr, iterator tgt_data, E_Int tgt_stride, E_Int nb_fields, E_Int nb_entries, const T* val)
  {
    if (nb_fields*nb_entries == 0)
      return;

    assert (nb_fields >=0);
    assert (nb_entries >=0);
    assert (nb_fields <= arr._rows);               // the source is a subpart of arr
    assert (nb_entries <= arr._cols);              // the source is a subpart of arr
    assert (arr._rows*arr._cols <= _allocated_sz); //required space is enough
    assert (tgt_stride >= arr._rows);                //shape is adapted
   
    if (!val || (tgt_stride == nb_fields)) //no room for filling or no value specified.
    {
      __mappingcopy(arr, tgt_data, tgt_stride, nb_fields, nb_entries);
      return;
    }

    T *start, *there;
    // Assign the values column by colum.
    for (size_type i = 0; i < nb_entries; ++i)
    {
      start = arr._data + i*arr._rowsMax;
      there = tgt_data + i*tgt_stride;
      __copy(start, start + nb_fields, there);
      __assign(there+nb_fields, there+tgt_stride, *val);
    }
  }

  template <typename T>
  void
    DynArray<T>::__assign (iterator begin, iterator end, const value_type& v){

      for (iterator it = begin; it != end; ++it)
        *it = v;

      //std::fill_n (begin, end, v);
  }

  template <typename T>
  template <typename U>
  void
    DynArray<T>::__importFldArray (const U& a, value_type shift){

      size_type rows(a.getNfld()), cols(a.getSize());
      const T* *p =  new const T*[rows];

      reserve(rows, cols);

      for (size_type i = 0; i < rows; ++i)
        p[i] = a.begin(i+1);

      for (size_type j = 0; j < cols; ++j){
        for (size_type i = 0; i < rows; ++i){
          *(_data + j*rows +i) = *(p[i]+j) + shift;
        }
      }

      delete[] p;

      _rows = _rowsMax = rows;
      _cols = _colsMax = cols;
  }
  
  template <>
  template <typename U>
  void
    DynArray<E_Float>::__importFldArray (const U& a, const E_Float** p, E_Int npos){

      size_type rows(a.getNfld()), cols(a.getSize());
            
      assert (npos <= rows);
      
      reserve(npos, cols);

      for (size_type j = 0; j < cols; ++j){
        for (size_type i = 0; i < npos; ++i){
          *(_data + j*rows +i) = *(p[i]+j);
        }
      }

      delete[] p;

      _rows = _rowsMax = rows;
      _cols = _colsMax = cols;
  }
  
  ///FldArrayF --> FloatArray
  template <> inline
    DynArray<E_Float>::DynArray(const FldArray<E_Float>& a, E_Int posx, E_Int posy, E_Int posz)
    :_allocated_sz(0), _rows(0), _cols(0), _rowsMax (0), _colsMax (0),_data(0), _calloc(false)
  {
    E_Int npos = (posz != -1) ? 3 : 2;
    const E_Float* *p =  new const E_Float*[npos];
    
    p[0] = a.begin(posx+1);
    p[1] = a.begin(posy+1);
    if (npos == 3) p[2] = a.begin(posz+1);
    
    __importFldArray(a, p, npos);
  }

  template <typename T>
  template <typename FldArrayType>
  void
    DynArray<T>::__exportFldArray (FldArrayType& a, value_type shift) const {

      FldArrayType out;
      out.resize(_cols, _rows);

      T* *p =  new T*[_rows];

      for (size_type i = 0; i < _rows; ++i)
        p[i] = (T*)out.begin(i+1);

      for (size_type j = 0; j < _cols; ++j){
        for (size_type i = 0; i < _rows; ++i){
          *(p[i]+j) = *(_data + j*_rows +i) + shift;
        }
      }

      delete[] p;

      a = out;
  }

  template <typename T>
  std::ostream & operator<<(std::ostream& out, const DynArray<T>& arr){

    out << "####################################" << std::endl;

    // Print out the matrix.
    for (E_Int i = 0; i < arr.rows(); ++i){
      for (E_Int j = 0; j < arr.cols(); ++j)
        out << ((arr(i,j) == IDX_NONE) ? -1 : arr(i,j)) << " ";
      out << std::endl;
    }

    out << std::endl;

    //Print out the attibutes.
    out << "rows: " << arr._rows << "(" << arr._rowsMax << ")" << std::endl;
    out << "cols: " << arr._cols << "(" << arr._colsMax << ")" << std::endl;
    out << std::endl;

    // Print out the raw data.
    //for (E_Int i = 0; i < array._rowsMax*array._colsMax; ++i)
    //  out << array._data[i] << " ";
    //out << std::endl;

    out << "####################################" << std::endl;

    return out;
  }

  template <typename T> inline
  void DynArray<T>::inverse2(self_type& A) 
  {
    E_Float detA = A(0,0)*A(1,1) - A(1,0)*A(0,1);
    detA = 1./detA;

    E_Float a0 =   A(1,1) * detA;
    
    A(0,1)     = - A(0,1) * detA;
    A(1,0)     = - A(1,0) * detA;
    A(1,1)     =   A(0,0) * detA;
    A(0,0)     =   a0;
  }

  template <typename T> inline
  void DynArray<T>::inverse3(self_type& A) 
  {
    E_Float a = A(0,0);
    E_Float b = A(0,1);
    E_Float c = A(0,2);
    E_Float d = A(1,0);
    E_Float e = A(1,1);
    E_Float f = A(1,2);
    E_Float g = A(2,0);
    E_Float h = A(2,1);
    E_Float i = A(2,2);

    E_Float detA = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    detA = 1./detA;

    A(0,0) = (e*i - f*h) * detA;
    A(0,1) = (c*h - b*i) * detA;
    A(0,2) = (b*f - c*e) * detA;
    A(1,0) = (f*g - d*i) * detA;
    A(1,1) = (a*i - c*g) * detA;
    A(1,2) = (c*d - a*f) * detA;
    A(2,0) = (d*h - e*g) * detA;
    A(2,1) = (b*g - a*h) * detA;
    A(2,2) = (a*e - b*d) * detA;
  }

} // namespace K_FLD

#endif // Header

