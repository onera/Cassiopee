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
#ifndef _KCORE_FLDARRAY_H_
#define _KCORE_FLDARRAY_H_

#include "Def/DefTypes.h"
#include <assert.h>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>

#include "Fld/FldFortranVecP.h"
#include "Def/DefFunction.h"
#include "parallel.h"

#define TEMPLATE_T template<typename T>
#define SPECIALISE template<>

// Triggers or not copy by memcpy instead of setallvaluesf
#define MEMCOPY

extern "C" void k6fldcopyfrom_(const E_Int& deb,
                               const E_Int& sizerhs,
                               const E_Int& sizelhs,
                               const E_Int& nfld,
                               const E_Float* rhs,
                               E_Float* lhs);
extern "C" void k6fldintcopyfrom_(const E_Int& deb,
                                  const E_Int& sizerhs,
                                  const E_Int& sizelhs,
                                  const E_Int& nfld,
                                  const E_Int* rhs,
                                  E_Int* lhs);
extern "C" void k6setallbvaluesatf_(const E_Int& sizeLoc,
                                    const E_Int& nfldLoc,
                                    E_Float* _data,
                                    const E_Float& val,
                                    const E_Int& ni,
                                    const E_Int& nj,
                                    const E_Int& nk,
                                    const E_Int& border);

// Set rake for a compact array
#define SETRAKE                                                         \
  if (_stride == 1)                                                     \
  { for (E_Int i = 0; i < _nfldLoc; i++) _rake[i] = _data+i*_sizeLoc; } \
  else                                                              \
  { for (E_Int i = 0; i < _nfldLoc; i++) _rake[i] = _data+i; }

namespace K_FLD
{
#define NUMFIELD0 1

// ============================================================================
// @Name FldArray
// @Memo Array of basic type elements (int, float,...)
/* @Text

> Caution: indices for elements go from 0 to size-1
>          indices for number of fields go from NUMFIELD0 to NUMFIELD0+nfld-1

*/
// ============================================================================
TEMPLATE_T class ArrayAccessor; //forward dec for frienship (gcc)

TEMPLATE_T
class FldArray
{
  public:
    template<typename U> friend class ArrayAccessor;
  public:
    // Typedefs for iterators
    /* Template type. */
    typedef T value_type;

    ///+ 1- Constructors / Destructor
    /** Default Constructor */
    FldArray();
    /** Constructor with number of elements (size) and
        number of field (nfld) */
    explicit FldArray(E_Int size, E_Int nfld=NUMFIELD0,
                      E_Boolean compact=true, E_Boolean fortranOrdered=true);
    /** Copy constructor. We take the size of
        the right hand FldArrayI to build the new FldArray.
        Constructed object has new (duplicated) datas */
    FldArray(const FldArray& rhs);
    /** Copy specified fields of a FldArray.
        begin: first field copied,
        end  : last  field copied */
    FldArray(const FldArray& rhs, E_Int begin, E_Int end);
    /** Constructor from T*.
        if shared=false, data are copied from listofvalues
        if shared=true, data is just a pointer on listofvalues,
        _shared attribute is set. In this case, FldArray is just a handle
        on data.
        Forcement compact. */
    FldArray(E_Int size, E_Int nfld,
             const value_type* listofvalues, E_Boolean shared=false,
             E_Boolean fortranOrdered=true);
    /** Constructor from T**.
        if shared=false, data are copied from listofvalues
        if shared=true, data is just a pointer on listofvalues,
        _shared attribute is set. In this case, FldArray is just a handle
        on data.
        Forcement non compact. */
    FldArray(E_Int size, E_Int nfld,
             value_type** listofvalues, E_Boolean shared=false,
             E_Boolean fortranOrdered=true);
    /** Constructor for NGON connectivity (shared, array1). */
    FldArray(E_Int nfaces, E_Int nelts, E_Int sizeNGon, E_Int sizeNFace,
             value_type* ngonCompact);
    /** Constructor for NGON connectivity (shared, array2). */
    FldArray(E_Int nfaces, E_Int nelts,
             value_type* ngon, value_type* nface,
             value_type* indPG, value_type* indPH, 
             E_Int sizeNGon, E_Int sizeNFace, value_type* PE=NULL);
    /** Constructor for ME connectivities (shared, array3). */
    FldArray(std::vector< FldArray<T>* >& list);
    
    /** Destructor */
    ~FldArray();
    ///-

    ///+ 2- Memory
    /* Alloc memory */
    void acquireMemory();
    /* Delete memory */
    void releaseMemory(void);
    ///-

    ///+ 3- Operators
    /** Result of = operator has new (duplicated) data. */
    FldArray& operator=  (const FldArray& rhs);
    /** Init a already dimensionned array to a value */
    FldArray& operator=  (T val);
    /** Copy from indice beg the field rhs as a subpart of the
        current array (static interface only). */
    void copyFrom(E_Int beg, const FldArray& rhs);
    /** Print */
    void print();
    ///-

    ///+ 4- Iterators
    /** Starting iterator on a field. Default field is the first one. */
    inline T* begin(E_Int fld=NUMFIELD0);
    /** Starting const T* on a field.
        Default field is the first one. */
    inline const T* begin(E_Int fld=NUMFIELD0) const;
    /** Ending iterator on the last field. */
    T* end();
    /** Ending const T* on the last field. */
    inline const T* end() const;
    /** Ending iterator on a field. Default field is the last one. */
    T* end(E_Int fld);
    /** Ending const T* on a field. Default field is the last one. */
    inline const T* end(E_Int fld) const;
    ///-

    ///+ 5- Accessors
    /** Get size of one field. */
    inline E_Int getSize() const { return _sizeLoc; }
    /** Get number of fields. */
    inline E_Int getNfld() const { return _nfldLoc; }
    /** Get stride */
    inline E_Int getStride() const { return _stride; }
    /** Return true if array is compact */
    inline E_Boolean getCompact() const { return _compact; }
    /** GetApi (1: compact, 2: rake) */
    inline E_Int getApi() { if (_compact == false) return 2; else return 1; }

    /** Get/Set NGon */
    inline E_Int isNGon() const { return _ngon; }
    // 1: compact array1 CGNSv3, 2: rake CGNSv3, 3: rake CGNSv4
    void setNGon(E_Int ngon) { _ngon = ngon; };
    /** Only if  ngon */
    inline E_Int getNFaces();
    inline E_Int getNElts();
    inline E_Int* getNGon();
    inline E_Int* getNFace();
    inline E_Int* getIndPG();
    inline E_Int* getIndPH();
    inline E_Int* getPE();
    inline E_Int getSizeNGon();
    inline E_Int getSizeNFace();
    inline E_Int* getFace(E_Int no, E_Int& size);
    inline E_Int* getFace(E_Int no, E_Int& size, E_Int* ngon, E_Int* indPG); // faster
    inline E_Int* getElt(E_Int no, E_Int& size);
    inline E_Int* getElt(E_Int no, E_Int& size, E_Int* nface, E_Int* indPH); // faster

    /** ME */
    inline size_t getNConnect();
    inline FldArray<T>* getConnect(E_Int i=0);
    
    /** Get l-th element of array of values (0<=l < _nfldLoc*_sizeMax) */
    inline value_type  operator[]  (E_Int l) const;
    /** Get/set l-th element of array of values (0<=l < _nfldLoc*_sizeMax) */
    inline value_type& operator[] (E_Int l);
    /** Get l-th element of fld-th field with operator () */
    inline value_type  operator() (E_Int l, E_Int fld) const;
    /** Access in write mode of the operator (). */
    inline value_type& operator() (E_Int l, E_Int fld);
    /** Index of min value in field nfld. We take 1 element over inc. */
    E_Int indMin(E_Int inc=1, E_Int nfld=NUMFIELD0) const;
    /** Index of max value in field nfld. We take one element over inc. */
    E_Int indMax(E_Int inc=1, E_Int nfld=NUMFIELD0) const;
    ///-

  public:
    ///+ 6- Modifiers
    /** Set all values to those of a given FldArray. */
    void setAllValuesAt(const FldArray& tableofvalues);
    /** Set an uniform value on each field. valuesoffield is an array of size nfld. */
    void setByField(const FldArray& valuesoffields);
    /** Set the value of the numfieldTo field with the numFieldFrom field of fromValArray */
    void setOneField(const FldArray& fromValArray,
                     E_Int numfieldFrom=NUMFIELD0,
                     E_Int numFieldTo=NUMFIELD0);
    /** Set an uniform value on all fields. */
    void setAllValuesAt(value_type val);
    /** Set all values at ZERO. */
    void setAllValuesAtNull();
    /** Set all border values at a given value. The field must be defined on
        interfaces. Ni, nj, nk is the grid number of nodes,
        border is the number of init layers. */

  public:
    ///+ 7- Memory Management
    // STATIC INTERFACE
    /** Resize with new array allocation WITHOUT copy. */
    void malloc(E_Int nsize, E_Int nnfld=NUMFIELD0);
    /** Resize with new array allocation.
     if shared=true, listofvalues is shared
     if shared=false, listofvalues is copied. */
    void malloc(E_Int nsize, E_Int nnfld,
              const T* listofvalues, E_Boolean shared=false);
    /** Resize with new array allocation WITH copy. */
    void reAlloc(E_Int nsize, E_Int nnfld=NUMFIELD0);
    /** Resize with new array allocation WITH copy, keeping
        elements at the same position in the matrix. */
    void reAllocMat(E_Int l, E_Int q);
    /** Resize with new array allocation WITH copy, keeping
        elements at the same position in the matrix. 
        No openMP in this version */
    void reAllocMatSeq(E_Int l, E_Int q);
    /** Resize the array. No new allocation if the new dimension
        if smaller than oldest dimension and WITHOUT copy. */
    void resize(E_Int nsize, E_Int nnfld=NUMFIELD0);
    /** Set the pointer to actual data. The given data pointer is not copied,
        it is used as is shared with the proprietary object. */
    static FldArray<T>* setMemoryCompactZone(value_type* data, E_Int size, E_Int nfld=NUMFIELD0);
    static FldArray<T>* setMemoryNonCompactZone(value_type** data, E_Int size, E_Int nfld=NUMFIELD0);
    // DYNAMIC INTERFACE (only for shared=false, compact=true)
    void reserve(E_Int nfld, E_Int size);
    void resize(E_Int nfld, E_Int size, const T& val);
    void pushBack(const FldArray& a);
    void pushBack(const FldArray& a, E_Int i);
    /// pushBack an entry defined by values between begin and end.
    template <typename Iterator> void pushBack(Iterator begin, Iterator end);
    /// pushBack a std::vector for a one-row FldArray
    void pushBack(const std::vector<value_type>& a);
    /// Set the sizeLoc to 0. Do not clear the memory.
    void clear(){_sizeLoc = _nfldLoc= 0;}
    // copy (static and dynamic interface)
    void copy(const FldArray& rhs, E_Int beg, E_Int end);
    /// Compacts the array using the input flag vector.
    /** Any column with a flag to false(true) is removed(kept). */
    template <typename Vector1, typename Vector2>
    static E_Int compact (FldArray& a, const Vector1& keep, Vector2& new_Ids);

  public:
    ///+ 8- Specific to FldArrayF
    ///
    /*
    void setAllBorderValuesAt(value_type val,
                              E_Int ni, E_Int nj, E_Int nk,
                              E_Int border);
    */
    /* */
    static value_type badValue();
    static value_type Zero();
    /** Array Multiplication (only for small arrays) */
    /* void matMult(FldArray& array); */
    /** Array Multiplication with transposed matrix (only for small arrays) */
    /* void matMultTrans(FldArray& array); */
    /** For all rows, multiplication between matrix coefficients and
        vector coefficients (it is not a true matrix vector product) */
    /* void matMultVec(const FldArray& array); */
    /** Sqrt of each element */
    void sqrt();
    ///-

  public:
    ///+ 9- Specific to FldArrayI
    ///- Adds val to all values.
    void shift (T val);

  public:
    T& getitems(int i, int j);
    T& getitem(int i);

    void setitem(int i, T value);
    void setitems(int i, int j, T value);
    void setNGonType(int i) {
      if (i < 1 && i > 3) return;
      _ngon = i;
    };

  private:
    // set all elements to val
    void __setAllValuesAt(E_Int sizetotloc, T* to, T val);
    //
    void __setAllValues(E_Int sizetotloc, T* to, const T* from);
    // set a sub array'elements all to val
    void __setSubArrValuesAt(T* data, E_Int first, E_Int last, T val);
    //copy from
    void __cpyFrom(const E_Int& deb, const E_Int& sizerhs, const E_Int& sizelhs,
                   const E_Int& nfld, const T* rhs, T* lhs);

  private:

    /* The allocated size. At initialisation,  _sizeTot = _sizeMax * _nfldMax */
    E_Int _sizeTot;
    /* The used size: _sizeLoc <= _sizeMax */
    E_Int _sizeLoc;
    /* The allocated size: when using the static interface, we always have _sizeMax == _sizeLoc */
    E_Int _sizeMax;
    /* The number of fields: _nfldLoc <= _nfldMax. */
    E_Int _nfldLoc;
    /* The number of allocated fields: when using the static interface, we always have _nfldMax == _nfldLoc */
    E_Int _nfldMax;
    /* data is shared (owned by someone else)? */
    E_Boolean _shared;
    /* is data compact? */
    E_Boolean _compact;
    /* stride */
    E_Int _stride;

    /* si _ngon==1, connect compact NGONv3 array1
       si _ngon==2, connect NGONv3 (stockage rake, NGON, NFACE, indPG, indPH, PE)
       si _ngon==3, connect NGONv4 (stockage rake, NGON, NFACE, indPG, indPH, PE) */
    E_Int _ngon;
    E_Int _nfaces; // nbre de faces (NGON)
    E_Int _nelts; // nbre d'elements (NGON)
    E_Int _sizeNGon; // taille de la connectivite NGon (NGON)
    E_Int _sizeNFace; // taille de la connectivite NFace (NGON)

    /* si multiple element connectivities */
    std::vector< FldArray<T>* > _BEConnects; 

    /* The data array (deprecated) */
    value_type* _data;
    /* The data rake */
    value_type* _rake[512];
};

typedef FldArray<E_Float> FldArrayF;
typedef FldArray<E_Int> FldArrayI;
typedef FldArray<short> FldArrayIS;
typedef FldArray<int> FldArrayI4;
typedef FldArray<E_Boolean> FldArrayB;

//INLINING

#define TEMPLATE_T template<typename T>

//OK============================================================================
TEMPLATE_T
inline T FldArray<T>::operator[] (E_Int l) const
{
  assert (_rake != NULL);
  assert(_ngon == 0);
  assert (l >= 0);
  assert (l < _sizeMax*_nfldLoc);
  E_Int n = E_Int(l/_sizeLoc);
  return _rake[n][(l-n*_sizeLoc)*_stride];
}

//OK============================================================================
TEMPLATE_T
inline T& FldArray<T>::operator[] (E_Int l)
{
  assert (_rake != NULL);
  assert (l >= 0);
  assert (l < _sizeMax*_nfldLoc);
  //return _data[l];
  E_Int n = E_Int(l/_sizeLoc);
  return _rake[n][(l-n*_sizeLoc)*_stride];
}

//OK============================================================================
TEMPLATE_T
inline T FldArray<T>::operator()(E_Int l, E_Int fld) const
{
  assert (_rake  !=  NULL);
  assert (l >= 0);
  assert (l < _sizeLoc);
  assert (fld >= NUMFIELD0);
  assert (fld <= _nfldLoc);
  //return (_data[(fld-NUMFIELD0)*_sizeMax + l]);
  return _rake[fld-NUMFIELD0][l*_stride];
}

//OK============================================================================
TEMPLATE_T
inline T& FldArray<T>::operator()(E_Int l, E_Int fld)
{
  assert (_rake != NULL);
  assert (l >= 0);
  assert (l < _sizeLoc);
  assert (fld >= NUMFIELD0);
  assert (fld <= _nfldLoc);
  //return (_data[(fld-NUMFIELD0)*_sizeMax + l]);
  return _rake[fld-NUMFIELD0][l*_stride];
}

//OK============================================================================
TEMPLATE_T
inline
T* FldArray<T>::begin(E_Int fld)
{
  assert (fld >= NUMFIELD0);
  // either we ask for a existing field in a non-empty array 
  // OR the array is empty and we ask only for the first field i.e ask for &rake[0] 
  assert ( (_nfldLoc && (fld <= _nfldLoc)) || (_nfldLoc == 0 && fld == NUMFIELD0) );
  //return (_data+((fld-NUMFIELD0)*_sizeMax));
  return _rake[fld-NUMFIELD0];
}

//OK============================================================================
TEMPLATE_T
inline
const T* FldArray<T>::begin(E_Int fld) const
{
  assert (fld >= NUMFIELD0);
  assert (fld <= _nfldLoc);
  //return (_data+((fld-NUMFIELD0)*_sizeMax));
  return _rake[fld-NUMFIELD0];
}

//OK============================================================================
TEMPLATE_T
inline
T* FldArray<T>::end()
{
  //return (_data + _sizeMax*(_nfldLoc-1) + _sizeLoc);
  return _rake[_nfldLoc-1]+_sizeLoc*_stride;
}

//OK============================================================================
TEMPLATE_T
inline
const T* FldArray<T>::end() const
{
  //return (_data + _sizeMax*(_nfldLoc-1) + _sizeLoc);
  return _rake[_nfldLoc-1]+_sizeLoc*_stride;
}

//OK============================================================================
TEMPLATE_T
inline
T* FldArray<T>::end(E_Int fld)
{
  assert (fld >= NUMFIELD0);
  assert (fld <= _nfldLoc);
  //return (_data + _sizeMax*(fld-1) + _sizeLoc);
  return _rake[fld-1]+_sizeLoc*_stride;
}

//OK===========================================================================
TEMPLATE_T
inline
const T* FldArray<T>::end(E_Int fld) const
{
  assert (fld >= NUMFIELD0);
  assert (fld <= _nfldLoc);
  //return (_data + _sizeMax*(fld-1) + _sizeLoc);
  return _rake[fld-1]+_sizeLoc*_stride;
}
//OK============================================================================
TEMPLATE_T
inline T& FldArray<T>::getitems(int i, int j)
{
  return this->operator()(i,j);
}

//OK============================================================================
TEMPLATE_T
inline T& FldArray<T>::getitem(int i)
{
  return this->operator[](i);
}

//OK============================================================================
TEMPLATE_T
void FldArray<T>::setitem(int i, T value)
{
  this->operator[](i) = value;
}

//OK============================================================================
TEMPLATE_T
void FldArray<T>::setitems(int i, int j, T value)
{
  this->operator()(i,j) = value;
}

//==============================================================================
TEMPLATE_T
E_Int FldArray<T>::getNFaces()
{
  if (_ngon >= 2) // Array2/3
  {
    return _nfaces;
  }
  else // Array1
  {
    return _rake[0][0];
  }
}

//==============================================================================
TEMPLATE_T
E_Int FldArray<T>::getNElts()
{
  if (_ngon >= 2) // Array2/3
  {
    return _nelts;
  }
  else // Array1
  {
    E_Int sizeNGon = _rake[0][1];
    return _rake[0][sizeNGon+2]; 
  }
  /*
  else // suppose a BE connectivity
  {
      return _sizeLoc;
  }*/
}

//==============================================================================
TEMPLATE_T
E_Int* FldArray<T>::getNGon()
{
  if (_ngon >= 2) // Array2/3
  {
    return _rake[0];
  }
  else // Array1
  {
    return _rake[0]+2;
  }
}

//==============================================================================
TEMPLATE_T
E_Int* FldArray<T>::getNFace()
{
  if (_ngon >= 2) // Array2/3
  {
    return _rake[1];
  }
  else // Array1
  {
    E_Int sizeNGon = _rake[0][1];
    return _rake[0]+sizeNGon+4;
  }
}

//==============================================================================
TEMPLATE_T
E_Int* FldArray<T>::getIndPG()
{
  if (_ngon >= 2) // Array 2/3
  {
    return _rake[2];
  }
  else // Array1
  {
    // if _ngon == 0; -> ngon=1 et cree le rake
    if (_rake[2] == NULL)
    {
      _ngon = 1;
      E_Int nfaces = _rake[0][0];
      E_Int* ptrf = _rake[0]+2;
      E_Int c = 0;
      _rake[2] = new E_Int [nfaces];
      E_Int* indPG = _rake[2];
      for (E_Int i = 0; i < nfaces; i++) { indPG[i] = c; c += ptrf[c]+1; }
    }
    return _rake[2];
  }
}

//==============================================================================
TEMPLATE_T
E_Int* FldArray<T>::getIndPH()
{
  if (_ngon >= 2) // Array 2/3
  {
    return _rake[3];
  }
  else // Array1
  {
    if (_rake[3] == NULL)
    {
      _ngon = 1;
      E_Int sizeNGon = _rake[0][1];
      E_Int nelts = _rake[0][2+sizeNGon];
      E_Int* ptre = _rake[0]+4+sizeNGon;
      E_Int c = 0;
      _rake[3] = new E_Int [nelts];
      E_Int* indPH = _rake[3];
      for (E_Int i = 0; i < nelts; i++) { indPH[i] = c; c += ptre[c]+1; }
    }
    return _rake[3];
  }
}

//==============================================================================
TEMPLATE_T
E_Int* FldArray<T>::getPE()
{
  if (_ngon >= 2) // Array 2/3
  {
    return _rake[3];
  }
  else // Array1
  {
    return NULL; // inexistant en Array1
    // peut etre le cree dans un rake?
  }
}

//==============================================================================
TEMPLATE_T
E_Int FldArray<T>::getSizeNGon()
{
  if (_ngon >= 2) // Array2/3
  {
    return _sizeNGon;
  }
  else // Array1
  {
    return _rake[0][1];
  }
}

//==============================================================================
TEMPLATE_T
E_Int FldArray<T>::getSizeNFace()
{
  if (_ngon >= 2) // Array2/3
  {
    return _sizeNFace;
  }
  else // Array1
  {
    return _data[3+_data[1]];
    //return _rake[0][_rake[0][1]+3];
  } 
}

//==============================================================================
TEMPLATE_T
E_Int* FldArray<T>::getFace(E_Int no, E_Int& size)
{
  if (_ngon == 3) // Array3
  {
    E_Int* ngon = _rake[0];
    E_Int* indPG = _rake[2];
    E_Int pos = indPG[no];
    size = indPG[no+1]-pos;
    return ngon+pos;
  }
  else if (_ngon == 2) // Array 2
  {
    E_Int* ngon = _rake[0];
    E_Int* indPG = _rake[2];
    E_Int pos = indPG[no];
    size = ngon[pos];
    return ngon+pos+1;
  }
  else // Array1
  {
    E_Int* ngon = getNGon();
    E_Int* indPG = getIndPG();
    E_Int pos = indPG[no];
    size = ngon[pos];
    return ngon+pos+1;
  }
}

TEMPLATE_T
E_Int* FldArray<T>::getFace(E_Int no, E_Int& size, E_Int* ngon, E_Int* indPG)
{
  E_Int pos = indPG[no];
  if (_ngon == 3) // Array3
  {
    size = indPG[no+1]-pos;
    return ngon+pos;
  }
  else // Array 1 or 2
  {
    size = ngon[pos];
    return ngon+pos+1;
  }
}

//==============================================================================
TEMPLATE_T
E_Int* FldArray<T>::getElt(E_Int no, E_Int& size)
{
  if (_ngon == 3) // Array3
  {
    E_Int* nface = _rake[1];
    E_Int* indPH = _rake[3];
    E_Int pos = indPH[no];
    size = indPH[no+1]-pos;
    return nface+pos;
  }
  else if (_ngon == 2) // Array 2
  {
    E_Int* nface = _rake[1];
    E_Int* indPH = _rake[3];
    E_Int pos = indPH[no];
    size = nface[pos];
    return nface+pos+1;
  }
  else // Array1
  {
    E_Int* nface = getNFace();
    E_Int* indPH = getIndPH();
    E_Int pos = indPH[no];
    size = nface[pos];
    return nface+pos+1;
  }
}

TEMPLATE_T
E_Int* FldArray<T>::getElt(E_Int no, E_Int& size, E_Int* nface, E_Int* indPH)
{
  E_Int pos = indPH[no];
  if (_ngon == 3) // Array3
  {
    size = indPH[no+1]-pos;
    return nface+pos;
  }
  else // Array 1 or 2
  {
    size = nface[pos];
    return nface+pos+1;
  }
}

//OK---------------------------------------------------------------------------------
TEMPLATE_T
FldArray<T>::FldArray()
  : _sizeTot(0), _sizeLoc(0), _sizeMax(0), _nfldLoc(1), _nfldMax(1),
    _shared(false), _compact(true), _stride(1), _ngon(0), _nfaces(0),
    _nelts(0), _sizeNGon(0),
    _sizeNFace(0), _data(NULL)
{
  _rake[2] = NULL; _rake[3] = NULL;
}

//OK-----------------------------------------------------------------------------
// Constructeur pour des champs (size, nfld) (compact+stride)
TEMPLATE_T
FldArray<T>::FldArray(E_Int size, E_Int nfld, E_Boolean compact,
    E_Boolean fortranOrdered)
  : _sizeTot(size*nfld), _sizeLoc(size), _sizeMax(size),
    _nfldLoc(nfld), _nfldMax(nfld), _shared(false), _compact(compact),
    _stride(1), _ngon(0), _nfaces(0), _nelts(0), 
    _sizeNGon(0), _sizeNFace(0), _data(NULL)
{
  _rake[2] = NULL; _rake[3] = NULL;
  if (fortranOrdered == false) _stride = _nfldLoc;
  acquireMemory();
}

//OK----------------------------------------------------------------------
TEMPLATE_T
FldArray<T>* FldArray<T>::setMemoryCompactZone(value_type* data, E_Int size, E_Int nfld)
{
  return new FldArray<T>(size, nfld, data, true, true);
}

//OK----------------------------------------------------------------------
TEMPLATE_T
FldArray<T>* FldArray<T>::setMemoryNonCompactZone(value_type** data, E_Int size, E_Int nfld)
{
  return new FldArray<T>(size, nfld, data, true, true);
}

//OK---------------------------------------------------------------------------
// Constructeur de copie
TEMPLATE_T
FldArray<T>::FldArray(const FldArray& rhs)
 : _sizeTot(rhs._sizeLoc * rhs._nfldLoc),
   _sizeLoc(rhs._sizeLoc),
   _sizeMax(rhs._sizeLoc),
   _nfldLoc(rhs._nfldLoc),
   _nfldMax(rhs._nfldLoc),
   _shared(false),
   _compact(rhs._compact),
   _stride(rhs._stride),
   _ngon(rhs._ngon),
   _nfaces(rhs._nfaces),
   _nelts(rhs._nelts),
   _data(NULL)
{
  _rake[2] = NULL; _rake[3] = NULL;
  acquireMemory();
  copy(rhs, 1, _nfldLoc);
}

//OK---------------------------------------------------------------------------
TEMPLATE_T
FldArray<T>::FldArray(const FldArray& rhs, E_Int begin, E_Int end)
 : _sizeTot(rhs._sizeLoc * (end-begin+1)),
   _sizeLoc(rhs._sizeLoc),
   _sizeMax(rhs._sizeLoc),
   _nfldLoc(end-begin+1),
   _nfldMax(end-begin+1),
   _shared(false),
   _compact(rhs._compact),
   _stride(rhs._stride),
   _ngon(rhs._ngon),
   _nfaces(rhs._nfaces),
   _nelts(rhs._nelts),
   _data(NULL)
{
  _rake[2] = NULL; _rake[3] = NULL;
  assert (rhs._nfldLoc >= (end-begin+1));
  acquireMemory();
  copy(rhs, begin, end);
}

//OK forcement compact---------------------------------------------------------
TEMPLATE_T
FldArray<T>::FldArray(E_Int size, E_Int nfld,
                      const T* listofvalues, E_Boolean shared,
                      E_Boolean fortranOrdered)
 : _sizeTot(size*nfld),
   _sizeLoc(size),
   _sizeMax(size),
   _nfldLoc(nfld),
   _nfldMax(nfld),
   _shared(shared),
   _compact(true),
   _stride(1),
   _ngon(0),
   _nfaces(0),
   _nelts(0),
   _data(NULL)
{
  _rake[2] = NULL; _rake[3] = NULL;
  if (fortranOrdered == false) _stride = _nfldLoc;

  if (shared == false)
  { // copie
    acquireMemory();
    const T* ptr = listofvalues;

    if (_sizeLoc*_nfldLoc < __MIN_SIZE_LIGHT__)
      memcpy(_data, ptr, _sizeLoc*_nfldLoc*sizeof(T));
    else
      __setAllValues(_sizeLoc*_nfldLoc, _data, ptr);
  }
  else
  { // shared
    //_rake = new T* [nfld];
    _data = (T*)listofvalues; SETRAKE;
  }
}

//OK forcement non compact------------------------------------------------------
TEMPLATE_T
FldArray<T>::FldArray(E_Int size, E_Int nfld,
                      T** listofvalues, E_Boolean shared,
                      E_Boolean fortranOrdered)
 : _sizeTot(size*nfld),
   _sizeLoc(size),
   _sizeMax(size),
   _nfldLoc(nfld),
   _nfldMax(nfld),
   _shared(shared),
   _compact(false),
   _stride(1),
   _ngon(0),
   _nfaces(0),
   _nelts(0),
   _data(NULL)
{
  _rake[2] = NULL; _rake[3] = NULL;
  if (fortranOrdered == false) _stride = _nfldLoc;

  if (shared == false)
  { // copie
    acquireMemory();
    for (E_Int i = 0; i < nfld; i++)
      __setAllValues(_sizeLoc, _rake[i], listofvalues[i]);
  }
  else
  { // shared
    //_rake = new T* [nfld];
    for (E_Int i = 0; i < nfld; i++)
      _rake[i] = listofvalues[i];
  }
}

//OK forcement NGON (array1) ----------------------------------------------------------------------
// Ne pas utiliser: on utilise le constructeur standard dans le cas Array1
TEMPLATE_T
FldArray<T>::FldArray(E_Int nfaces, E_Int nelts, E_Int sizeNGon, E_Int sizeNFace, T* ngonc)
 : _sizeTot(4+sizeNGon+sizeNFace),
   _sizeLoc(4+sizeNGon+sizeNFace),
   _sizeMax(4+sizeNGon+sizeNFace),
   _nfldLoc(1),
   _nfldMax(1),
   _shared(true),
   _compact(true),
   _stride(1),
   _ngon(1),
   _nfaces(nfaces),
   _nelts(nelts),
   _data(NULL)
{
  _data = ngonc;
  _data[0] = nfaces;
  _data[1] = sizeNGon;
  _data[sizeNGon+2] = nelts;
  _data[sizeNGon+3] = sizeNFace;
  _rake[0] = _data; // NGON
  _rake[1] = NULL;
  _rake[2] = NULL; // a ajouter
  _rake[3] = NULL;
}

//OK forcement NGON (array2 ou array3) ---------------------------------------------------------------------
TEMPLATE_T
FldArray<T>::FldArray(E_Int nfaces, E_Int nelts, T* ngon, T* nface, T* indPG, T* indPH, 
                      E_Int sizeNGon, E_Int sizeNFace, T* PE)
 : _sizeTot(nfaces*nelts),
   _sizeLoc(nfaces), // nbre de faces
   _sizeMax(nelts), // nbre d'elements
   _nfldLoc(1),
   _nfldMax(1),
   _shared(true),
   _compact(false),
   _stride(1),
   _ngon(2),
   _nfaces(nfaces),
   _nelts(nelts),
   _data(NULL)
{
  _data = ngon;
  _rake[0] = ngon;
  _rake[1] = nface;
  _rake[2] = indPG;
  _rake[3] = indPH;
  _rake[4] = PE;
  _sizeNGon = sizeNGon;
  _sizeNFace = sizeNFace;
}

//OK===========================================================================
TEMPLATE_T
FldArray<T>::~FldArray()
{
  releaseMemory();
}

//OK===========================================================================
TEMPLATE_T
E_Int FldArray<T>::indMin(E_Int inc, E_Int nfld) const
{
  T val = K_CONST::E_MAX_FLOAT;
  E_Int ind  = 0;
  T* pt = _rake[nfld-1];
  for (E_Int i = 0; i < _sizeLoc; i= i+inc)
  {
    if (pt[i*_stride] < val)
    {
      ind = i;
      val = pt[i*_stride];
    }
  }
  return ind;
}

//OK===========================================================================
TEMPLATE_T
E_Int FldArray<T>::indMax(E_Int inc, E_Int nfld) const
{
  T val = -K_CONST::E_MAX_FLOAT;
  E_Int ind  = 0;
  T* pt = _rake[nfld-1];
  for (E_Int i = 0; i < _sizeLoc; i = i+inc)
  {
    if (pt[i*_stride] > val)
    {
      ind = i;
      val = pt[i*_stride];
    }
  }
  return ind;
}

//----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::sqrt()
{
  #pragma omp parallel
  {
    for (E_Int n = 0; n < _nfldLoc; n++)
    {
        T* pt = _rake[n];
        #pragma omp for
        for (E_Int i = 0; i < _sizeLoc; i++) pt[i*_stride] = std::sqrt(pt[i*_stride]);
    }
  }
}

//OK---------------------------------------------------------------------------
TEMPLATE_T
FldArray<T>& FldArray<T>::operator= (const FldArray<T>& rhs)
{
  if (this != &rhs)
  {
    if (_sizeTot >= (rhs._sizeLoc*rhs._nfldLoc) && rhs._compact == true)
    {
      //if (rhs._nfldMax > _nfldMax) { delete [] _rake; _rake = new T* [rhs._nfldMax]; }
      _sizeMax = _sizeLoc = rhs._sizeLoc;
      _nfldMax = _nfldLoc = rhs._nfldLoc;
      _compact = rhs._compact;
      _stride = rhs._stride;
      _ngon = rhs._ngon;
      SETRAKE;
    }
    else
    {
      releaseMemory();
      _sizeMax = _sizeLoc = rhs._sizeLoc;
      _nfldMax = _nfldLoc = rhs._nfldLoc;
      _sizeTot = _sizeMax * _nfldMax;
      _shared = false;
      _compact = rhs._compact;
      _stride = rhs._stride;
      _ngon = rhs._ngon;
      acquireMemory();
    }
    copy(rhs, 1, _nfldLoc);
  }
  return (*this);
}

SPECIALISE inline void FldArray<E_Float>::__setSubArrValuesAt(E_Float* data, E_Int first, E_Int last, E_Float val) {k6operatoregalf_(data, first, last, val);}
SPECIALISE inline void FldArray<E_Int>::__setSubArrValuesAt(E_Int* data, E_Int first, E_Int last, E_Int val) {for (E_Int i=first; i<last; ++i)data[i] = val;}

//OK---------------------------------------------------------------------------
TEMPLATE_T
FldArray<T>& FldArray<T>::operator= (T val)
{
  if (_compact == true) __setAllValuesAt(_sizeMax*_nfldLoc, _data, val);
  else
  {
    for (E_Int i = 0; i < _nfldLoc; i++) __setAllValuesAt(_sizeLoc, _rake[i], val);
  }
  return (*this);
}

// ----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::copyFrom(E_Int deb, const FldArray<T>& rhs)
{
  assert( rhs.getNfld() == getNfld() );
  __cpyFrom(deb, rhs.getSize(), _sizeMax, _nfldLoc, rhs.begin(), _data);
}

SPECIALISE inline void FldArray<E_Float>::__cpyFrom(
  const E_Int& deb, const E_Int& sizerhs, const E_Int& sizelhs,
  const E_Int& nfld, const E_Float* rhs, E_Float* lhs)
{k6fldcopyfrom_(deb, sizerhs, sizelhs, nfld, rhs, lhs);}
SPECIALISE inline void FldArray<E_Int>::__cpyFrom(
  const E_Int& deb, const E_Int& sizerhs, const E_Int& sizelhs,
  const E_Int& nfld, const E_Int* rhs, E_Int* lhs)
{k6fldintcopyfrom_(deb, sizerhs, sizelhs, nfld, rhs, lhs);}
//SPECIALISE inline void FldArray<short>::__cpyFrom(const E_Int& deb, const E_Int& sizerhs, const E_Int& sizelhs,
//                               const E_Int& nfld, const short* rhs, short* lhs)
//                               {}

//==============================================================================
TEMPLATE_T
void FldArray<T>::print()
{
  printf("------------------------------------------------------------\n");
  printf("_sizeLoc=%d, _sizeMax=%d, _sizeTot=%d\n", _sizeLoc, _sizeMax, _sizeTot);
  printf("_nfldLoc=%d, _nfldMax=%d\n", _nfldLoc, _nfldMax);
  printf("_shared=%d\n", _shared);
  printf("_compact=%d\n", _compact);
  printf("_stride=%d\n", _stride);
  printf("_ngon=%d\n", _ngon);
  printf("_nfaces=%d\n", _nfaces);
  printf("_nelts=%d\n", _nelts);
}

//OK---------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::copy(const FldArray& rhs, E_Int beg, E_Int end)
{
  assert (_sizeMax >= rhs._sizeLoc);
  assert (_nfldMax >= (end-beg+1));
  assert (beg >= 1);
  assert (end <= rhs._nfldLoc);
  // do not thread

  size_t j0(0);
  if (_compact == true)
  {
    for (E_Int j = beg-1; j < end; ++j, ++j0)
    {
      for (E_Int i = 0; i < rhs._sizeLoc; ++i)
        _data[i+j0*_sizeMax] = rhs._data[i+j*rhs._sizeMax];
    }
  }
  else
  {
    for (E_Int j = beg-1; j < end; ++j, ++j0)
    {
      T* pt1 = _rake[j0]; T* pt2 = rhs._rake[j];
      for (E_Int i = 0; i < rhs._sizeLoc; ++i) pt1[i] = pt2[i];
    }
  }
}

//OK---------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::acquireMemory()
{
  if (_compact == true)
  {
    _data = new T[_sizeTot];
    if (_stride == 1)
    {
      for (E_Int i = 0; i < _nfldLoc; i++) _rake[i] = &_data[i*_sizeLoc];
    }
    else
    {
      for (E_Int i = 0; i < _nfldLoc; i++) _rake[i] = &_data[i];
    }
  }
  else // rake
  {
    assert((_nfldLoc > 0) && "_nfldLoc must be greater than zero !");
    assert((_sizeLoc > 0) && "_sizeLoc must be greater than zero !");
    _rake[0] = new T[_sizeLoc];
    for (E_Int i = 1; i < _nfldLoc; i++) _rake[i] = new T[_sizeLoc];
    _data = _rake[0];
  }
}

//OK---------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::releaseMemory(void)
{
  if (_shared == false)
  {
    if (_compact == true) delete [] _data;
    else
    {
      for (E_Int i = 0; i < _nfldMax; i++) delete [] _rake[i];
    }
  }
  _data = NULL;
  for (E_Int i = 0; i < _nfldMax; i++) _rake[i] = NULL;
  _sizeTot = _sizeMax = _sizeLoc = _nfldMax = _nfldLoc = 0;
  if (_ngon == 1) { delete [] _rake[2]; delete [] _rake[3]; }
  //delete [] _rake; _rake = NULL;
  for (size_t i = 0; i < _BEConnects.size(); i++) delete _BEConnects[i];
}

// ---------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::setAllValuesAt(const FldArray<T>& table)
{
  assert (_sizeMax*_nfldLoc <= table.getNfld()*table.getSize());
  if (_sizeMax*_nfldLoc < 100)
    memcpy(_data, table._data, _sizeMax*_nfldLoc*sizeof(T));
  else
    __setAllValues(_sizeMax*_nfldLoc, _data, table._data);
}

TEMPLATE_T
void FldArray<T>::__setAllValues(E_Int size, T* to, const T* from)
{
  for (E_Int i = 0; i < size; i++) to[i] = from[i];
}

#ifndef PURE_C
SPECIALISE inline void FldArray<E_Float>::__setAllValues(E_Int size, E_Float* to, const E_Float* from){k6setallvaluesf_(size, to, from);}
SPECIALISE inline void FldArray<E_Int>::__setAllValues(E_Int size, E_Int* to, const E_Int* from){k6setallvaluesi_(size, to, from);}
#endif

// ---------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::setByField(const FldArray<T>& valuesoffields)
{
  assert(valuesoffields.getSize() == _nfldLoc);
  assert(valuesoffields.getNfld() == 1);

#pragma omp parallel default(shared)
  {
    T vect;
    for (E_Int n = 0; n < _nfldLoc; n++)
    {
      vect = valuesoffields[n];
#pragma omp for
      for (E_Int lsize = _sizeMax*n; lsize < (_sizeMax*n)+_sizeLoc; lsize++)
        _data[lsize] = vect;
    }
  }
}

// ----------------------------------------------------------------------------
// Deprecated: use FldArray<T>::operator=(T) instead
// ----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::setAllValuesAt(T val)
{
  if (_compact == true)
    __setAllValuesAt(_sizeMax*_nfldLoc, _data, val);
  else
  {
    for (E_Int j = 0; j < _nfldLoc; j++) __setAllValuesAt(_sizeLoc, _rake[j], val);
  }
}

TEMPLATE_T
void FldArray<T>::__setAllValuesAt(E_Int size, T* to, T val)
{
//#pragma omp parallel for
  for (E_Int i = 0; i < size; ++i) to[i] = val;
}

#ifndef PURE_C
SPECIALISE inline void FldArray<E_Float>::__setAllValuesAt(E_Int size, E_Float* to, E_Float val) {k6setallvaluesatf_(size, to, val);}
SPECIALISE inline void FldArray<E_Int>::__setAllValuesAt(E_Int size, E_Int* to, E_Int val) {k6setallvaluesati_(size, to, val);}
#endif

// ----------------------------------------------------------------------------
/*
SPECIALISE
inline void FldArray<E_Float>::setAllBorderValuesAt(E_Float val,
                                     E_Int ni, E_Int nj, E_Int nk,
                                     E_Int border)
{
  assert(ni*(nj-1)*(nk-1)+(ni-1)*nj*(nk-1)+(ni-1)*(nj-1)*nk == _sizeLoc);
  k6setallbvaluesatf_(_sizeMax, _nfldLoc, _data, val, ni, nj, nk, border);
}
*/
// ----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::setOneField(const FldArray<T>& fromValArray,
                              E_Int noFieldFrom,
                              E_Int noFieldTo)
{
  assert (noFieldFrom >= NUMFIELD0 &&
          noFieldFrom <= NUMFIELD0+fromValArray.getNfld()-1);
  assert (noFieldTo   >= NUMFIELD0 &&
          noFieldTo   <= NUMFIELD0+_nfldLoc-1);

  E_Int beg = (noFieldTo-1)*_sizeMax;

  const T* ptvalues = fromValArray.begin(noFieldFrom);

#ifdef MEMCOPY
  T* start = _data+beg;
  E_Int size = K_FUNC::E_min(_sizeLoc, fromValArray.getSize());
  memcpy(start, ptvalues, size*sizeof(T));
#else
  E_Int end = noFieldTo*_sizeMax;
  for (E_Int l = beg; l < end; l++) _data[l] = ptvalues[l-beg];
#endif
}

//-----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::resize(E_Int size, E_Int nfld)
{
  if ((size*nfld) <= _sizeTot && _compact == true)
  {
    //if (nfld > _nfldMax) { delete [] _rake; _rake = new T* [nfld]; }
    _sizeMax = _sizeLoc = size;
    _nfldMax = _nfldLoc = nfld;
    SETRAKE;
    // From now on, _sizeLoc*_nfldLoc may be different from _sizeTot!
  }
  else
  {
    malloc(size, nfld);
    _shared = false;
    _compact = true;
  }
}
//-----------------------------------------------------------------------------
SPECIALISE inline E_Float FldArray<E_Float>::badValue() {return K_CONST::E_BADVALUE_F;}
SPECIALISE inline E_Int FldArray<E_Int>::badValue() {return K_CONST::E_BADVALUE_I;}
//SPECIALISE inline E_Boolean FldArray<E_Boolean>::badValue() {return K_CONST::E_BADVALUE_B;}

SPECIALISE inline E_Float FldArray<E_Float>::Zero() {return K_CONST::E_ZERO_FLOAT;}
SPECIALISE inline E_Int FldArray<E_Int>::Zero() {return K_CONST::E_ZERO_INT;}
SPECIALISE inline short FldArray<short>::Zero() {return 0;}

//-----------------------------------------------------------------------------
// Alloue seulement si changement de taille
// On garde les caract. (compact et ordered)
//------------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::malloc(E_Int size, E_Int nfld)
{
  assert(nfld > 0);
  assert(size >= 0);
  releaseMemory();
  _sizeTot = size*nfld;
  _sizeMax = _sizeLoc = size;
  _nfldMax = _nfldLoc = nfld;
  if (_stride > 1) _stride = nfld;
  if (_sizeTot > 0)
  {
    // keep compact
    acquireMemory();
    _shared = false;
    assert (_data != NULL);
  }
  else
  {
    //if (nfld > _nfldMax) { delete [] _rake; _rake = new T* [nfld]; }
    SETRAKE;
  }
}

//-----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::malloc(E_Int size, E_Int nfld, const T* listofvalues,
                         E_Boolean shared)
{
  assert(nfld > 0);
  assert(size >= 0);
  releaseMemory();
  //_rake = new T* [nfld];
  _shared = shared;
  _compact = true;
  _sizeTot = size*nfld;
  _sizeMax = _sizeLoc = size;
  _nfldMax = _nfldLoc = nfld;
  if (_sizeTot == 0) return;
  if (shared == false)
  {
    acquireMemory();
    assert (_data != NULL);
    const T* ptr = listofvalues;

    if (_sizeLoc*_nfldLoc < __MIN_SIZE_LIGHT__)
      memcpy(_data, ptr, _sizeLoc*_nfldLoc*sizeof(T));
    else
      setAllValues(_sizeLoc*_nfldLoc, _data, ptr);
  }
  else
  { // shared
    _data = (T*)listofvalues; SETRAKE;
  }
}

//-----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::reAlloc(E_Int size, E_Int nfld)
{
  T* newdata = NULL;
  E_Int sizeTot = size*nfld;
  assert(sizeTot >= 0);

#ifdef E_DEBUG_MEMORY
  newdata = reinterpret_cast<T*>(::malloc( sizeof(T)*sizeTot));
#else
  newdata = new T[sizeTot];
#endif

  assert(newdata != NULL);

  E_Int minSizeTot = K_FUNC::E_min(sizeTot, _sizeTot);

  if (minSizeTot < __MIN_SIZE_LIGHT__)
    memcpy(newdata, _data, minSizeTot*sizeof(T));
  else
    __setAllValues(minSizeTot, newdata, _data);

  releaseMemory();
  _shared = false;
  _compact = true;
  _sizeTot = size*nfld;
  _sizeMax = _sizeLoc = size;
  _nfldMax = _nfldLoc = nfld;
  _data = newdata;
  //_rake = new T* [nfld];
  SETRAKE;
}

//-----------------------------------------------------------------------------
// For a matrix, copy and keep elements at the same position in the matrix
// for ex: A(12,12) will be B(12,12) after reAlloc
// This is not the case with the previous method.
//-----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::reAllocMat(E_Int size, E_Int nfld)
{
  T* newdata = NULL;
  E_Int sizeTot = size*nfld;
  assert(sizeTot >= 0);

#ifdef E_DEBUG_MEMORY
    newdata = (T*) ::malloc( sizeof(T)*sizeTot);
#else
    newdata = new T[sizeTot];
#endif

  assert(newdata != NULL);

  E_Int sz = K_FUNC::E_min(size, _sizeLoc);
  E_Int nf = K_FUNC::E_min(nfld, _nfldLoc);

#pragma omp parallel default(shared) if (sizeTot > __MIN_SIZE_MEAN__)
  for (E_Int j = 0; j < nf; j++)
  {
#pragma omp for
    for (E_Int i = 0; i < sz; i++)
      newdata[i+j*size] = _data[i+j*_sizeLoc];
#pragma omp for
    for (E_Int i = sz; i < size; i++)
      newdata[i+j*size] = 0;
  }
#pragma omp for
  for (E_Int i = size*nf; i < sizeTot; i++) newdata[i] = 0;

  releaseMemory();
  _shared = false;
  _compact = true;
  _sizeTot = size*nfld;
  _sizeMax = _sizeLoc = size;
  _nfldMax = _nfldLoc = nfld;
  _data = newdata;
  //_rake = new T* [nfld];
  SETRAKE;
}
//-----------------------------------------------------------------------------
// For a matrix, copy and keep elements at the same position in the matrix
// for ex: A(12,12) will be B(12,12) after reAlloc
// This is not the case with the previous method.
// No open MP.
//-----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::reAllocMatSeq(E_Int size, E_Int nfld)
{
  T* newdata = NULL;
  E_Int sizeTot = size*nfld;
  assert(sizeTot >= 0);

#ifdef E_DEBUG_MEMORY
    newdata = (T*) ::malloc( sizeof(T)*sizeTot);
#else
    newdata = new T[sizeTot];
#endif

  assert(newdata != NULL);

  E_Int sz = K_FUNC::E_min(size, _sizeLoc);
  E_Int nf = K_FUNC::E_min(nfld, _nfldLoc);

  for (E_Int j = 0; j < nf; j++)
  {
    for (E_Int i = 0; i < sz; i++)
      newdata[i+j*size] = _data[i+j*_sizeLoc];
    for (E_Int i = sz; i < size; i++)
      newdata[i+j*size] = 0;
  }
  for (E_Int i = size*nf; i < sizeTot; i++) newdata[i] = 0;

  releaseMemory();
  _shared = false;
  _compact = true;
  _sizeTot = size*nfld;
  _sizeMax = _sizeLoc = size;
  _nfldMax = _nfldLoc = nfld;
  _data = newdata;
  //_rake = new T* [nfld];
  SETRAKE;
}
//-----------------------------------------------------------------------------
/*
SPECIALISE
inline void FldArray<E_Float>::matMultVec(const FldArray<E_Float>& array)
{
  const E_Float* tab = array.begin();
  k6multmatvecf_(_sizeLoc, _nfldLoc, _data, tab);
}
*/
//-----------------------------------------------------------------------------
/*
SPECIALISE
inline void FldArray<E_Float>::matMult(FldArray<E_Float>& array)
{
#ifdef DEBUG
  E_Int sizeOfarray = array._sizeLoc;
  E_Int nfldOfarray = array._nfldLoc;
  assert(sizeOfarray == nfldOfarray);
  assert(_sizeLoc   == _nfldLoc);
  assert(_sizeLoc   == sizeOfarray);
#endif

  E_Int il, ic, lc;
  E_Float* dataProv = new value_type[_sizeLoc*_nfldLoc];
  E_Float resProv;
  for (E_Int l = 0; l < _sizeLoc; l++)
  for (E_Int c = 1; c <= _nfldLoc; c++)
  {
    resProv = 0.;
    lc = c-1 + l*_nfldLoc;
    for (E_Int i = 1; i <= _sizeLoc; i++)
    {
      ic = (i-1)+l*_sizeLoc;
      il = c-1 + (i-1)*_sizeLoc;
      resProv = resProv + _data[ic]*array._data[il];
    }
    dataProv[lc] = resProv;
  }

  for (E_Int l = 0; l < _sizeLoc*_nfldLoc; l++)
    _data[l] = dataProv[l];

  delete [] dataProv;
}
*/
//-----------------------------------------------------------------------------
/*
SPECIALISE
inline void FldArray<E_Float>::matMultTrans(FldArray<E_Float>& array)
{
#ifdef DEBUG
  E_Int sizeOfarray = array._sizeLoc;
  E_Int nfldOfarray = array._nfldLoc;
  assert (sizeOfarray == nfldOfarray);
  assert (_sizeLoc   == _nfldLoc);
  assert (_sizeLoc   == sizeOfarray);
#endif

  E_Int il, ic, lc;
  E_Float* dataProv = new value_type[_sizeLoc*_nfldLoc];
  E_Float resProv;
  for (E_Int l = 0; l < _sizeLoc; l++)
  for (E_Int c = 1; c <= _nfldLoc; c++)
  {
    resProv = 0.;
    lc = c-1 + l*_nfldLoc;
    for (E_Int i = 1; i <= _sizeLoc; i++)
    {
      ic = (i-1)+l*_sizeLoc;
      il = i-1 + (c-1)*_sizeLoc;
      resProv = resProv + _data[ic]*array._data[il];
    }
    dataProv[lc] = resProv;
  }

  for (E_Int l = 0; l < _sizeLoc*_nfldLoc; l++)
    _data[l] = dataProv[l];

  delete [] dataProv;
}
*/
//OK===========================================================================
TEMPLATE_T
void FldArray<T>::setAllValuesAtNull()
{
  if (_compact == true)
  {
    if (_sizeMax*_nfldLoc < __MIN_SIZE_LIGHT__)
      memset(_data, 0, _sizeMax*_nfldLoc*sizeof(T));
    else
      __setAllValuesAt(_sizeMax*_nfldLoc, _data, Zero());
  }
  else
  {
    for (E_Int i = 0; i < _nfldLoc; i++)
      __setAllValuesAt(_sizeLoc, _rake[i], Zero());
  }
}

//------------------------------------------------------------------------------
// DYNAMIC INTERFACE (Il faut shared=false, compact=true, stride=1)
//------------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::reserve(E_Int nfld, E_Int size)
{
  if ((size*nfld) <= _sizeTot) return;

  T* newdata = NULL;
  E_Int sizeTot = size*nfld;
  assert(sizeTot >= 0);

  newdata = new T[sizeTot];

  assert(newdata != NULL);

  E_Int sz = K_FUNC::E_min(size, _sizeLoc);
  E_Int nf = K_FUNC::E_min(nfld, _nfldLoc);

  for (E_Int j = 0; j < nf; j++)
  {
    for (E_Int i = 0; i < sz; i++)
      newdata[i+j*size] = _data[i+j*_sizeMax];
  }

  releaseMemory();
  _shared = false;
  _compact = true;
  _sizeTot = sizeTot;
  _sizeLoc = sz;
  _sizeMax = size;
  _nfldLoc = nf;
  _nfldMax = nfld;
  _data    = newdata;
  //_rake = new T* [nfld];
  SETRAKE;
}

//-----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::resize(E_Int nfld, E_Int size, const T& val)
{
  // WARNING ; structure (_nfldMax, _sizeMax) are changed only when reallocation is required.

  nfld = std::max(E_Int(1), nfld);

  E_Int sizeTot = size*nfld;
  E_Int sz = K_FUNC::E_min(size, _sizeLoc);
  E_Int nf = K_FUNC::E_min(nfld, _nfldLoc);

  bool realloc = (nfld > _nfldMax || size > _sizeMax);

  T* data = _data;
  if (realloc)
  {
    data = new T[sizeTot];
    assert(data != NULL);
    for (E_Int j = 0; j < nf; j++)
    {
      for (E_Int i = 0; i < sz; i++) // pass the relevant exisiting data
        data[i+j*size] = _data[i+j*_sizeMax];
    }
  }

  for (E_Int j = 0; j < nf; j++)
  {
    for (E_Int i = sz; i < size; i++) // Fills the stride's tail
      data[i+j*size] = val;
  }

  for (E_Int j = nf; j < nfld; j++) // Fills the array's tail
    for (E_Int i = 0; i < size; i++)
      data[i+j*size] = val;

  if (realloc)
  {
    releaseMemory();
    _data = data;
    _shared = false;
    _compact = true;
    _sizeTot = sizeTot;
    _sizeMax = size;
    _nfldMax = nfld;
    //_rake = new T* [nfld];
  }
  else
  {
    //if (nfld > _nfldMax) { delete [] _rake; _rake = new T* [nfld]; }
  }

  _sizeLoc = size;
  _nfldLoc = nfld;
  SETRAKE;
}

// ----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::pushBack(const FldArray<T>& a)
{
  if (_sizeLoc == 0)// Empty array
  {
    *this = a;
    return;
  }
  else if (_nfldLoc != a._nfldLoc) return;

  if (a._sizeLoc> (_sizeMax - _sizeLoc))     // If the spare room is too tight.
    reserve(_nfldLoc, 2*(_sizeLoc+a._sizeLoc));

  for (E_Int j = 0; j < _nfldLoc; j++)
  {
    for (E_Int i = 0; i < a._sizeLoc; i++)
      _data[i + (j*_sizeMax) + _sizeLoc] = a._data[i + j*a._sizeMax];
  }

  _sizeLoc += a._sizeLoc;
}

// ----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::pushBack(const FldArray<T>& a, E_Int i)
{
  if ((_nfldLoc != a._nfldLoc) && (_sizeLoc != 0))
    return;

  if (_sizeMax == _sizeLoc) {    // If the spare room is too tight.
    reserve(a._nfldLoc, 2*(_sizeLoc+1));
    if (_sizeLoc == 0)_nfldLoc=a._nfldLoc;
  }

  for (E_Int j = 0; j < _nfldLoc; j++)
    _data[(j*_sizeMax) + _sizeLoc] = a._data[i + j*a._sizeMax];

  ++_sizeLoc;
}

// ----------------------------------------------------------------------------
TEMPLATE_T
template <typename Iterator>
void
FldArray<T>::pushBack(Iterator begin, Iterator end)
{
  E_Int nf = end-begin;
  if (_sizeLoc == 0)// Empty array
    _nfldLoc = nf;
  else if (_nfldLoc != nf)
    return;

  if (_sizeLoc == _sizeMax) // If the spare room is too tight.
    reserve(_nfldLoc, 2*(_sizeLoc+1));    // Double the capacity.

  for (E_Int j = 0; j < _nfldLoc; j++)
    _data[j*_sizeMax + _sizeLoc] = *(begin++);

  ++_sizeLoc;
}

// ----------------------------------------------------------------------------
template <typename T>
void FldArray<T>::pushBack(const std::vector<T>& a)
{
  if (_nfldLoc > 1)  return;
  _nfldLoc = 1;

  E_Int sz = a.size();

  if (sz > (_nfldMax*_sizeMax - _sizeLoc)) // If the spare room is too tight.
    reserve(1, 2*(_sizeLoc+sz));

  for (E_Int i = 0; i < sz; i++)
    _data[i+_sizeLoc] = a[i];

  _sizeLoc += a.size();
}

// ----------------------------------------------------------------------------
TEMPLATE_T
void FldArray<T>::shift(T val)
{
  for (E_Int i = 0; i < _sizeLoc*_nfldLoc; ++i)
    _data[i] += val;
}

//------------------------------------------------------------------------------
TEMPLATE_T
template <typename Vector1, typename Vector2>
E_Int FldArray<T>::compact(FldArray& a, const Vector1& keep, Vector2& new_Ids)
{
  E_Int   cols(a._sizeLoc);
  // Fast returns
  if (cols == 0)   return 0;

  assert (keep.size() == (size_t)cols);

  bool    carry_on(false);
  E_Int    i1(0), i2(cols-1);

  new_Ids.clear();
  new_Ids.resize(cols, E_IDX_NONE);

  do {

    while ((i1 <= i2) && keep[i1]){new_Ids[i1] = i1; ++i1;}  // Get the first empty column.
    while ((i1 <= i2) && !keep[i2]){--i2;} // Get the first column to move from the tail.

    carry_on = (i1 < i2);

    if (carry_on)
    { // Copy column i2 in column i1.
      new_Ids[i2] = i1;

      for (E_Int j = 0; j < a._nfldLoc; j++)
        a._data[j*a._sizeMax + i1] =  a._data[j*a._sizeMax + i2];

      i2--; i1++;
    }
  }
  while (carry_on);

  a.resize (a._nfldLoc, i1, 0);

  return (cols - i1);
}

// ME connectivities
TEMPLATE_T
FldArray<T>::FldArray(std::vector< FldArray<T>* >& list)
{
  _ngon = 0;
  _compact = false;
  _nfldMax = 0;
  _BEConnects = list;
}

TEMPLATE_T
FldArray<T>* FldArray<T>::getConnect(E_Int i)
{
  if (_BEConnects.size() == 0) return this;
  else return _BEConnects[i];
}

TEMPLATE_T
size_t FldArray<T>::getNConnect()
{
  E_Int size = (E_Int)_BEConnects.size();
  if (size == 0) return 1;
  else return size;
}
} // End namespace K_FLD
#endif

