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

#include "kcore.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

// Les differents tests
# define EXTARITH       0
# define TESTMEMORY     0
# define TESTLOGGER     0
#define  TEST1          0
#define  TEST2          0
#define  TEST3          0
#define  TESTARRAY2     0
#define  TESTFLD        0
#define  TESTNUMPY      0
#define  TESTHIGHORDER  2

#if EXTARITH == 1
#include "ExtArith/quad_double.hpp"
#endif
#if TESTMEMORY == 1
#include "Memory/shared_ptr.hpp"
#include "Memory/vector_view.hpp"
#include "Memory/unique_ptr.hpp"
namespace {
  class TestMemory {
  public:
    static int nbObjs;
    TestMemory( int v ) : val(v)  { nbObjs += 1; }
    ~TestMemory() { nbObjs -= 1; }
    int val;
  private:
    TestMemory( const TestMemory& );
  };
  int TestMemory::nbObjs = 0;
}
# if __cplusplus < 201103L
#   define nullptr NULL
# endif
#endif
#if TESTLOGGER == 1
# include "Logger/logger.hpp"
# include "Logger/log_to_file.hpp"
# include "Logger/log_from_root_output.hpp"
# include "Logger/log_to_std_output.hpp"
# include "Logger/log_from_distributed_file.hpp"
#endif
//==============================================================================
PyObject* K_KCORE::tester(PyObject* self, PyObject* args)
{

#if TESTNUMPY == 1
  K_FLD::FldArrayF f(1,3);
  for (E_Int i = 0; i < 1; i++) f(i,1) = 1.;
    for (E_Int i = 0; i < 1; i++) f(i,2) = 2.;
      for (E_Int i = 0; i < 1; i++) f(i,3) = 3.;

  //f.setAllValuesAt(2.);
  PyObject* a = K_NUMPY::buildNumpyArray(f, 1);
  K_FLD::FldArrayF* out;
  K_NUMPY::getFromNumpyArray(a , out, true);
  return a;
#endif

#if TEST1 == 1
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  E_Float* t = new E_Float[1]; t[0] = 12.;
  PyObject* p = K_PYTREE::createChild(o, "coucou", "DataArray_t", 
                                      t, 1, 1);
  delete [] t;
  return p;
#endif

#if TEST2 == 1
  PyObject* o; char* path;
  if (!PyArg_ParseTuple(args, "Os", &o, &path)) return NULL;
  PyObject* p = K_PYTREE::getNodeFromPath(o, path);
  return p;
#endif

#if TEST3 == 1
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  f->print();
  E_Float* x = f->begin();
  x[0] = -0.05;
  if (ret == 2) c->print();
  RELEASESHAREDB(ret, o, f, c);
  Py_INCREF(Py_None);
  return Py_None;
#endif

#if TESTARRAY2 == 1
  // Structured array1 - build - ni=5,nj=5,nk=5
  PyObject* o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 5,5,5, 1);
  // Structured array1 - get
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType; E_Int ni,nj,nk;
  E_Int ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // Structured array1 - usage
  f->setAllValuesAt(1);
  // Getting information from Fld
  E_Int nfld = f->getNfld(); // nbre de champs
  E_Int size = f->getSize(); // nbre de pts (ni x nj x nk)
  E_Float* x = f->begin(1); // ptrs sur les champs
  E_Float* y = f->begin(2);

  // Structured array1 - free
  RELEASESHAREDB(ret, o, f, c);

  // Structured array2 - build - ni=2,nj=2,nk=2
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 2,2,2, 2);
  // Classic array2 - get
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // Structured array2 - usage
  f->setAllValuesAt(2.2);
  // Getting information from Fld f
  nfld = f->getNfld(); // nbre de champs
  size = f->getSize(); // nbre de pts
  x = f->begin(1);
  y = f->begin(2);

  // Structured array2 - free
  RELEASESHAREDB(ret, o, f, c);

  // Hexa array1 - build - nvertex=10, nelts=5
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 10, 5,-1, "HEXA", false, 0, 0, 0, 1);
  // Hexa array1 - get
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // Hexa array1 - usage
  f->setAllValuesAt(1);
  c->setAllValuesAt(1);
  // Getting information from Fld c
  E_Int ne = c->getSize(); // nbre d'elements
  E_Int nv = c->getNfld(); // nbre de vertex par element
  // Hexa array1 - free
  RELEASESHAREDB(ret, o, f, c);

  // Hexa array2 - build - nvertex=10, nelts=5
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 10, 5,-1, "HEXA", false, 0, 0, 0, 2);
  // Hexa array2 - get
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // Hexa array1 - usage
  f->setAllValuesAt(1);
  c->setAllValuesAt(1);
  // Getting information from Fld c
  ne = c->getSize(); // nbre d'elements
  nv = c->getNfld(); // nbre de vertex par element
  // Hexa array2 - free
  RELEASESHAREDB(ret, o, f, c);

  // NGon array1 - build - nvertex=10, nelts=5, nface=6, sizeNGon=20, sizeNFace=30
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 10, 5,-1, "NGON", false, 20, 30, 6, 1);
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // NGon array1 - usage
  f->setAllValuesAt(1);
  c->setAllValuesAt(1);
  E_Int nelts = c->getNElts();
  E_Int nfaces = c->getNFaces();
  E_Int* ngon = c->getNGon();
  E_Int* nface = c->getNFace();
  E_Int* indPG = c->getIndPG();
  E_Int* indPH = c->getIndPH();
  RELEASESHAREDB(ret, o, f, c);
  //printf("nelts=%d, nfaces=%d\n", nelts, nfaces);

  // NGon array2 - build
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 10, 5,-1, "NGON", false, 20, 30, 6, 2);
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // NGon array1 - usage
  f->setAllValuesAt(1);
  nelts = c->getNElts();
  nfaces = c->getNFaces();
  ngon = c->getNGon();
  nface = c->getNFace();
  indPG = c->getIndPG();
  indPH = c->getIndPH();
  //printf("nelts=%d, nfaces=%d\n", nelts, nfaces);
  RELEASESHAREDB(ret, o, f, c);
  return o;
#endif

#if TESTHIGHORDER == 1
  char eltType[128];
  E_Int loc, nvpe, typeId;
  K_ARRAY::eltString2TypeId((char*)"TETRA_20*", eltType, nvpe, loc, typeId);
  printf("%s %d %d\n", eltType, nvpe, loc);
  K_ARRAY::eltString2TypeId((char*)"TETRA", eltType, nvpe, loc, typeId);
  printf("%s %d %d\n", eltType, nvpe, loc);
  
  K_ARRAY::eltString2TypeId((char*)"NODE", eltType, nvpe, loc, typeId);
  printf("%s %d %d\n", eltType, nvpe, loc);
  
  K_ARRAY::eltString2TypeId((char*)"QUAD_8", eltType, nvpe, loc, typeId);
  printf("%s %d %d\n", eltType, nvpe, loc);
  
  K_ARRAY::eltString2TypeId((char*)"QUAD", eltType, nvpe, loc, typeId);
  printf("%s %d %d\n", eltType, nvpe, loc);
    
  K_ARRAY::eltString2TypeId((char*)"BAR_3", eltType, nvpe, loc, typeId);
  printf("%s %d %d\n", eltType, nvpe, loc);
  K_ARRAY::eltString2TypeId((char*)"TRI_6", eltType, nvpe, loc, typeId);
  printf("%s %d %d\n", eltType, nvpe, loc);
  
  K_ARRAY::eltString2TypeId((char*)"NGON", eltType, nvpe, loc, typeId);
  printf("%s %d %d\n", eltType, nvpe, loc);
#endif

#if TESTHIGHORDER == 2
  printf("test 2\n");
  // Construction vide
  PyObject* o = K_ARRAY::buildArray2(3, "x,y,z", 6, 1, -1, "TRI_6", false, 0, 0, 0, 2);
  // Recuperation des pointeurs
  //E_Int ni,nj,nk; char* eltType; K_FLD::FldArrayF* f; K_FLD::FldArrayI* c; char* varString; 
  //E_Int ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);

  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c; 
  E_Int ret = K_ARRAY::getFromArray2(o, f, c);
  
  // Champs --
  // Acces direct Fld
  (*f)(0,1) = 0.; (*f)(0,2) = -2.; (*f)(0,3) = 0.;
  (*f)(1,1) = 2.; (*f)(1,2) =  0.; (*f)(1,3) = 0.;
  (*f)(2,1) = -2.; (*f)(2,2) =  0.; (*f)(2,3) = 0.;
  (*f)(3,1) =  1.; (*f)(3,2) =  -1.; (*f)(3,3) = 0.;
  (*f)(4,1) =  0.; (*f)(4,2) =   0.; (*f)(4,3) = 0.2;
  (*f)(5,1) = -1.; (*f)(5,2) =  -1.2; (*f)(5,3) = 0.;

  (*c)(0,1) = 1; (*c)(0,2) = 2; (*c)(0,3) = 3; (*c)(0,4) = 4; (*c)(0,5) = 5; (*c)(0,6) = 6;

  // Acces par champ Fld
  //E_Float* f3 = f->begin(3);
  //f3[0] = 0.; f3[1] = 1.;

  // Connectivite
  // Acces direct 
  //(*c)(0,1) = 1; (*c)(0,2) = 2; (*c)(0,3) = 3; (*c)(0,4) = 4; (*c)(0,5) = 5; (*c)(0,6) = 6;
  // Parcours avec stride
  //E_Int* cp = c->begin();
  //E_Int nelts = c->getSize(); E_Int nvpe = c->getNfld(); E_Int s = c->getStride();
  //for (E_Int i = 0; i < nelts; i++) { cp[0] = i; cp += s; }


  // Acces direct pour decoupage de la boucle
  
  RELEASESHAREDB(ret, o, f, c);
  return o;
#endif

#if TESTFLD == 1
  // test nouveau FldArray
  K_FLD::FldArrayF t1; 

  //=========================
  // compact, fortranOrdered
  //========================
  K_FLD::FldArrayF t2(120, 3, true);
  t2 = 0.;
  t2.setAllValuesAt(1);
  t2.setAllValuesAtNull();
  for (E_Int i = 0; i < 120; i++)
  {
    t2(i,1) = i;
    t2(i,2) = 0.1*i;
    t2(i,3) = 0.01*i;
  }
  K_FLD::FldArrayF t2p(t2);
  printf("indmin %d %d %d\n", t2p.indMin(1,1), t2p.indMin(1,2), t2p.indMin(1,3));
  printf("indmax %d %d %d\n", t2p.indMax(1,1), t2p.indMax(1,2), t2p.indMax(1,3));
  printf("1: %f %f %f\n", t2(10,1), t2(10,2), t2(10,3));
  printf("1: %f %f %f\n", t2[10], t2[10+120], t2[10+240]);
  E_Float* pt1 = t2.begin(1);
  E_Float* pt2 = t2.begin(2);
  E_Float* pt3 = t2.begin(3);
  printf("1: %f %f %f\n", pt1[10], pt2[10], pt3[10]);
  
  K_FLD::FldArrayF t2t(120, 3, t2.begin(), true); // shared with t2
  t2t(10,2) = 12.;

  // Build shared
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z", 10, 10, 10);
  E_Float* sp = K_ARRAY::getFieldPtr(tpl);
  K_FLD::FldArrayF s(10*10*10, 3, sp, true);
  
  for (E_Int n = 1; n <= 3; n++)
  {
    E_Float* spi = s.begin(n);
    for (E_Int i = 0; i < 10*10*10; i++) spi[i] = 12.;
  }

  //============================
  // non compact, fortranOrdered
  //===========================
  K_FLD::FldArrayF t3(120, 3, false);
  t3 = 0.;
  t3.setAllValuesAt(1);
  t3.setAllValuesAtNull();
  for (E_Int i = 0; i < 120; i++)
  {
    t3(i,1) = i;
    t3(i,2) = 0.1*i;
    t3(i,3) = 0.01*i;
  }
  K_FLD::FldArrayF t3p(t3);
  printf("indmin %d %d %d\n", t3p.indMin(1,1), t3p.indMin(1,2), t3p.indMin(1,3));
  printf("indmax %d %d %d\n", t3p.indMax(1,1), t3p.indMax(1,2), t3p.indMax(1,3));
  printf("2: %f %f %f\n", t3(10,1), t3(10,2), t3(10,3));
  printf("2: %f %f %f\n", t3[10], t3[10+120], t3[10+240]);
  pt1 = t3.begin(1);
  pt2 = t3.begin(2);
  pt3 = t3.begin(3);
  printf("2: %f %f %f\n", pt1[10], pt2[10], pt3[10]);

  //===================
  // compact, C ordered
  //===================
  K_FLD::FldArrayF t4(120, 3, true, false);
  t4 = 0.;
  t4.setAllValuesAt(1);
  t4.setAllValuesAtNull();
  for (E_Int i = 0; i < 120; i++)
  {
    t4(i,1) = i;
    t4(i,2) = 0.1*i;
    t4(i,3) = 0.01*i;
  }
  K_FLD::FldArrayF t4p(t4);
  printf("indmin %d %d %d\n", t4p.indMin(1,1), t4p.indMin(1,2), t4p.indMin(1,3));
  printf("indmax %d %d %d\n", t4p.indMax(1,1), t4p.indMax(1,2), t4p.indMax(1,3));
  printf("3: %f %f %f\n", t4(10,1), t4(10,2), t4(10,3));
  printf("3: %f %f %f\n", t4[10], t4[10+120], t4[10+240]);
  pt1 = t4.begin(1);
  pt2 = t4.begin(2);
  pt3 = t4.begin(3);
  printf("3: %f %f %f\n", pt1[10*3], pt2[10*3], pt3[10*3]);
#endif

#if TESTMEMORY == 1
  // Unique pointer tests
  // ================================================================================
  {
    printf("Testing unique_ptr implementation\n");
    K_MEMORY::unique_ptr<TestMemory> uptr(new TestMemory(314));
    assert((TestMemory::nbObjs == 1) && "Wrong number of TestMemory objects ! Would be one at creation of unique_ptr!");
#   if __cplusplus > 199711L  
    K_MEMORY::unique_ptr<TestMemory> uptr2 = std::move(uptr);
#   else
    K_MEMORY::unique_ptr<TestMemory> uptr2 = uptr;
#   endif
    assert((TestMemory::nbObjs == 1) && "Wrong number of TestMemory objects ! Would be one after move/copy of unique_ptr!");
    assert((uptr.get() == nullptr) && "uptr pointer would be null after moving operation !" );
    assert((!uptr) && "For null pointer, boolean conversion would be return false and not true !");
    assert((uptr2)&&"For not null pointer, boolean conversion would be return true and not false !");
    assert(((*uptr2).val == 314) && "Wrong value for the attribute of the object after deferencing !");
    assert((uptr2->val == 314) && "Wrong value for the attribute of the object after acccesing with -> operator !");
    uptr2.reset(new TestMemory(315));
    uptr.reset (new TestMemory(313));
    assert((TestMemory::nbObjs == 2) && "Wrong number of TestMemory objects ! Would be two after resetting two unique_ptr!");
    assert((uptr2->val == 315) && "Wrong value for the attribute of the object after resetting the pointer !");
    uptr2.swap(uptr);
    assert((uptr2->val == 313) && "Wrong value for the attribute of the object after swapping two unique_ptr !");
    K_MEMORY::unique_ptr<TestMemory>::pointer pt = uptr.release();
    assert((pt != nullptr) && "The released pointer would be not null in the test !");
    delete pt;
    assert((TestMemory::nbObjs == 1) && "Wrong number of TestMemory objects ! Would be one after releasing one unique_ptr!");
  }
  assert((TestMemory::nbObjs == 0) && "Memory leaks ! The number of TestMemory object would be zero after unique_ptr test !");
  // shared_ptr tests
  // ==================================================================================
  {
    printf("Testing shared_ptr implementation\n");
    K_MEMORY::shared_ptr<TestMemory> sptr(new TestMemory(314));
    assert((TestMemory::nbObjs == 1) && "Wrong number of TestMemory objects ! Would be one at creation of shared_ptr!");
    assert((sptr.use_count() == 1) && "Wrong number of reference for shared_ptr : at creation, would be only one !");
    K_MEMORY::shared_ptr<TestMemory> sptr2(sptr);
    assert((TestMemory::nbObjs == 1) && "Wrong number of TestMemory objects ! Would be still one after constructor copy of a shared_ptr!");
    assert((sptr.use_count() == 2) && "Wrong number of reference for shared_ptr : after constructor copy, would be two !");
    sptr2.reset();
    assert((sptr.use_count() == 1) && "Wrong number of reference for shared_ptr : after releasing on pointer, would be only one !");
    assert((sptr) && "Would be return true and not false for non null pointer");
    assert((!sptr2) && "Would be return false and not true for null pointer");
    TestMemory* pt = sptr.get();
    assert((pt != NULL) && "Non null shared_ptr would return a non null pointer !");
    pt = sptr2.get();
    assert((pt == NULL) && "Null shared_ptr would return a null pointer !");    
    sptr2 = sptr;
    assert((sptr.use_count() == 2) && "Wrong number of reference for shared_ptr : after operator copy, would be two !");
    assert((TestMemory::nbObjs == 1) && "Wrong number of TestMemory objects ! Would be still one after operator copy of a shared_ptr!");
    assert((sptr->val == 314) && "Wrong pointed value for shared_ptr !");
    sptr2.reset(new TestMemory(315));
    assert((TestMemory::nbObjs == 2) && "Wrong number of TestMemory objects ! Would be two after reseting with new value of a shared_ptr!");
    assert((sptr2->val == 315) && "Wrong pointed value for shared_ptr after resetting for new pointer !");
    assert(((*sptr2).val == 315) && "Wrong value when using deferencing operator for shared_ptr !");
    sptr.swap(sptr2);
    assert((sptr->val == 315) && "Wrong pointed value for shared_ptr after swapping with another shared_ptr !");
    assert((sptr2->val == 314) && "Wrong pointed value for shared_ptr after swapping with previous shared_ptr !");
    std::cout << "address of pointer : " << sptr << std::endl;
    sptr2.reset();
    assert((TestMemory::nbObjs == 1) && "Wrong number of TestMemory objects ! Would be two after reseting a shared_ptr with one reference count!");
  }
  assert((TestMemory::nbObjs == 0) && "Memory leaks ! The number of TestMemory object would be zero after shared_ptr test !");
  // vector_view tests
  // ==================================================================================
  {
    double tab[2][7] = { {1.,2.,3.,4.,5.,6.,7.},
                         {7.,6.,5.,4.,3.,2.,1.} };
    printf("Testing vector_view implementation\n");
    K_MEMORY::vector_view<double> view_u;
    assert((view_u.data() == nullptr) && "Default constructor of vector_view would return a null pointer !");
    assert((view_u.size() == 0) && "Default constructor of vector_view would return a zero sized vector !" );
    assert((view_u.empty()) && "Default constructor would be return an empty vector !");
    K_MEMORY::vector_view<double> view_v(&tab[0][0], 7);
    assert((view_v.data() == &tab[0][0]) && "The view memory of view_v would be the same as the static array !");
    assert((view_v.size() == 7) && "The size of view_v would be seven !");
    assert((!view_v.empty()) && "view_v would be not empty !");
#   ifndef NDEBUG
    for ( int i = 0; i < view_v.size(); ++i ) assert((view_v[i] == i+1.) && "Wrong value. Index problem in operator [] ?" );
#   endif
    view_v[3] = 8;
    assert((tab[0][3] == 8.) && "Wrong modification of the array through the view_vector view_v !");
    std::vector<double> w(10);
    for (size_t i = 0; i < w.size(); ++i ) w[i] = (3*i+1);
    K_MEMORY::vector_view<double> view_w(w);
    assert((view_w.data() == &w[0]) && "view_w build from w and has not same memory adress for the first data ?");
    assert((view_w.size() == w.size()) && "view_w build from w has not the same size than w ?" );
    K_MEMORY::vector_view<double> view_x(&tab[1][1],&tab[1][6]);
    assert((view_x.data() == &tab[1][1]) && "Wrong memory view adress in static table with partiel viewing");
    assert((view_x.size() == 5) && "Wrong size of view_x with partial viewing");
#   ifndef NDEBUG
    for ( int i = 0; i < view_x.size(); ++i ) assert((view_x[i] == 6.-i) && "Wrong value. Index problem in operator [] ?" );
#   endif
    K_MEMORY::vector_view<double> view_y(view_x);
    assert((view_y.data() == view_x.data()) && "view_y would have same memory view as view_x !");
    assert((view_y.size() == view_x.size()) && "view_y would have same size as view_x !");
#   ifndef NDEBUG
    int ind = 0;
    for ( K_MEMORY::vector_view<double>::iterator it = view_y.begin(); it != view_y.end(); ++it, ++ind ) {
      assert(((*it) == view_x[ind]) && "Wrong value. Problem in iterator ?" );
    } 
#   endif
    bool hasCapturedException = false;
    try {
      double x = view_y.at(10);
      std::cout << "Error, may don't execute this line with " << x << std::endl;
    } catch(std::out_of_range& e) {
      hasCapturedException = true;
    }
    assert( hasCapturedException && "None exception throwed with at() method and wrong index ?");
    assert((view_y.at(0) == view_y[0]) && "Wrong value. Problem in at method ?" );
    assert((view_y.front() == view_y[0]) && "Wrong value. Problem in front method ?");
    assert((view_y.back() == view_y[view_y.size()-1]) && "Wrong value. Problem in back method ?");
    view_y.clear();
    assert((view_y.data() == nullptr) && "view_y would have nullptr after clearing !");
    assert((view_y.size() == 0) && "view_y would have size zero after clearing !");
    view_y.swap(view_x);
    assert((view_x.data() == nullptr) && "view_x would have nullptr after swapping !");
    assert((view_x.size() == 0) && "view_x would have size zero after swapping !");
    assert((view_y.data() == &tab[1][1]) && "view_y has wrong memory viewing after swapping !");
    assert((view_y.size() == 5) && "Wrong size of view_y after swapping");
    std::vector<double> y = view_y;
    assert((y.size() == view_y.size()) && "Wrong size for y after converting view_y !");
    assert((&y[0] != view_y.data()) && "y must have different address than view_y for data !!");
#   ifndef NDEBUG
    for ( int i = 0; i < y.size(); ++i ) 
      assert((y[i] == view_y[i]) && "Wrong value. y must have same values than view_y after conversion ?" );
#   endif
  }
  printf("End of memory tests\n");
#endif

#if EXTARITH == 1
  {
    using namespace ExtendedArithmetics;
    printf("Beginning extended arithmetics tests\n");
    double sum0(0.);
    double sum1(0.);
    quad_double qsum0(0.);
    quad_double qsum1(0.);

    quad_double qone = quad_double(1.);
    for ( int i = 1; i <= 65535; i += 2)
    {
        sum0 += 1./double(i);
        qsum0 = qsum0 + qone/double(i);
    }
    for ( int i = 65535; i > 0; i -= 2)
    {
        sum1 += 1./double(i);
        qsum1 = qsum1 + qone/double(i);
    }
    double diff = sum1 - sum0;
    double qdiff= qsum1 - qsum0;
    printf("Erreur somme double : %lg et erreur somme quad_double : %lg\n", diff, double(qdiff));
    assert(std::abs(double(qdiff)) < std::abs(diff));
    double prod0(1.), prod1(1.);
    quad_double qprod0(1.), qprod1(1.);
    for ( int i = 1; i <= 63; i += 2 )
    {
      prod0 *= i;
      qprod0 = qprod0 * double(i);
    }
    for ( int i = 63; i > 0; i -= 2 )
    {
      prod1 *= i;
      qprod1 = qprod1 * double(i);
    }
    diff = prod1 - prod0;
    qdiff = qprod1 - qprod0;
    printf("Erreur produit double : %lg et erreur produit quad_double : %lg\n", diff, double(qdiff));
    assert(std::abs(double(qdiff)) < std::abs(diff));
    double div0(1.), div1(1.);
    quad_double qdiv0(1.), qdiv1(1.);
    for ( int i = 1; i <= 17; i += 2 )
    {
      div0 /= i;
      qdiv0 = qdiv0 / double(i);
    }
    for ( int i = 17; i > 0; i -= 2 )
    {
      div1 /= i;
      qdiv1 = qdiv1 / double(i);
    }
    diff = div1 - div0;
    qdiff = qdiv1 - qdiv0;
    printf("Erreur division double : %lg et erreur division quad_double : %lg\n", diff, double(qdiff));
    assert(std::abs(double(qdiff)) < std::abs(diff));
    double root = std::sqrt(prod0);
    quad_double qroot = sqrt(qprod0);
    printf("Square root of double : %lg et square root of quad_double : %lg\n", root, double(qroot));
  }
  printf("End of Extended arithmetics tests\n");
#endif

#if TESTLOGGER == 1
  {
    K_LOGGER::logger log;
    log.subscribe( new K_LOGGER::log_to_std_output(K_LOGGER::logger::all) );
    log.subscribe( new K_LOGGER::log_from_distributed_file(K_LOGGER::logger::information, "output", 0));
    log << "Hello from test with " << 234 << std::endl;
    log << LogError << "Aie ca fait mal...\n\n\n...\nMais non, je rigole :-)" << std::endl;
    log << LogAssert((1==0)) << "Mmmh, ne devrait pas s'afficher" << std::endl;
  }
#endif  

  Py_INCREF(Py_None);
  return Py_None;
}
