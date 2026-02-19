# include "intersector.h"
# include <string>
# include <sstream> 


# include "Nuga/include/ngon_t.hxx"

//#include <iostream>
#include <memory>

using namespace std;
using namespace K_FLD;


PyObject* K_INTERSECTOR::testmain(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  //K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;

  ngon_type ngi(cnt);

  // todo
  
  delete f; delete cn;
  return NULL;
}
