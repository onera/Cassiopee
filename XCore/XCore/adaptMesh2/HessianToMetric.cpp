#include "Proto.h"

// TODO(Imad): metric can be constructed in different ways

PyObject *K_XCORE::hessianToMetric(PyObject *self, PyObject *args)
{
  PyObject *HESSIAN;
  E_Float hmin, hmax, eps;
  if (!PYPARSETUPLE_(args, O_ RRR_, &HESSIAN, &hmin, &hmax, &eps)) {
    RAISE("Wrong input.");
    return NULL;
  }

  E_Int nfld, size;
  E_Float *H = NULL;
  K_NUMPY::getFromNumpyArray(HESSIAN, H, size, nfld, true);
  //assert(ret == 1 && nfld == 1);
  E_Int ncells = size / 6;
  assert(H);

  npy_intp dims[2];
  dims[0] = (npy_intp)(ncells*6);
  dims[1] = 1;

  PyArrayObject *METRIC = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  E_Float *M = (E_Float *)PyArray_DATA(METRIC);

  //E_Float cd_by_eps = 0.28125 / eps;
  //E_Float lmin = 1.0/(hmax*hmax);
  //E_Float lmax = 1.0/(hmin*hmin);

  for (E_Int i = 0; i < ncells; i++) {
    E_Float *pH = &H[6*i];

    E_Float L[3], v0[3], v1[3], v2[3];
    K_LINEAR::sym3mat_eigen(pH, L, v0, v1, v2);
    
    //E_Float L0 = sqrt(std::min(std::max(cd_by_eps*fabs(L[0]), lmin), lmax));
    //E_Float L1 = sqrt(std::min(std::max(cd_by_eps*fabs(L[1]), lmin), lmax));
    //E_Float L2 = sqrt(std::min(std::max(cd_by_eps*fabs(L[2]), lmin), lmax));
    
   
    E_Float L0 = sqrt(fabs(L[0]));
    E_Float L1 = sqrt(fabs(L[1]));
    E_Float L2 = sqrt(fabs(L[2]));

    E_Float *pM = &M[6*i];
    pM[0] = L0*v0[0]*v0[0] + L1*v1[0]*v1[0] + L2*v2[0]*v2[0];
    pM[1] = L0*v0[0]*v0[1] + L1*v1[0]*v1[1] + L2*v2[0]*v2[1];
    pM[2] = L0*v0[0]*v0[2] + L1*v1[0]*v1[2] + L2*v2[0]*v2[2];
    pM[3] = L0*v0[1]*v0[1] + L1*v1[1]*v1[1] + L2*v2[1]*v2[1];
    pM[4] = L0*v0[1]*v0[2] + L1*v1[1]*v1[2] + L2*v2[1]*v2[2];
    pM[5] = L0*v0[2]*v0[2] + L1*v1[2]*v1[2] + L2*v2[2]*v2[2];
  }

  return (PyObject *)METRIC;
}
