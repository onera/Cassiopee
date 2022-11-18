#ifndef __PDM_CONFIG_H__
#define __PDM_CONFIG_H__

// #define PDM_HAVE_PARMETIS
#define PDM_HAVE_PTSCOTCH
//#define PDM_HAVE_GETRUSAGE
//#define PDM_HAVE_GETTIMEOFDAY

// CBX
#ifdef E_DOUBLEINT
  #define PDM_LONG_G_NUM
  #define PDM_LONG_G_NUM_BOOL "1"

  //#define __G_NPY_T__ NPY.int64_t
  //#define __L_NPY_T__ NPY.int32_t
  //#define __G_T__ int64_t
  //#define __L_T__ int
  //#define __G_NPY_ENUM__ NPY.NPY_INT64
#else
  #define PDM_LONG_G_NUM_BOOL "0"
  
  //#define __G_NPY_T__ NPY.int32_t
  //#define __L_NPY_T__ NPY.int32_t
  //#define __G_T__ int32_t
  //#define __L__T__ int
  //#define __G_NPY_ENUM__ NPY.NPY_INT32
#endif
// END CBX

#define PDM_IN_PDMA_BOOL 0
#define PDM_VERSION_MAJOR "1"
#define PDM_VERSION_MINOR "6"
#define PDM_VERSION_PATCH "2"
#define PDM_VERSION       "1.6.2"

#endif /*__PDM_CONFIG_H__*/
