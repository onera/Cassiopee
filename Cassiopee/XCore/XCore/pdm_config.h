#ifndef __PDM_CONFIG_H__
#define __PDM_CONFIG_H__

// #define PDM_HAVE_PARMETIS
#define PDM_HAVE_PTSCOTCH
//#define PDM_HAVE_GETRUSAGE
//#define PDM_HAVE_GETTIMEOFDAY

// CBX
#ifdef G_DOUBLEINT
  #define PDM_LONG_G_NUM
  #define PDM_LONG_G_NUM_BOOL "1"
#else
  #define PDM_LONG_G_NUM_BOOL "0"  
#endif
// END CBX

#define PDM_IN_PDMA_BOOL 0
#define PDM_VERSION_MAJOR "1"
#define PDM_VERSION_MINOR "6"
#define PDM_VERSION_PATCH "2"
#define PDM_VERSION       "1.6.2"

#endif /*__PDM_CONFIG_H__*/
