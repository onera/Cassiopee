#include "pdm_configf.h"
#ifdef PDM_LONG_G_NUM
  integer, parameter :: pdm_g_num_s = 8
#else
  integer, parameter :: pdm_g_num_s = 4
#endif
  integer, parameter :: pdm_l_num_s = 4

  integer, parameter :: pdm_max_char_length = 100
