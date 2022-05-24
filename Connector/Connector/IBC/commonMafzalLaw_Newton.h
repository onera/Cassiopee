// Newton pour utau
// E_Int count = 0;
count = 0;

while (err == 0 && count < 500)
{
  err = 1;
#ifdef _OPENMP4
#pragma omp simd
#endif
  for (E_Int noind = 0; noind < ifin-ideb; noind++)  // test pour imiter stephanie
  {
    if ((K_FUNC::E_abs(utauv_vec[noind]) > newtoneps) || (utauv_vec[noind] != utauv_vec[noind]))
    {
# include "IBC/mafzalprime_vec.h"
      utau_vec[noind] = K_FUNC::E_abs(utau_vec[noind]-utauv_vec[noind]/fp);
# include "IBC/mafzal_vec.h"
    }
  }//loop point
  count++;
}  //loop newton

if (count == 500) // blindage si au moins un pt newton ne converge pas: utau = utau0
{
#ifdef _OPENMP4
#pragma omp simd
#endif
  for (E_Int noind = 0; noind < ifin-ideb; noind++)  // test pour imiter stephanie
  {
    // si pas de convergence, on utilise la solution de musker
    if ((K_FUNC::E_abs(utauv_vec[noind]) > newtoneps) || (utauv_vec[noind] != utauv_vec[noind])){
      utau_vec[noind] = utauOri_vec[noind];
      gradP_vec[noind] = 0.;
    }
  }
}
