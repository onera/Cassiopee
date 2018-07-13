// Newton pour utau
E_Int count = 0;
while (err == 0 && count < 50)
{
  err = 1;
#ifdef _OPENMP4
#pragma omp simd
#endif 
  for (E_Int noind = 0; noind < ifin-ideb; noind++)  // test pour imiter stephanie
  {
    if (K_FUNC::E_abs(utauv_vec[noind]) > newtoneps)
    {
      skip = 0; 
      //fp
# include "IBC/muskerprime_vec.h"
      if (K_FUNC::E_abs(fp) >= newtonepsprime) 
      { 
        utau_vec[noind] = K_FUNC::E_abs(utau_vec[noind]-utauv_vec[noind]/fp); 
      }
      else
      {
        skip = 1; 
      } 
      
# include "IBC/musker_vec.h"
    }
  }//loop point
  count++;
}  //loop newton

//for (E_Int noind = 0; noind < ifin-ideb; noind++)
//  printf("newton %d %d %f\n", count, noind, utau_vec[noind]);

if (count == 50) // blindage si au moins un pt newton ne converge pas: utau = utau0
{
#ifdef _OPENMP4
#pragma omp simd
#endif 
  for (E_Int noind = 0; noind < ifin-ideb; noind++)  // test pour imiter stephanie
  { 
    if (K_FUNC::E_abs(utauv_vec[noind]) > newtoneps) utau_vec[noind] = nutcible_vec[noind];  
  }
}
