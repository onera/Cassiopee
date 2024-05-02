        // Newton pour mut - out : aa_vec
        count = 0;
        while (err == 0 && count < 30)
        {
         err = 1;
#ifdef _OPENM4
         #pragma omp simd
#endif 
         for (E_Int noind = 0; noind < ifin-ideb; noind++)
           {
             if (K_FUNC::E_abs(utauv_vec[noind]) > newtonepsnutilde)
              {
                skip = 0; 
                //fp
                nutilde  = aa_vec[noind];
#               include "IBC/fnutildeprime_vec.h"
                if (K_FUNC::E_abs(f1p) >= newtonepsprime) { aa_vec[noind] = K_FUNC::E_abs(nutilde-utauv_vec[noind]/f1p); }
                else                                      {  skip = 1; }

                nutilde  =  aa_vec[noind];
#               include "IBC/fnutilde_vec.h"
              }
           }         //loop point
         count++;
        }            //loop newton

       if (count == 30) // blindage si newton ne converge pas: utau = utau0
        {
#ifdef _OPENMP4
          #pragma omp simd
#endif 
          for (E_Int noind = 0; noind < ifin-ideb; noind++)  // test pour imiter stephanie
            { if (K_FUNC::E_abs(utauv_vec[noind]) > newtonepsnutilde) aa_vec[noind]  = ut_vec[noind];  
            }
        }
