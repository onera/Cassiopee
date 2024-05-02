#     include "IBC/commonLaws1.h" 
      // loi de Musker 
      press_vec[noind]     = pext;
      yplus_vec[noind]     = roext*yibc/muext;
      aa_vec[noind]       = roext*yext/muext;
      uext_vec[noind]     = uext;
      utau_vec[noind]     = utau0;
      nutcible_vec[noind] = utau0;
      sign_vec[noind]     = signibc/uext;
      ucible_vec[noind]   = alphasbeta*un;
      vcible_vec[noind]   = alphasbeta*vn;
      wcible_vec[noind]   = alphasbeta*wn;
      ut_vec[noind]       = ut;
      vt_vec[noind]       = vt;
      wt_vec[noind]       = wt;
      ro_vec[noind]       = roext;
      mu_vec[noind]       = muext;
      alpha_vec[noind]    = alpha;
      tcible_vec[noind]   = text;
#     include "IBC/logf_vec.h"
      // out= utau  et err 
