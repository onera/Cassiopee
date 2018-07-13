utauv_vec[noind] = kappainv*log( aa_vec[noind]*utau_vec[noind]) - uext_vec[noind]/utau_vec[noind] + cc;

if (K_FUNC::E_abs( utauv_vec[noind] ) > newtoneps && skip == 0) err = 0;
