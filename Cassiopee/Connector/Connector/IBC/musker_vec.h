ax = aa_vec[noind]*utau_vec[noind];  
l1 = pow(ax+10.6,9.6);
l2 = ax*(ax-8.15) + 86.;
l3 = (2.*ax-8.15)/16.7;

utauv_vec[noind] = 5.424*atan(l3) + log10(l1/(l2*l2)) - uext_vec[noind]/utau_vec[noind] - 3.52;

if (K_FUNC::E_abs( utauv_vec[noind] ) > newtoneps && skip == 0) err = 0;
