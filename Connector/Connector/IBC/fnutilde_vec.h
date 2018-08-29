//cv1cube=7.1^3=  357.911 
l1               = (nutilde/mu_vec[noind])*ro_vec[noind];       //chi
l2               = l1*l1*l1;                                    //chicube
utauv_vec[noind] = l2/(l2+357.911)*nutilde-nutcible_vec[noind]; // delta de nut

if (K_FUNC::E_abs(utauv_vec[noind]) > newtonepsnutilde  && skip == 0) err = 0; 
