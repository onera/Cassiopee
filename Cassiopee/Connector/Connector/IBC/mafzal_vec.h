ax = aa_vec[noind]*utau_vec[noind];
px = gradP_vec[noind]/pow(utau_vec[noind],3);

l11 = pow(ax+10.6,9.6);
l12 = ax*(ax-8.15) + 86.;
l13 = (2.*ax-8.15)/16.7;

l1 = 5.424*atan(l13) + log10(l11/(l12*l12)) - 3.52; //MUSKER
l2 = 0.;
l3 = 0.;

if (px > 0.){
  if (MafzalMode == 1){
    l2 = -2.*kappainv*log((sqrt(1.+px*ax)+1.)/2.); //PRESSURE
    l3 = 2.*kappainv*(sqrt(1.+px*ax)-1.); //PRESSURE
  }
  else{
    l2 = 2.*kappainv*log((sqrt(1.+px*ax)+1.)/2.); //PRESSURE
    l3 = 0.; //PRESSURE
  }
}
else{
  if (MafzalMode == 3){
    px = -px;
    l2 = -2.*kappainv*log((sqrt(1.+px*ax)+1.)/2.); //PRESSURE
    l3 = 0.; //PRESSURE
  }
  else{
    l2 = 0.; //PRESSURE
    l3 = 0.; //PRESSURE
  }
}

utauv_vec[noind] = l1 + l2 + l3 - uext_vec[noind]/utau_vec[noind];

if ((K_FUNC::E_abs(utauv_vec[noind]) > newtoneps) || (utauv_vec[noind] != utauv_vec[noind])) err = 0;
