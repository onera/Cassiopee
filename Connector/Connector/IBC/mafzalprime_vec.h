ax = aa_vec[noind]*utau_vec[noind];
px = gradP_vec[noind]/pow(utau_vec[noind],3.);

l11 = pow(ax+10.6,9.6);
l12 = ax*(ax-8.15) + 86.;
l13 = (2.*ax-8.15)/16.7;
tp =  aa_vec[noind]*( 9.6*pow((ax + 10.6),8.6)  - l11/l12*(  4.*aa_vec[noind]*utau_vec[noind] - 16.30 ) );

l1 = 0.649580838323*aa_vec[noind]/(1. + l13*l13 )  +  tp/(l11*2.302585093); //MUSKER
l2 = 0.;
l3 = 0.;

if (px > 0.){
  if (MafzalMode == 1){
    l2 = 2.*px*ax*kappainv/(utau_vec[noind]*sqrt(1.+px*ax)*(sqrt(1.+px*ax)+1.)); //PRESSURE
    l3 = -2.*px*ax*kappainv/(utau_vec[noind]*sqrt(1.+px*ax)); //PRESSURE
  }
  else{
    l2 = -2.*px*ax*kappainv/(utau_vec[noind]*sqrt(1.+px*ax)*(sqrt(1.+px*ax)+1.)); //PRESSURE
    l3 = 0.; //PRESSURE
  }
}
else{
  if (MafzalMode == 3){
    px = -px;
    l2 = 2.*px*ax*kappainv/(utau_vec[noind]*sqrt(1.+px*ax)*(sqrt(1.+px*ax)+1.)); //PRESSURE
    l3 = 0.; //PRESSURE
  }
  else{
    l2 = 0.; //PRESSURE
    l3 = 0.; //PRESSURE
  }
}

fp = l1 + l2 + l3 + uext_vec[noind]/(utau_vec[noind]*utau_vec[noind]);