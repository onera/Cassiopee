ax = aa_vec[noind]*utau_vec[noind];  
l1 = pow(ax+10.6,9.6);
l2 = ax*(ax-8.15) + 86.;
l3 = (2.*ax-8.15)/16.7;

tp =  aa_vec[noind]*( 9.6*pow((ax + 10.6),8.6)  - l1/l2*(  4.*aa_vec[noind]*utau_vec[noind] - 16.30 ) );

//fp =  5.424*2./16.7*aa/(1. + ((2.*ax - 8.15)/16.7)*((2.*ax - 8.15)/16.7) )+ tp/(l1*log(10.)) + bb/(utau*utau);

fp = 0.649580838323*aa_vec[noind]/(1. + l3*l3 )  +  tp/(l1*2.302585093) + uext_vec[noind]/(utau_vec[noind]*utau_vec[noind]);
