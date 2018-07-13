l1 = nutilde*nutilde*nutilde;                  //nutilde^3
l2 = l1*l1;                                    //nutile^6
ax = 7.1*mu_vec[noind]/ro_vec[noind];          //fnu = cv1*muext/roext;
l3 = l1 + ax*ax*ax;            //g
//  E_Float g2 = g*g;
//  E_Float fnu = 7.1*muext/roext;
//  E_Float fp = 4*l1;                              
//  E_Float gp = 3*nutilde*nutilde;

f1p =  (4.*l1*l3 - 3.*l2)/(l3*l3);
