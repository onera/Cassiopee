c_fd=coef*rho_[indR];

if((int)cellNIBCx[indR]==0){
  fd_a   = rho_u[l1]/rho_[l1];
  fd_b   = rho_u[l2]/rho_[l2];
  a1fd_1 = c_fd*0.5*(fd_b-fd_a);

  fd_a    = rho_v[l1]/rho_[l1];
  fd_b    = rho_v[l2]/rho_[l2];
  fd_prt1 = 0.5*(fd_b-fd_a);

  fd_a    = rho_w[l1]/rho_[l1];
  fd_b    = rho_w[l2]/rho_[l2];
  fd_prt12= 0.5*(fd_b-fd_a);
 }
 else if((int)cellNIBCx[indR]==1) {
   fd_a   =  rho_u[l  ]/rho_[l  ];
   fd_b   = -rho_u[l2 ]/rho_[l2 ];
   fd_c   =  rho_u[l22]/rho_[l22];
   a1fd_1 = c_fd*0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a    =  rho_v[l  ]/rho_[l  ];
   fd_b    = -rho_v[l2 ]/rho_[l2 ];
   fd_c    =  rho_v[l22]/rho_[l22];
   fd_prt1 = 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a    =  rho_w[l  ]/rho_[l  ];
   fd_b    = -rho_w[l2 ]/rho_[l2 ];
   fd_c    =  rho_w[l22]/rho_[l22];
   fd_prt12= 0.5*(3*fd_a+4*fd_b+fd_c);
 }	
 else if((int)cellNIBCx[indR]==-1) {
   fd_a   = -rho_u[l  ]/rho_[l  ];
   fd_b   =  rho_u[l1 ]/rho_[l1 ];
   fd_c   = -rho_u[l12]/rho_[l12];
   a1fd_1 = c_fd*0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a    = -rho_v[l  ]/rho_[l  ];
   fd_b    =  rho_v[l1 ]/rho_[l1 ];
   fd_c    = -rho_v[l12]/rho_[l12];
   fd_prt1 = 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a    = -rho_w[l  ]/rho_[l  ];
   fd_b    =  rho_w[l1 ]/rho_[l1 ];
   fd_c    = -rho_w[l12]/rho_[l12];
   fd_prt12= 0.5*(3*fd_a+4*fd_b+fd_c);
 }

if((int)cellNIBCy[indR]==0) {
  fd_a    = rho_u[l3]/rho_[l3];
  fd_b    = rho_u[l4]/rho_[l4];
  fd_prt2 = 0.5*(fd_b-fd_a);

  fd_a   = rho_v[l3]/rho_[l3];
  fd_b   = rho_v[l4]/rho_[l4];
  a1fd_5 = c_fd*0.5*(fd_b-fd_a);

  fd_a    = rho_w[l3]/rho_[l3];
  fd_b    = rho_w[l4]/rho_[l4];
  fd_prt22= 0.5*(fd_b-fd_a);
 }
 else if ((int)cellNIBCy[indR]==1) {
   fd_a    =  rho_u[l  ]/rho_[l  ];
   fd_b    = -rho_u[l4 ]/rho_[l4 ];
   fd_c    =  rho_u[l42]/rho_[l42];
   fd_prt2 = 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a   =  rho_v[l  ]/rho_[l  ];
   fd_b   = -rho_v[l4 ]/rho_[l4 ];
   fd_c   =  rho_v[l42]/rho_[l42];
   a1fd_5 = 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a    =  rho_w[l  ]/rho_[l  ];
   fd_b    = -rho_w[l4 ]/rho_[l4 ];
   fd_c    =  rho_w[l42]/rho_[l42];
   fd_prt22= 0.5*(3*fd_a+4*fd_b+fd_c);
 }
 else if ((int)cellNIBCy[indR]==-1) {
   fd_a    = -rho_u[l  ]/rho_[l  ];
   fd_b    =  rho_u[l3 ]/rho_[l3 ];
   fd_c    = -rho_u[l32]/rho_[l32];
   fd_prt2 = 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a   = -rho_v[l  ]/rho_[l  ];
   fd_b   =  rho_v[l3 ]/rho_[l3 ];
   fd_c   = -rho_v[l32]/rho_[l32];
   a1fd_5 = c_fd*0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a    = -rho_w[l  ]/rho_[l  ];
   fd_b    =  rho_w[l3 ]/rho_[l3 ];
   fd_c    = -rho_w[l32]/rho_[l32];
   fd_prt22= 0.5*(3*fd_a+4*fd_b+fd_c);
 }
                   
if((int)cellNIBCz[indR]==0) {
  fd_a    = rho_u[l5]/rho_[l5];
  fd_b    = rho_u[l6]/rho_[l6];
  fd_prt3 = 0.5*(fd_b-fd_a);

  fd_a    = rho_v[l5]/rho_[l5];
  fd_b    = rho_v[l6]/rho_[l6];
  fd_prt32= 0.5*(fd_b-fd_a);

  fd_a   = rho_w[l5]/rho_[l5];
  fd_b   = rho_w[l6]/rho_[l6];
  a1fd_9 = c_fd*0.5*(fd_b-fd_a);
 }
 else if ((int)cellNIBCz[indR]==1) {
   fd_a    =  rho_u[l  ]/rho_[l  ];
   fd_b    = -rho_u[l6 ]/rho_[l6 ];
   fd_c    =  rho_u[l62]/rho_[l62];
   fd_prt3 = 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a    =  rho_v[l  ]/rho_[l  ];
   fd_b    = -rho_v[l6 ]/rho_[l6 ];
   fd_c    =  rho_v[l62]/rho_[l62];
   fd_prt32= 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a   =  rho_w[l  ]/rho_[l  ];
   fd_b   = -rho_w[l6 ]/rho_[l6 ];
   fd_c   =  rho_w[l62]/rho_[l62];
   a1fd_9 = c_fd*0.5*(3*fd_a+4*fd_b+fd_c);
 }
 else if ((int)cellNIBCz[indR]==-1) {
   fd_a    = -rho_u[l  ]/rho_[l  ];
   fd_b    =  rho_u[l5 ]/rho_[l5 ];
   fd_c    = -rho_u[l52]/rho_[l52];
   fd_prt3 = 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a    = -rho_v[l  ]/rho_[l  ];
   fd_b    =  rho_v[l5 ]/rho_[l5 ];
   fd_c    = -rho_v[l52]/rho_[l52];
   fd_prt32= 0.5*(3*fd_a+4*fd_b+fd_c);
                   
   fd_a   = -rho_w[l  ]/rho_[l  ];
   fd_b   =  rho_w[l5 ]/rho_[l5 ];
   fd_c   = -rho_w[l52]/rho_[l52];
   a1fd_9 = c_fd*0.5*(3*fd_a+4*fd_b+fd_c);
 }
a1fd_2 = c_fd*0.5*(fd_prt2+fd_prt1 );
a1fd_3 = c_fd*0.5*(fd_prt3+fd_prt12);
a1fd_4 = a1fd_2;
a1fd_6 = c_fd*0.5*(fd_prt22+fd_prt32);
                   
a1fd_7=a1fd_3;
a1fd_8=a1fd_6;
     
a1hrr_1 = um_sigma*a1fd_1;
a1hrr_2 = um_sigma*a1fd_2;
a1hrr_3 = um_sigma*a1fd_3;
a1hrr_4 = um_sigma*a1fd_4;
a1hrr_5 = um_sigma*a1fd_5;
a1hrr_6 = um_sigma*a1fd_6;
a1hrr_7 = um_sigma*a1fd_7;
a1hrr_8 = um_sigma*a1fd_8;
a1hrr_9 = um_sigma*a1fd_9;
