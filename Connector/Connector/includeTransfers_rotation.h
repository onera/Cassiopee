// Extraction de l angle de rotation
E_Int dirR = 0; E_Float theta = 0.;
if      (K_FUNC::E_abs(angle[0]) > 0.){dirR=1; theta=angle[0];}
else if (K_FUNC::E_abs(angle[1]) > 0.){dirR=2; theta=angle[1];}
else if (K_FUNC::E_abs(angle[2]) > 0.){dirR=3; theta=angle[2];} 
// E_Float thetar = theta*K_CONST::E_PI/180.;// theta en radians pour sin/cos
E_Float ctheta = cos(theta); E_Float stheta = sin(theta);

//printf(" teta= %f, ithread= %d \n", theta , ithread ); 
E_Float v1, v2;

switch (dirR)
{
  case 1:  //  nuage de pts quelconque
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];

      v1 = vectOfRcvFields[2][indR]; //vy
      v2 = vectOfRcvFields[3][indR]; //vz

      vectOfRcvFields[2][indR] = ctheta*v1-stheta*v2;
      vectOfRcvFields[3][indR] = stheta*v1+ctheta*v2;
    }
    break;

  case 2:
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];
    
      v1 = vectOfRcvFields[1][indR]; //vx
      v2 = vectOfRcvFields[3][indR]; //vz

      vectOfRcvFields[1][indR] = ctheta*v1+stheta*v2;
      vectOfRcvFields[3][indR] =-stheta*v1+ctheta*v2;

    }
    break;
    
  case 3: // Structure Lineaire O2 par tetra
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];

      v1 = vectOfRcvFields[1][indR]; //vx
      v2 = vectOfRcvFields[2][indR]; //vy

      vectOfRcvFields[1][indR] = ctheta*v1-stheta*v2;
      vectOfRcvFields[2][indR] = stheta*v1+ctheta*v2;
    }
    break;
    
  default: ;
}
