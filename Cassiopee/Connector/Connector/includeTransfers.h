  // E_Float thetar = theta*K_CONST::E_PI/180.;// theta en radians pour sin/cos
  E_Float ctheta = cos(theta); E_Float stheta = sin(theta);
  if ( posvx > -1 && posvy > -1 && posvz > -1)
  {
    E_Float* ptrVx = vectOfRcvFields[posvx];
    E_Float* ptrVy = vectOfRcvFields[posvy];
    E_Float* ptrVz = vectOfRcvFields[posvz];
    for (E_Int noind = 0; noind < nbRcvPts; noind++)
    { 
      E_Int indR = rcvPts[noind];
      E_Float v1, v2;
      if ( dirR == 1)
      {
        v1 = ptrVy[indR]; v2 = ptrVz[indR];
        ptrVy[indR] = ctheta*v1-stheta*v2;
        ptrVz[indR] = stheta*v1+ctheta*v2;
      }
      else if ( dirR == 2)
      {
        v1 = ptrVx[indR]; v2 = ptrVz[indR];
        ptrVx[indR] = ctheta*v1+stheta*v2;
        ptrVz[indR] =-stheta*v1+ctheta*v2;
      }
      else if (dirR == 3)
      {
       v1 = ptrVx[indR]; v2 = ptrVy[indR];
       ptrVx[indR] = ctheta*v1-stheta*v2;
       ptrVy[indR] = stheta*v1+ctheta*v2; 
      }
    }
  }
  if ( posmx > -1 && posmy > -1 && posmz > -1)
  {
    E_Float* ptrVx = vectOfRcvFields[posmx];
    E_Float* ptrVy = vectOfRcvFields[posmy];
    E_Float* ptrVz = vectOfRcvFields[posmz];
    for (E_Int noind = 0; noind < nbRcvPts; noind++)
    { 
      E_Int indR = rcvPts[noind];
      E_Float v1, v2;
      if ( dirR == 1)
      {
        v1 = ptrVy[indR]; v2 = ptrVz[indR];
        ptrVy[indR] = ctheta*v1-stheta*v2;
        ptrVz[indR] = stheta*v1+ctheta*v2;
      }
      else if ( dirR == 2)
      {
        v1 = ptrVx[indR]; v2 = ptrVz[indR];
        ptrVx[indR] = ctheta*v1+stheta*v2;
        ptrVz[indR] =-stheta*v1+ctheta*v2;
      }
      else if (dirR == 3)
      {
       v1 = ptrVx[indR]; v2 = ptrVy[indR];
       ptrVx[indR] = ctheta*v1-stheta*v2;
       ptrVy[indR] = stheta*v1+ctheta*v2; 
     }
   }
 }
