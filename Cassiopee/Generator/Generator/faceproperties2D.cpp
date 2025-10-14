      xt1 = structF[v1]->begin(posxt[v1]);
      yt1 = structF[v1]->begin(posyt[v1]);
      zt1 = structF[v1]->begin(poszt[v1]);
      ni1 = nit[v1]; nj1 = njt[v1]; nk1 = nkt[v1]; 
      indA1 = 0; indB1 = ni1-1; indD1 = (nj1-1)*ni1; indC1 = indB1 + indD1;

      if (njt[v1] !=1){ s1 = (ymaxp[v1]-yminp[v1])/(njt[v1]-1);}
      else s1= 1.e-12; //en 2d existance de zone lineique effacerr plus tard. on met cutoo pour eviter / zero

      ipt_dx[v1]=s1;

      /* facette AD ou i = 1 */
      minB[0] = xt1[indA1]; minB[1] = yt1[indA1]; minB[2] = zt1[indA1];
      maxB[0] = xt1[indD1]; maxB[1] = yt1[indD1]; maxB[2] = zt1[indD1];

      // elts intersectant la facette i = 1
      indicesBB.clear();
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
      dhmin = 1.e10;// dh min des grilles adjacentes

      found1 = 0; found2 = 0; found3 = 0; found4 = 0;
      // facette opposee en i = imax : B'C'G'F'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        v2 = indicesBB[noe]; 
        printf("v1: %d , v2: %d ; candidat: %d nijk1: %d %d %d nijk2: %d %d %d , dx1: %f \n",  v1, v2, noe,  ni1, nj1, nk1,  nit[v2], njt[v2], nkt[v2], s1);
        if (v2 == v1 ) 
          goto finimin2d;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        face(0,1) = xt2[indB2]; face(0,2) = yt2[indB2]; face(0,3) = zt2[indB2];
        face(1,1) = xt2[indC2]; face(1,2) = yt2[indC2]; face(1,3) = zt2[indC2];

        //printf("posX %f \n", xt1[indA1]-xt2[indB2]);

        if ( K_FUNC::fEqualZero(xt1[indA1]-xt2[indB2],tol) == false ) goto finimin2d;

        if (njt[v2] !=1){ s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1);}
        else s2= 1.e-12; //en 2d existance de zone lineique effacerr plus tard. on met cutoo pour eviter / zero
        
        ymax2 = yt2[indC2]; ymin2 = yt2[indA2]; ymax1 = yt1[indC1]; ymin1 = yt1[indA1]; 
        
        if ( ymax2 > ymin1  && ymin2 < ymax1 && 
            K_FUNC::fEqualZero(ymax2-ymin1,tol) == false && K_FUNC::fEqualZero(ymax1-ymin2,tol) == false
            )
        { iptvois[v1*6*mx_sf+ 0*mx_sf +found1]=v2;
          found1 +=1; 
          dhmax = K_FUNC::E_max(dhmax,s2); dhmin = K_FUNC::E_min(dhmin,s2);
          printf("imin: Nr: %d , nijkR: %d %d %d , s2: %f \n", v2, ni2, nj2, nk2, s2); }

        finimin2d:;
          
      }// fin parcours de ts les elts intersectant
       
      if ( found1 > 0) 
        {  diff = s1-dhmax; diff1 = dhmax-dhmin;
           E_Int res=-100;
           if      (dhmax/dhmin> 2.5 && dhmax/s1 > 1.5) res =3;
           else if (dhmax/dhmin> 1.5)
               {if     (dhmax/s1 > 1.4) res = 2;
                else if(dhmin/s1 < 0.7) res =-2;
               }
           else
              { if      (dhmax/s1 > 1.4) res = 2;
                else if (dhmax/s1 < 0.7) res =-2;
                else if (dhmax/s1 < 1.1 && dhmax/s1 > 0.9) res = 1;
              }

           sf[v1*6 + 0]   = found1;
           resol[v1*6 +0]= res;
           printf("imin: Nd: %d , sous face: %d , Resol: %d , diff: %f Hmax/s1 %f , Hmin/s1: %f , hmax/hmin: %f \n", v1, found1, res, diff, dhmax/s1, dhmin/s1, dhmax/dhmin);
        }
      // fin test facette i = 1

      /* facette BC ou i = imax */
      //faceimax:
      minB[0] = xt1[indB1]; minB[1] = yt1[indB1]; minB[2] = zt1[indB1];
      maxB[0] = xt1[indC1]; maxB[1] = yt1[indC1]; maxB[2] = zt1[indC1];
      // elts intersectant la facette i = imax
      indicesBB.clear(); 
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.; dhmin=1.e10;// dh max des grilles adjacentes

      found1 = 0; found2 = 0; found3 = 0; found4 = 0;
      // facette opposee en i = 1: A'D'H'E'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        v2 = indicesBB[noe]; 
        if ( v2 == v1 ) 
          goto finimax2d;

        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        face(0,1) = xt2[indA2]; face(0,2) = yt2[indA2]; face(0,3) = zt2[indA2];
        face(1,1) = xt2[indD2]; face(1,2) = yt2[indD2]; face(1,3) = zt2[indD2];

        if ( K_FUNC::fEqualZero(xt1[indB1]-xt2[indA2],tol) == false ) goto finimax2d;
        
        if (njt[v2] !=1){ s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1);}
        else s2= 1.e-12; //en 2d existance de zone lineique effacerr plus tard. on met cutoo pour eviter / zero
        
        ymax2 = yt2[indD2]; ymin2 = yt2[indA2]; ymax1 = yt1[indC1]; ymin1 = yt1[indA1];

          printf("Yfen %f %f \n",  ymax2-ymin1, ymax1-ymin2 );

        if ( ymax2 > ymin1  && ymin2 < ymax1 &&
            K_FUNC::fEqualZero(ymax2-ymin1,tol) == false && K_FUNC::fEqualZero(ymax1-ymin2,tol) == false
           )
        { iptvois[v1*6*mx_sf+ 1*mx_sf +found1]=v2;
          found1 += 1;
          dhmax = K_FUNC::E_max(dhmax,s2); dhmin = K_FUNC::E_min(dhmin,s2);
          printf("imax: Nr: %d , nijkR: %d %d %d , s2: %f \n", v2, ni2, nj2, nk2, s2); }

        finimax2d:;
          
      }// fin parcours de ts les elts intersectant 
      if ( found1 > 0) 
        {  diff = s1-dhmax; diff1 = dhmax-dhmin;
           E_Int res=-100;
           if      (dhmax/dhmin> 2.5 && dhmax/s1 > 1.5) res =3;
           else if (dhmax/dhmin> 1.5)
               {if     (dhmax/s1 > 1.4) res = 2;
                else if(dhmin/s1 < 0.7) res =-2;
               }
           else
              { if      (dhmax/s1 > 1.4) res = 2;
                else if (dhmax/s1 < 0.7) res =-2;
                else if (dhmax/s1 < 1.1 && dhmax/s1 > 0.9) res = 1;
              }

           sf[v1*6 + 1]   = found1;
           resol[v1*6 +1]= res;
           printf("imax: Nd: %d , sous face: %d , Resol: %d , diff: %f Hmax/s1 %f , Hmin/s1: %f , hmax/hmin: %f \n", v1, sf[v1*6+1], res, diff, dhmax/s1, dhmin/s1, dhmax/dhmin);
        }
      // fin test facette i = imax

      s1 = (xmaxp[v1]-xminp[v1])/(nit[v1]-1);
      /* facette j =jmin ou AB */
      //facejmin:;    
      minB[0] = xt1[indA1]; minB[1] = yt1[indA1]; minB[2] = zt1[indA1];
      maxB[0] = xt1[indB1]; maxB[1] = yt1[indB1]; maxB[2] = zt1[indB1];

      // elts intersectant la facette j=1
      indicesBB.clear(); 
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
      dhmin = 1.e10;// dh max des grilles adjacentes
      found1 = 0; found2 = 0; found3 = 0; found4 = 0;
      // facette opposee en j = jmax: D'C'G'H'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        v2 = indicesBB[noe]; 
        if ( v2 == v1 ) 
          goto finjmin2d;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        face(0,1) = xt2[indD2]; face(0,2) = yt2[indD2]; face(0,3) = zt2[indD2];
        face(1,1) = xt2[indC2]; face(1,2) = yt2[indC2]; face(1,3) = zt2[indC2];

        if ( K_FUNC::fEqualZero(yt1[indA1]-yt2[indD2],tol) == false ) goto finjmin2d;
        
        if (njt[v2] !=1){ s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1);}
        else s2= 1.e-12; //en 2d existance de zone lineique effacerr plus tard. on met cutoo pour eviter / zero

        xmax2 = xt2[indC2]; xmin2 = xt2[indD2]; xmax1 = xt1[indC1]; xmin1 = xt1[indD1];

        if ( xmax2 > xmin1  && xmin2 < xmax1 &&
            K_FUNC::fEqualZero(xmax2-xmin1,tol) == false && K_FUNC::fEqualZero(xmax1-xmin2,tol) == false
           )
        { iptvois[v1*6*mx_sf+ 2*mx_sf +found1]=v2;
          found1 += 1;
          dhmax = K_FUNC::E_max(dhmax,s2); dhmin = K_FUNC::E_min(dhmin,s2);
          printf("jmin: Nr: %d , nijkR: %d %d %d , s2: %f \n", v2, ni2, nj2, nk2, s2); }


        finjmin2d:;
          // goto facejmax;
      }// fin parcours de ts les elts intersectant 
      if ( found1 > 0) 
        {  diff = s1-dhmax; diff1 = dhmax-dhmin;
           E_Int res=-100;
           if      (dhmax/dhmin> 2.5 && dhmax/s1 > 1.5) res =3;
           else if (dhmax/dhmin> 1.5)
               {if     (dhmax/s1 > 1.4) res = 2;
                else if(dhmin/s1 < 0.7) res =-2;
               }
           else
              { if      (dhmax/s1 > 1.4) res = 2;
                else if (dhmax/s1 < 0.7) res =-2;
                else if (dhmax/s1 < 1.1 && dhmax/s1 > 0.9) res = 1;
              }

           sf[v1*6 + 2]   = found1;
           resol[v1*6 +2]= res;
           printf("jmin: Nd: %d , sous face: %d , Resol: %d , diff: %f Hmax/s1 %f , Hmin/s1: %f , hmax/hmin: %f \n", v1, sf[v1*6+2], res, diff, dhmax/s1, dhmin/s1, dhmax/dhmin);
        }
      // fin test facette j=1
      
      //facejmax:;
      //facette j =jmax ou DC
      minB[0] = xt1[indD1]; minB[1] = yt1[indD1]; minB[2] = zt1[indD1];
      maxB[0] = xt1[indC1]; maxB[1] = yt1[indC1]; maxB[2] = zt1[indC1];

      // elts intersectant la facette j=jmax
      indicesBB.clear(); 
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
      dhmin = 1.e10;// dh max des grilles adjacentes
      // facette opposee en i = 1: A'B'F'E'
      found1 = 0; found2 = 0; found3 = 0; found4 = 0;
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        v2 = indicesBB[noe]; 
        if ( v2 == v1 ) 
         goto finjmax2d;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        face(0,1) = xt2[indA2]; face(0,2) = yt2[indA2]; face(0,3) = zt2[indA2];
        face(1,1) = xt2[indB2]; face(1,2) = yt2[indB2]; face(1,3) = zt2[indB2];

        if ( K_FUNC::fEqualZero(yt1[indD1]-yt2[indA2],tol) == false ) goto finjmax2d;
        
        if (njt[v2] !=1){ s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1);}
        else s2= 1.e-12; //en 2d existance de zone lineique effacerr plus tard. on met cutoo pour eviter / zero

        xmax2 = xt2[indB2]; xmin2 = xt2[indA2]; xmax1 = xt1[indC1]; xmin1 = xt1[indD1];

        if ( xmax2 > xmin1  && xmin2 < xmax1  &&
            K_FUNC::fEqualZero(xmax2-xmin1,tol) == false && K_FUNC::fEqualZero(xmax1-xmin2,tol) == false
           )
        { iptvois[v1*6*mx_sf+ 3*mx_sf +found1]=v2;
          found1 += 1;
          dhmax = K_FUNC::E_max(dhmax,s2); dhmin = K_FUNC::E_min(dhmin,s2);
          printf("jmax: Nr: %d , nijkR: %d %d %d , s2: %f \n", v2, ni2, nj2, nk2, s2); }

        finjmax2d:;
          // goto facekmin;
      }// fin parcours de ts les elts intersectant 
      if ( found1 > 0) 
        {  diff = s1-dhmax; diff1 = dhmax-dhmin;
           E_Int res=-100;
           if      (dhmax/dhmin> 2.5 && dhmax/s1 > 1.5) res =3;
           else if (dhmax/dhmin> 1.5)
               {if     (dhmax/s1 > 1.4) res = 2;
                else if(dhmin/s1 < 0.7) res =-2;
               }
           else
              { if      (dhmax/s1 > 1.4) res = 2;
                else if (dhmax/s1 < 0.7) res =-2;
                else if (dhmax/s1 < 1.1 && dhmax/s1 > 0.9) res = 1;
              }

           sf[v1*6 + 3]   = found1;
           resol[v1*6 +3]= res;
           printf("jmax: Nd: %d , sous face: %d , Resol: %d , diff: %f Hmax/s1 %f , Hmin/s1: %f , hmax/hmin: %f \n", v1, sf[v1*6+3], res, diff, dhmax/s1, dhmin/s1, dhmax/dhmin);
        }
        
      // fin test facette j=jmax
