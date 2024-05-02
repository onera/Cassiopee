      
        npass = npass + 1;
         
        for (E_Int j = 0; j < nbptslinelets; ++j)
        {
          matp[j] = 0.0;
          matm[j] = 0.0;
          mat[j]  = 0.0;
        }
        
        mat[0] = 1.0;
        mat[nbptslinelets-1] = 1.0;

        ipt_u1dold[0]               = 0;
        ipt_u1dold[nbptslinelets-1] = uext_vec[noind];

//       Remplissage matrice 

        #pragma omp simd
        for (E_Int j = 1; j < nbptslinelets-1; ++j)
        {
//         Cell Vertex
//         calcul coordonnées des noeuds de la cellule
//         et des coeffs pour les fluxs
          ynm = 0.5*(yline1d[j-1] + yline1d[j]);
          ynp = 0.5*(yline1d[j]   + yline1d[j+1]);
          dy = ynp - ynm;
          dym = yline1d[j] - yline1d[j-1];
          dyp = yline1d[j+1] - yline1d[j];

          xim = pow((nutilde1d[j-1]/nu),3);
          xi  = pow((nutilde1d[j  ]/nu),3);
          xip = pow((nutilde1d[j+1]/nu),3);

          nutrm = nutilde1d[j-1]*xim/(xim + Cv1cube);
          nutr  = nutilde1d[j  ]*xi/(xi + Cv1cube);
          nutrp = nutilde1d[j+1]*xip/(xip + Cv1cube);

          nutrm=max(nutrm,0.0);
          nutr=max(nutr,0.0);
          nutrp=max(nutrp,0.0);
          
          nutm = 0.5*(nutrm + nutr);
          nutp = 0.5*(nutr  + nutrp);
//         Schéma centré ordre 2
          mat[j] = - (nu + nutm)/(dy*dym)
                   - (nu + nutp)/(dy*dyp);          

          matm[j] = (nu + nutm)/(dy*dym);
          matp[j] = (nu + nutp)/(dy*dyp);

          // Source term
          ipt_u1dold[j] = (gradP_vec[noind] + psild[j]*conv_vec[noind]);  
        }

        for (E_Int j = 1; j < nbptslinelets; ++j)
        {
          m      = matm[j]/mat[j-1];
          mat[j] = mat[j] - m*matp[j-1];         
          ipt_u1dold[j] = ipt_u1dold[j] - m*ipt_u1dold[j-1]; 
        }

        ipt_u1dold[nbptslinelets-1] = ipt_u1dold[nbptslinelets-1]/mat[nbptslinelets-1];

        for (E_Int j = nbptslinelets-2; j > -1; --j)          
        {
          ipt_u1dold[j] = (ipt_u1dold[j] - matp[j]*ipt_u1dold[j+1])/mat[j];          
        }

        L2norm = 0.0;

        E_Int j_cvg;

        for (E_Int j = 1; j < nbptslinelets-1; ++j)
        {
          if (L2norm<=K_FUNC::E_abs((ipt_u1dold[j] - u1d[j])/u1d[j]))
          {
            L2norm = K_FUNC::E_abs((ipt_u1dold[j] - u1d[j])/u1d[j]);
            j_cvg = j;
          }
          u1d[j] = ipt_u1dold[j];                 
        }

        if (npass == 2)
        {
          L2norm0 = L2norm;
        }
        else if (npass < 2)
        {
          L2norm = 1.0;
          L2norm0=1.0;
        }