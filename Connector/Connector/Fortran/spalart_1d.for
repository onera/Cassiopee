C  calcul du champ de viscosité turbulente sur le maillage TBLE
c     par le modèle de Spalart-Allmaras simplifié sous les hypothèses
c     TBLE et convection negligée
c
      subroutine spalart_1d(ithread,y,matm,mat,matp,nutilde,utble,
     &                        pdtc,nu,nutildeext,jmax,kappa)

      implicit none
C       include 'dom.h'
C       include 'temps.h'

      INTEGER_E jmax,ithread
      REAL_E nu,nutildeext,kappa
      REAL_E y(jmax)
      REAL_E utble(jmax)
C       REAL_E mata(jmax,jmax)
      REAL_E matm(1:jmax),mat(1:jmax),matp(1:jmax)
      REAL_E nutilde(jmax)
      REAL_E pdtc,m

      INTEGER_E j,npass,nmax,indstatio
      REAL_E nutdn(1:jmax)
      REAL_E dy,dym,dyp,ynm,ynp,nutildem,nutildep,Pv,Dv,Mv1,Mv2,dt
      REAL_E Fo,Focible,L2norm,L2norm0
      REAL_E A,B,Aprime,Bprime,P,D,Pprime,Dprime,Source
      REAL_E M1m,M1d,M1p,M2m,M2d,M2p
      REAL_E fw,g,r,fv2,fv1,xi,vorticity,St,Sb
      REAL_E fwprime,gprime,rprime,fv2prime,fv1prime,Stprime
      REAL_E sigma,Cb2,Cw1,Cb1,Cw3,Cw2,Cv1
      REAL_E inpout(1:jmax)
      REAL_E rlim
      logical statio

      statio  = .true.

      rlim =10.
      if (statio) then
        indstatio = 1
c       Pas de sous-itérations
        nmax = 1
c       nmax = 20
      else
        indstatio = 0
c       Nombre de sous-itérations maximum
        nmax = 20
c       nmax = 1
      endif
c     write(6,*) nmax

c     Pseudo pas de temps si constant
c     dt = 1.d-8
c     write(6,*) pdtc,y(2)**2/(nu)

c     Fourier cible si pas de temps local
      Focible = 10.
c     Focible = 10.
c
      Cb1 = 0.1355
      Cb2 = 0.622
      sigma = 2./3.
      Cw1 = Cb1/kappa**2 + (1. + Cb2)/sigma
      Cw2 = 0.3
      Cw3 = 2.
      Cv1 = 7.1
c
c     Copie nutilde itération temporelle précédente
C       do j=1,jmax
C         nutdn(j) = nutilde(j)
C       enddo
c
c     CL nutilde = 0 à la paroi
      nutilde(1) = 0.
c     CL externe : match avec nutilde maillage externe
      nutilde(jmax) = nutildeext
c
      npass = 0
      L2norm = 1.
      L2norm0 = 1.
c     Boucle sous-itérations
      DO WHILE ((L2norm/L2norm0.ge.1.d-2).and.(npass.lt.nmax))

        npass = npass + 1

        matm =0.0
        mat  =0.0
        matp =0.0
        mat(1)    = 1.0
        mat(jmax) = 1.0

C         mata(1:,1:) = 0.
c       CL delta_nutilde = 0 à la paroi
C         mata(1,1) = 1.
        inpout(1) = 0.
c       CL externe delta_nutilde = 0
C         mata(jmax,jmax) = 1.
        inpout(jmax) = 0.

      do j=2,jmax-1
c       Cell Vertex
c       calcul coordonnées des noeuds de la cellule
c       et des coeffs pour les fluxs
        ynm = 0.5*(y(j-1) + y(j))
        ynp = 0.5*(y(j) + y(j+1))
        dy = ynp - ynm
        dym = y(j) - y(j-1)
        dyp = y(j+1) - y(j)
c       nutilde aux interfaces
        nutildem = 0.5*(nutilde(j-1) + nutilde(j))
        nutildep = 0.5*(nutilde(j) + nutilde(j+1))
c       vorticité
        vorticity = abs((utble(j+1) - utble(j-1))/(2.*dy))
c       fonctions du modèle
        xi = nutilde(j)/nu
        fv1 = xi**3/(xi**3 + Cv1**3)
        fv2 = 1. - xi/(1. + xi*fv1)
        Sb = nutilde(j)/(kappa**2*y(j)**2)*fv2
        St = vorticity + Sb
        r = min(nutilde(j)/(St*kappa**2*y(j)**2), rlim)
        g = r + Cw2*(r**6 - r)
c  MT : On limite aussi g sinon ca claque
        g=min(g,10000.)
        fw = g*((1. + Cw3**6)/(g**6 + Cw3**6))**(1./6.)
c
c       Fourier local pour pas de temps local
        dt = Focible*dy**2/(nu + nutilde(j))
c       write(6,*) dt,pdtc
c
c       Calcul du second membre RHS
c       fluxs diffusifs, schéma centré ordre 2
c       terme en divergence de gradient
        Mv1 = (1. + Cb2)/sigma*1./dy*
     &        ((nu + nutildep)*(nutilde(j+1) - nutilde(j))/dyp
     &      - (nu + nutildem)*(nutilde(j) - nutilde(j-1))/dym)
c       terme en laplacien
        Mv2 = - Cb2/sigma*(nu + nutilde(j))*1./dy*
     &          ((nutilde(j+1) - nutilde(j))/dyp
     &        - (nutilde(j) - nutilde(j-1))/dym)
c       terme de production
        Pv = Cb1*St*nutilde(j)
c       terme de destruction
        Dv = Cw1*fw*(nutilde(j)/y(j))**2
c       RHS : somme
C         inpout(j) = - (1 - indstatio)*(nutilde(j) - nutdn(j))/pdtc
C      &              + Mv1 + Mv2 + Pv - Dv

        inpout(j) = Mv1 + Mv2 + Pv - Dv


c       Calcul des jacobiennes des termes implicites
c       fluxs diffusifs, schéma centré ordre 2
c       terme en divergence de gradient :
        M1m =   (1. + Cb2)/sigma*(nu + nutildem)/(dy*dym)
        M1d = - (1. + Cb2)/sigma*(nu + nutildem)/(dy*dym)
     &        - (1. + Cb2)/sigma*(nu + nutildep)/(dy*dyp)
        M1p =   (1. + Cb2)/sigma*(nu + nutildep)/(dy*dyp)
c       terme en laplacien
        M2m = - Cb2/sigma*(nu + nutilde(j))/(dy*dym)
        M2d =   Cb2/sigma*(nu + nutilde(j))/(dy*dym)
     &        + Cb2/sigma*(nu + nutilde(j))/(dy*dyp)
        M2p = - Cb2/sigma*(nu + nutilde(j))/(dy*dyp)

c       termes source
        fv1prime = (3.*nutilde(j)**2/nu**3*(xi**3 + Cv1**3)
     &             - xi**3*3.*nutilde(j)**2/nu**3)/(xi**3 + Cv1**3)**2
        fv2prime = (nutilde(j)/nu*(fv1/nu + nutilde(j)/nu*fv1prime)
     &             - 1./nu*(1. + nutilde(j)/nu*fv1))
     &             /(1. + nutilde(j)/nu*fv1)**2
        Stprime = 1./(kappa**2*y(j)**2)*(fv2 + nutilde(j)*fv2prime)
        P = Cb1*St
        Pprime = Cb1*Stprime
        if (r.ge.rlim) then
          rprime = 0.
        else
          rprime = (St - nutilde(j)*Stprime)/(St**2*kappa**2*y(j)**2
     &                                        +1.e-10)
        endif
        gprime = (1. + Cw2*(6.*r**5 - 1.))*rprime
        A = ((1. + Cw3**6)/(g**6 + Cw3**6))**(1./6.)
        B = (1. + Cw3**6)/(g**6 + Cw3**6)
        Bprime = -(1. + Cw3**6)*gprime*6*g**5/(g**6 + Cw3**6)**2
        Aprime = Bprime*1./6.*B**(-5./6.)
        fwprime = gprime*A + g*Aprime
        D = Cw1/y(j)**2*fw*nutilde(j)
        Dprime = Cw1/y(j)**2*(fwprime*nutilde(j) + fw)
c
        source = max(D - P, 0.) + max(Dprime - Pprime, 0.)*nutilde(j)

c       Remplissage matrice
C         mata(j,j-1) = - (M1m + M2m)
C         mata(j,j)   =   1./dt + (1 - indstatio)*1./pdtc
C      &                - (M1d + M2d - source)
C         mata(j,j+1) = - (M1p + M2p)

        matm(j) = - (M1m + M2m)
        mat(j)   =   1./dt + (1 - indstatio)*1./pdtc
     &                - (M1d + M2d - source)
        matp(j) = - (M1p + M2p)


      enddo


      do j=2,jmax 
          m      = matm(j)/mat(j-1)
          mat(j) = mat(j) - m*matp(j-1)
          inpout(j) = inpout(j) - m*inpout(j-1)
      enddo

        inpout(jmax) = inpout(jmax)/mat(jmax)

      
        do j=jmax-1,1,-1
          inpout(j) = (inpout(j) - matp(j)*inpout(j+1))/mat(j)         
        enddo
c     Resolution du systeme tridiag : calcul de nutilde
c     Méthode de Thomas
C       call tdma(jmax,mata(1:jmax,1:jmax),inpout(1:jmax))

      L2norm = 0.
      do j=2,jmax-1
        L2norm = L2norm + (inpout(j)/(nu + nutilde(j)))**2
        nutilde(j) = nutilde(j) + inpout(j)
      enddo
      L2norm = sqrt(L2norm)

C       if (npass .eq. 2)then        
C           L2norm0 = L2norm      
C       elseif (npass .le. 2)then
        
C           L2norm = 1.0
C           L2norm0=1.0
C       endif 

      ENDDO

C       if(ithread .eq. 1)then
C         if(npass.ge.nmax)write(*,*)"SA doesn't converge"
C C         write(*,*)"SA convergence",L2norm/L2norm0,npass
C       endif 


      end subroutine