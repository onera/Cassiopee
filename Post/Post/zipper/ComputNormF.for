C  
C    Copyright 2013-2019 Onera.
C
C    This file is part of Cassiopee.
C
C    Cassiopee is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    Cassiopee is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.

      SUBROUTINE computnorm(NbofNodes, coord,
     &                      NbofElements, connectEA,
     &                      NbofEdges, connectAE,
     &                      refelmt, refpts,
     &                      dejaVuA, norm,
     &                      connectEN2, dejaVuE)
C
      IMPLICIT NONE
C

C
C_IN
      INTEGER_E NbofNodes                        ! number of nodes
      REAL_E coord(0:NbofNodes-1,3)              ! nodes coordinates
      INTEGER_E NbofElements                     ! number of elements
      INTEGER_E connectEA(0:NbofElements-1,3)    ! connectivity Elements -> Edges
      INTEGER_E NbofEdges                        ! number of edges
      INTEGER_E connectAE(0:NbofEdges-1,2)       ! connectivity Edges -> Elements
      INTEGER_E refelmt                          ! element de reference determine dans Integ.C a partir des noeuds de reference
      INTEGER_E refpts(0:2)                      ! numeros ORDONNES des noeuds de l element de reference
C_OUT
      INTEGER_E dejaVuA(0:NbofEdges-1)           ! mark on edges : 0->interior ; 3->exterior (initially)
      REAL_E norm(0:NbofElements-1,3)            ! coordinates of the normales
      INTEGER_E connectEN2(0:NbofElements-1,3)   ! connectivity Elements -> Nodes ; built in order not to modify connectEN
      INTEGER_E dejaVuE(0:NbofElements-1)        ! mark on elements : 0->not treated ; 1->treated
C_LOCAL
      INTEGER_E i                                ! compteur
      INTEGER_E j                                ! compteur
      REAL_E u(0:2)                              ! vecteur 0 du triangle etudie
      REAL_E v(0:2)                              ! vecteur 1 du triangle etudie
      REAL_E w(0:2)                              ! w = u vectoriel v
      REAL_E wl                                  ! longueur de w
      INTEGER_E e(0:1)                           ! 
      INTEGER_E ind(0:1)                         ! 
      INTEGER_E test                             ! 
      INTEGER_E eref                             ! 
      INTEGER_E enew                             ! 
      INTEGER_E ptref(0:2)                       ! 
      INTEGER_E ptnew(0:2)                       ! 
      INTEGER_E tmp                              ! 
      REAL_E coordloc(0:2,3)                     ! coordinates of the triangle studied
      INTEGER_E c(0:2)                           !

C To make the process stop if it fails 
      INTEGER_E  maxloop, loop
C to avoid an infinite loop in case of norms of opposite sense anyway
      maxloop = NbofElements
      loop = 0 
       
C Determination des points 0,1,2 de chaque element de proche en proche a 
C partir d elements deja traites (dejaVuE = 1)
C Identification des noeuds de l'element deja vu avec l'element a traiter
C Inversion de 2 points si necessaire (cad si les noeuds sont ranges dans 
C l ordre oppose)

C Initialisations de tableaux
      DO i = 0, nbOfElements-1 
         dejaVuE(i) = 0
      ENDDO
C
C On initialise avec l'element de reference
      connectEN2(refelmt, 1) = refpts(0) 
      connectEN2(refelmt, 2) = refpts(1)
      connectEN2(refelmt, 3) = refpts(2)
C On marque que l element de reference a ete traite
      dejaVuE(refelmt) = 1      

      c(0) = connectEA(refelmt, 1)
      c(1) = connectEA(refelmt, 2)
      c(2) = connectEA(refelmt, 3)

C On incremente le marquage des aretes traitees
      dejaVuA(c(0)) = dejaVuA(c(0)) + 1          
      dejaVuA(c(1)) = dejaVuA(c(1)) + 1
      dejaVuA(c(2)) = dejaVuA(c(2)) + 1

C Compte le nombre d'elements traites, ie dont on a ordonne les noeuds
      test = 1          
C     Calcul de proche en proche
      
      DO WHILE (test .LT. NbofElements .AND.  loop .LE. maxloop)
         loop = loop + 1

         DO i = 0, NbofEdges-1
            j = dejaVuA(i)
C Recherche d'une arete interne traitee une fois et une seule
            IF (j .EQ. 1 .AND. connectAE(i,1).GE.0 
     &           .AND. connectAE(i,2) .GE. 0) THEN  

               e(0) = connectAE(i,1)
               e(1) = connectAE(i,2)
               ind(0) = dejaVuE(e(0))
               ind(1) = dejaVuE(e(1))
C determination de l'element qui a deja ete traite : element de reference
               IF (ind(1) .EQ. 1 ) THEN        
                  eref = e(1)
                  enew = e(0)
               ELSE IF (ind(0) .EQ. 1 ) THEN 
                  eref = e(0)
                  enew = e(1)
               ELSE
C                  WRITE(*,*) 'ERROR :
C     &                 dejaVu must be equal to 1 for one element'
C                  WRITE(*,*) ' STOP !'
C                  STOP
                  GOTO 1
               END IF
               ptref(0) = connectEN2(eref,1)
               ptref(1) = connectEN2(eref,2)
               ptref(2) = connectEN2(eref,3)
               ptnew(0) = connectEN2(enew,1)
               ptnew(1) = connectEN2(enew,2)
               ptnew(2) = connectEN2(enew,3)

               IF (ptnew(0) .EQ. ptref(0)) THEN
                  IF (ptnew(1) .eq. ptref(1)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ELSE IF (ptnew(2) .eq. ptref(2)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ENDIF
               ELSE IF (ptnew(0) .eq. ptref(1)) THEN
                  IF (ptnew(1) .eq. ptref(2)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ELSE IF (ptnew(2) .eq. ptref(0)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ENDIF
               ELSE IF (ptnew(0) .eq. ptref(2)) THEN
                  IF (ptnew(1) .eq. ptref(0)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ELSE IF (ptnew(2) .eq. ptref(1)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ENDIF
               ELSE IF (ptnew(1) .eq. ptref(0)) THEN
                  IF (ptnew(2) .eq. ptref(1)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ENDIF
               ELSE IF (ptnew(1) .eq. ptref(1)) THEN
                  IF (ptnew(2) .eq. ptref(2)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ENDIF
               ELSE IF (ptnew(1) .eq. ptref(2)) THEN
                  IF (ptnew(2) .eq. ptref(0)) THEN
                     tmp = ptnew(0)
                     ptnew(0) = ptnew(1)
                     ptnew(1) = tmp
                  ENDIF
               ENDIF

               connectEN2(enew,1) = ptnew(0)
               connectEN2(enew,2) = ptnew(1)
               connectEN2(enew,3) = ptnew(2)
C marquage a 1 de l'element nouvellement traite
               dejaVuE(enew) = 1 

               c(0) = connectEA(enew,1)
               c(1) = connectEA(enew,2)
               c(2) = connectEA(enew,3)
               
               dejaVuA(c(0)) = dejaVuA(c(0))+ 1
               dejaVuA(c(1)) = dejaVuA(c(1))+ 1
               dejaVuA(c(2)) = dejaVuA(c(2))+ 1
               test = test+1
            ENDIF
 1          CONTINUE
         END DO
      END DO
      IF ( loop .GE. maxloop) THEN
         WRITE(*,*) 
     &        'WARNING : computenorm : CANNOT DIRECT NORMALS',
     &        ' IN THE SAME SENSE'
         WRITE(*,*) 'Number of elements not treated : ',
     &        NbofElements-test
         
      ENDIF 
C Calcul de la normale de chaque element
      DO i = 0, NbofElements-1
         c(0) = connectEN2(i,1) ! connect(i,1) = noeud 1 de l'element i
         c(1) = connectEN2(i,2) ! coord(A,1) = coordonnee x ("1") du noeud A
         c(2) = connectEN2(i,3)

         coordloc(0,1) = coord(c(0),1)              ! x du point 0
         coordloc(0,2) = coord(c(0),2)              ! y du point 0
         coordloc(0,3) = coord(c(0),3)              ! z du point 0
         coordloc(1,1) = coord(c(1),1)              ! x du point 1
         coordloc(1,2) = coord(c(1),2)              ! y du point 1
         coordloc(1,3) = coord(c(1),3)              ! z du point 1
         coordloc(2,1) = coord(c(2),1)              ! x du point 2
         coordloc(2,2) = coord(c(2),2)              ! y du point 2
         coordloc(2,3) = coord(c(2),3)              ! z du point 2

         u(0) = coordloc(1,1)-coordloc(2,1)         ! x du vecteur 0
         u(1) = coordloc(1,2)-coordloc(2,2)         ! y du vecteur 0
         u(2) = coordloc(1,3)-coordloc(2,3)         ! z du vecteur 0
         v(0) = coordloc(0,1)-coordloc(2,1)         ! x du vecteur 1
         v(1) = coordloc(0,2)-coordloc(2,2)         ! y du vecteur 1
         v(2) = coordloc(0,3)-coordloc(2,3)         ! z du vecteur 1
         w(0) = u(1)*v(2) - u(2)*v(1)
         w(1) = u(2)*v(0) - u(0)*v(2)
         w(2) = u(0)*v(1) - u(1)*v(0)
         wl = sqrt(w(0)*w(0) + w(1)*w(1) + w(2)*w(2))
         wl = MAX(wl, 1.e-12)
         norm(i,1) = w(0)/wl
         norm(i,2) = w(1)/wl
         norm(i,3) = w(2)/wl
      ENDDO

      END
