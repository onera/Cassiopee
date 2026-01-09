C  
C    Copyright 2013-2026 ONERA.
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

      SUBROUTINE constrconnect(NbofElements, NbofNodes,
     &                         connectEN,
     &                         NbofEdges, NbofEdgesInt, NbofEdgesExt,
     &                         connectEA, connectAE, connectAN,
     &                         connectAN2, connectAoldAnew, 
     &                         dejaVuA)
C
      IMPLICIT NONE
C_IN
      INTEGER_E NbofElements                      ! number of elements
      INTEGER_E NbofNodes                         ! number of nodes
      INTEGER_E connectEN(0:NbofElements-1,3)     ! connectivity Elts -> Nodes
C
C_OUT
      INTEGER_E NbofEdges                         ! number of Edges
      INTEGER_E NbofEdgesInt                      ! number of Interior Edges
      INTEGER_E NbofEdgesExt                      ! number of Exterior Edges
      INTEGER_E connectEA(0:NbofElements-1,3)     ! connectivity Elts ->Edges    --- "A" for "Arete"
      INTEGER_E connectAE(0:3*NbofElements-1,2)   ! connectivity Edges -> Elts   --- unknown size <= 3*NbofElements
      INTEGER_E connectAN(0:3*NbofElements-1,3)   ! connectivity Edges -> Nodes  --- unknown size <= 3*NbofElements
      INTEGER_E dejaVuA(0:3*NbofElements-1)       ! mark on edges = 0,1,2,3,4
C   
C_LOCAL
      INTEGER_E i, j, ind, ind1, ind2, c                                
      INTEGER_E test, tmp, deb                  
      INTEGER_E connectAN2(0:3*NbofElements-1,6)  ! connectivity Edges -> Nodes  --- unknown size <= 3*NbofElements
      INTEGER_E connectAoldAnew(0:3*NbofElements-1) ! connectivity initial number of Edge -> final number of Edge
                                !(after deleting redundant edges and reordering by writing first interior edges)
                                ! unknown size <= 3*NbofElements
      INTEGER_E tmpNbOfEdges
      tmpNbOfEdges = 3* NbofElements
      
C Initialisation des variables
      DO j = 1, 6
         DO i = 0, 3*NbofElements-1
            connectAN2(i, j) = -1
         ENDDO
      ENDDO
      DO i = 0, 3*NbofElements-1
         connectAoldAnew(i) = -1
         dejaVuA(i) = 5
      ENDDO
      DO j = 1, 2
         DO i = 0, 3*NbofElements-1
            connectAE(i, j) = -1
         ENDDO
      ENDDO

C_Construction des tables de connectivite:
C Aretes -> Noeuds (connectAN) 
C Elts   -> Aretes (connectEA)
      DO i = 0, NbofElements-1
         j = 3*i
         connectAN(j,1) = MIN(connectEN(i,1),connectEN(i,2))
         connectAN(j,2) = MAX(connectEN(i,1),connectEN(i,2))
         connectAN(j,3) = j
         connectEA(i,1) = j

         j = 3*i+1
         connectAN(j,1) = MIN(connectEN(i,1),connectEN(i,3))
         connectAN(j,2) = MAX(connectEN(i,1),connectEN(i,3))
         connectAN(j,3) = j
         connectEA(i,2) = j

         j = 3*i+2
         connectAN(j,1) = MIN(connectEN(i,2),connectEN(i,3))
         connectAN(j,2) = MAX(connectEN(i,2),connectEN(i,3))
         connectAN(j,3) = j
         connectEA(i,3) = j
      END DO

C Classement des aretes : tri ameliore
      deb = 0
      DO test = 0, NbofNodes-1
         DO i = deb, tmpNbofEdges-1
            IF ( connectAN(i,1) .EQ. test ) THEN
               tmp = connectAN(deb,1)
               connectAN(deb,1) = connectAN(i,1)
               connectAN(i,1) = tmp

               tmp = connectAN(deb,2)
               connectAN(deb,2) = connectAN(i,2)
               connectAN(i,2) = tmp

               tmp = connectAN(deb,3)
               connectAN(deb,3) = connectAN(i,3)
               connectAN(i,3) = tmp
               deb = deb + 1
            END IF
         END DO
      END DO

C Classement des aretes : tri croissant sur la deuxieme ligne en cas 
C d'egalite sur la premiere
      test = 1
      DO WHILE (test .EQ. 1)
         test = 0
         DO i = 0, tmpNbofEdges-2
            IF ( (connectAN(i,1) .EQ. connectAN(i+1,1) )
     &      .and.(connectAN(i,2) .GT. connectAN(i+1,2) ) ) THEN
               tmp = connectAN(i,2)
               connectAN(i,2) = connectAN(i+1,2)
               connectAN(i+1,2) = tmp

               tmp = connectAN(i,3)
               connectAN(i,3) = connectAN(i+1,3)
               connectAN(i+1,3) = tmp
               test = 1
            END IF
         END DO
      END DO

C Marquage des aretes redondantes et ecriture des aretes definies de facon 
C unique dans connectAN2
      ind = 0
      connectAN2(ind,1) = connectAN(0,1)
      connectAN2(ind,2) = connectAN(0,2)
      connectAN2(ind,3) = connectAN(0,3)

      c = 0
      DO i = 1, tmpNbofEdges-1
         IF ( (connectAN(i,1) .EQ. connectAN(i-1,1) )
     &        .AND.(connectAN(i,2) .EQ. connectAN(i-1,2) ) ) THEN
            IF (c.GE.3) THEN
               WRITE(*,*) 'Error : constrconnect : ',
     &              'an edge is defined more than 3 times'
               STOP
            ENDIF
            connectAN2(ind, 4+c) = connectAN(i, 3)
            c = c + 1
         ELSE
            ind = ind+1
            connectAN2(ind,1) = connectAN(i,1)
            connectAN2(ind,2) = connectAN(i,2)
            connectAN2(ind,3) = connectAN(i,3)
            c = 0
         END IF
      END DO
      NbofEdges = ind+1

C Comptage du nombre d aretes internes et externes
      DO i = 0, NbofEdges-1
         IF (connectAN2(i,4) .NE. -1) THEN
            NbofEdgesInt = NbofEdgesInt+1 ! aretes redondantes
         END IF
      END DO
      NbofEdgesExt = NbofEdges - NbofEdgesInt

C Tri des aretes : rangement des aretes internes au debut et ecriture 
C dans connectAN
C + Correspondance numero initial des aretes -> numero final des aretes : 
C ecriture du tableau conectAoldAnew
C + Marquage des aretes : 0->arete interne ; 3->arete externe ; le marquage 
C sera egal a 1, 2, 4 dans ComputNormF.for
      ind1 = -1
      ind2 = NbofEdgesInt-1
      DO i = 0, NbofEdges-1
         IF (connectAN2(i,4) .NE. -1) THEN ! arete interne
            ind1 = ind1+1
        
            connectAN(ind1,1) = connectAN2(i,1)
            connectAN(ind1,2) = connectAN2(i,2)
            dejaVuA(ind1) = 0
            connectAoldAnew(connectAN2(i,3)) = ind1
            connectAoldAnew(connectAN2(i,4)) = ind1
            IF (connectAN2(i,5) .NE. -1) THEN
               connectAoldAnew(connectAN2(i,5)) = ind1  
            ENDIF
            IF (connectAN2(i,6) .NE. -1) THEN
               connectAoldAnew(connectAN2(i,6)) = ind1  
            ENDIF
         ELSE                   ! arete externe
            ind2 = ind2+1
            connectAN(ind2,1) = connectAN2(i,1)
            connectAN(ind2,2) = connectAN2(i,2)
            dejaVuA(ind2) = 3
            connectAoldAnew(connectAN2(i,3)) = ind2
         END IF
      END DO

C Mise a jour de la table de connectivite Elements -> Aretes avec la 
C nouvelle numerotation des aretes
      DO i = 0, NbofElements-1
         connectEA(i,1) = connectAoldAnew(connectEA(i,1))
         connectEA(i,2) = connectAoldAnew(connectEA(i,2))
         connectEA(i,3) = connectAoldAnew(connectEA(i,3))
      END DO

C_Construction de la table de connectivite Aretes -> Elements : connectAE
      DO i = 0, NbofElements-1
         DO j = 1, 3
            IF (connectAE(connectEA(i,j), 1) .EQ. -1) THEN
               connectAE(connectEA(i,j), 1) = i
            ELSE
               connectAE(connectEA(i,j), 2) = i
            END IF  
         END DO
      END DO

      END
