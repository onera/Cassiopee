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

C==============================================================================
C Compute xmincell, ymincell, zmincell, xmaxcell, ymaxcell, zmaxcell
C xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4,...
C of cell (i,j,k) in 3D
C==============================================================================
               ip = MIN( i+1, nic )
               jp = MIN( j+1, nj-1 )
               kp = MIN( k+1, nk-1 )
C     
C     Cell (i,j,k)
C     
               ind = i + j*ni + k*nij
               xp1 = meshX(ind)
               yp1 = meshY(ind)
               zp1 = meshZ(ind)
               xmincell = xp1
               ymincell = yp1
               zmincell = zp1
               xmaxcell = xp1
               ymaxcell = yp1
               zmaxcell = zp1
               
C     Cell (i+1,j,k)
C
               ind  = ip + j*ni + k*nij
               xp2 = meshX(ind)
               yp2 = meshY(ind)
               zp2 = meshZ(ind)
               xmincell = MIN( xmincell, xp2)
               ymincell = MIN( ymincell, yp2)
               zmincell = MIN( zmincell, zp2)
               xmaxcell = MAX( xmaxcell, xp2)
               ymaxcell = MAX( ymaxcell, yp2)
               zmaxcell = MAX( zmaxcell, zp2)              
C     
C     Cell (i+1, j+1, k)
C     
               ind  = ip + jp*ni + k*nij                 
               xp3 = meshX(ind)
               yp3 = meshY(ind)
               zp3 = meshZ(ind)
               xmincell = MIN( xmincell, xp3)
               ymincell = MIN( ymincell, yp3)
               zmincell = MIN( zmincell, zp3)
               xmaxcell = MAX( xmaxcell, xp3)
               ymaxcell = MAX( ymaxcell, yp3)
               zmaxcell = MAX( zmaxcell, zp3)              
C
C     Cell (i,j+1, k)
C     
               ind  = i + jp*ni + k* nij                  
               xp4 = meshX(ind)
               yp4 = meshY(ind)
               zp4 = meshZ(ind)
               xmincell = MIN( xmincell, xp4)
               ymincell = MIN( ymincell, yp4)
               zmincell = MIN( zmincell, zp4)
               xmaxcell = MAX( xmaxcell, xp4)
               ymaxcell = MAX( ymaxcell, yp4)
               zmaxcell = MAX( zmaxcell, zp4)               
C     
C     Cell (i,j,k+1)
C     
               ind = i + j*ni + kp*nij
               xp5 = meshX(ind)
               yp5 = meshY(ind)
               zp5 = meshZ(ind)
               xmincell = MIN( xmincell, xp5)
               ymincell = MIN( ymincell, yp5)
               zmincell = MIN( zmincell, zp5)
               xmaxcell = MAX( xmaxcell, xp5)
               ymaxcell = MAX( ymaxcell, yp5)
               zmaxcell = MAX( zmaxcell, zp5) 
C
C     Cell (i+1,j,k+1)
C     
               ind  = ip + j*ni + kp*nij
               xp6 = meshX(ind)
               yp6 = meshY(ind)
               zp6 = meshZ(ind)
               xmincell = MIN( xmincell, xp6)
               ymincell = MIN( ymincell, yp6)
               zmincell = MIN( zmincell, zp6)
               xmaxcell = MAX( xmaxcell, xp6)
               ymaxcell = MAX( ymaxcell, yp6)
               zmaxcell = MAX( zmaxcell, zp6) 
C     
C     Cell (i+1, j+1, k+1)
C     
               ind  = ip + jp*ni + kp*nij    
               xp7 = meshX(ind)
               yp7 = meshY(ind)
               zp7 = meshZ(ind)           
               xmincell = MIN( xmincell, xp7)
               ymincell = MIN( ymincell, yp7)
               zmincell = MIN( zmincell, zp7)
               xmaxcell = MAX( xmaxcell, xp7)
               ymaxcell = MAX( ymaxcell, yp7)
               zmaxcell = MAX( zmaxcell, zp7)             
C     
C     Cell (i,j+1,k+1)
C     
               ind  = i + jp*ni + kp*nij   
               xp8 = meshX(ind)
               yp8 = meshY(ind)
               zp8 = meshZ(ind)            
               xmincell = MIN( xmincell, xp8)
               ymincell = MIN( ymincell, yp8)
               zmincell = MIN( zmincell, zp8)
               xmaxcell = MAX( xmaxcell, xp8)
               ymaxcell = MAX( ymaxcell, yp8)
               zmaxcell = MAX( zmaxcell, zp8) 
