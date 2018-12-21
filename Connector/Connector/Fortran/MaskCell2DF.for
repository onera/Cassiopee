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
C xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4
C of cell (i,j,k) in 2D
C==============================================================================
            ip = MIN( i+1, nic )
            jp = MIN( j+1, nj-1 )
C     
C     Cell (i,j)
C
            ind = i + j*ni  
            xp1 = meshX(ind)
            yp1 = meshY(ind)
            xmincell = xp1
            ymincell = yp1
            xmaxcell = xp1
            ymaxcell = yp1
C
C     Cell (i+1,j)
C
            ind  = ip + j*ni
            xp2 = meshX(ind)
            yp2 = meshY(ind)
            xmincell = MIN( xmincell,xp2)
            ymincell = MIN( ymincell,yp2)
            xmaxcell = MAX( xmaxcell,xp2)
            ymaxcell = MAX( ymaxcell,yp2)               
C     
C     Cell (i+1, j+1)
C     
            ind  = ip + jp*ni       
            xp3 = meshX(ind)
            yp3 = meshY(ind)
            xmincell = MIN( xmincell,xp3)
            ymincell = MIN( ymincell,yp3)
            xmaxcell = MAX( xmaxcell,xp3)
            ymaxcell = MAX( ymaxcell,yp3)  
C
C     Cell (i,j+1)
C     
            ind  = i + jp*ni      
            xp4 = meshX(ind)
            yp4 = meshY(ind)
            xmincell = MIN( xmincell,xp4)
            ymincell = MIN( ymincell,yp4)
            xmaxcell = MAX( xmaxcell,xp4)
            ymaxcell = MAX( ymaxcell,yp4)           
