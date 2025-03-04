C  
C    Copyright 2013-2025 Onera.
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

C Common part to blanking applied to TETRA zones for a given TETRA element et
C     
         ind = cn1(et)-1
         xp1 = xt(ind)
         yp1 = yt(ind)
         zp1 = zt(ind)
         xmincell = xp1
         ymincell = yp1
         zmincell = zp1
         xmaxcell = xp1
         ymaxcell = yp1
         zmaxcell = zp1
C     
         ind = cn2(et)-1
         xp2 = xt(ind)
         yp2 = yt(ind)
         zp2 = zt(ind)
         xmincell = MIN(xmincell, xp2)
         ymincell = MIN(ymincell, yp2)
         zmincell = MIN(zmincell, zp2)
         xmaxcell = MAX(xmaxcell, xp2)
         ymaxcell = MAX(ymaxcell, yp2)
         zmaxcell = MAX(zmaxcell, zp2)        
C
         ind = cn3(et)-1
         xp3 = xt(ind)
         yp3 = yt(ind)
         zp3 = zt(ind)
         xmincell = MIN(xmincell, xp3)
         ymincell = MIN(ymincell, yp3)
         zmincell = MIN(zmincell, zp3)
         xmaxcell = MAX(xmaxcell, xp3)
         ymaxcell = MAX(ymaxcell, yp3)
         zmaxcell = MAX(zmaxcell, zp3)   
C
         ind = cn4(et)-1
         xp4 = xt(ind)
         yp4 = yt(ind)
         zp4 = zt(ind)
         xmincell = MIN(xmincell, xp4)
         ymincell = MIN(ymincell, yp4)
         zmincell = MIN(zmincell, zp4)
         xmaxcell = MAX(xmaxcell, xp4)
         ymaxcell = MAX(ymaxcell, yp4)
         zmaxcell = MAX(zmaxcell, zp4)    
