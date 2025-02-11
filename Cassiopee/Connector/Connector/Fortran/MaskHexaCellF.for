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

C Common part to blanking applied to PENTA zones for a given PENTA element et
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
         xmincell = MIN( xmincell, xp2)
         ymincell = MIN( ymincell, yp2)
         zmincell = MIN( zmincell, zp2)
         xmaxcell = MAX( xmaxcell, xp2)
         ymaxcell = MAX( ymaxcell, yp2)
         zmaxcell = MAX( zmaxcell, zp2)        
C
         ind = cn3(et)-1
         xp3 = xt(ind)
         yp3 = yt(ind)
         zp3 = zt(ind)
         xmincell = MIN( xmincell, xp3)
         ymincell = MIN( ymincell, yp3)
         zmincell = MIN( zmincell, zp3)
         xmaxcell = MAX( xmaxcell, xp3)
         ymaxcell = MAX( ymaxcell, yp3)
         zmaxcell = MAX( zmaxcell, zp3)   
C
         ind = cn4(et)-1
         xp4 = xt(ind)
         yp4 = yt(ind)
         zp4 = zt(ind)
         xmincell = MIN( xmincell, xp4)
         ymincell = MIN( ymincell, yp4)
         zmincell = MIN( zmincell, zp4)
         xmaxcell = MAX( xmaxcell, xp4)
         ymaxcell = MAX( ymaxcell, yp4)
         zmaxcell = MAX( zmaxcell, zp4)    
C
         ind = cn5(et)-1
         xp5 = xt(ind)
         yp5 = yt(ind)
         zp5 = zt(ind)
         xmincell = MIN( xmincell, xp5)
         ymincell = MIN( ymincell, yp5)
         zmincell = MIN( zmincell, zp5)
         xmaxcell = MAX( xmaxcell, xp5)
         ymaxcell = MAX( ymaxcell, yp5)
         zmaxcell = MAX( zmaxcell, zp5)    
C
         ind = cn6(et)-1
         xp6 = xt(ind)
         yp6 = yt(ind)
         zp6 = zt(ind)
         xmincell = MIN( xmincell, xp6)
         ymincell = MIN( ymincell, yp6)
         zmincell = MIN( zmincell, zp6)
         xmaxcell = MAX( xmaxcell, xp6)
         ymaxcell = MAX( ymaxcell, yp6)
         zmaxcell = MAX( zmaxcell, zp6)    
C
         ind = cn7(et)-1
         xp7 = xt(ind)
         yp7 = yt(ind)
         zp7 = zt(ind)
         xmincell = MIN( xmincell, xp7)
         ymincell = MIN( ymincell, yp7)
         zmincell = MIN( zmincell, zp7)
         xmaxcell = MAX( xmaxcell, xp7)
         ymaxcell = MAX( ymaxcell, yp7)
         zmaxcell = MAX( zmaxcell, zp7)    
C
         ind = cn8(et)-1
         xp8 = xt(ind)
         yp8 = yt(ind)
         zp8 = zt(ind)
         xmincell = MIN( xmincell, xp8)
         ymincell = MIN( ymincell, yp8)
         zmincell = MIN( zmincell, zp8)
         xmaxcell = MAX( xmaxcell, xp8)
         ymaxcell = MAX( ymaxcell, yp8)
         zmaxcell = MAX( zmaxcell, zp8)    
