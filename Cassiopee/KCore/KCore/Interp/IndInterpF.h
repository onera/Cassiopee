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

c***********************************************************************
c
c  G.Desquesnes
c
c  ACT
c_A  Fonction d indexation des points d interpolation
c
c  IN
c_I  npts_1d : nombre de points d interpolation en 1D 
c_I  i, j, k : coordonnees du point d interpolation dans l espace de
c_I            reference
c
c  OUT
c_O  ind_interp : id de l emplacement memoire du point specifie
c
c  I/O
c_/
c
c***********************************************************************

c 
c ---- declaration de la fonction
c 
      INTEGER_E ind_interp

c 
c ---- declaration des variables de la fonction
c 
      INTEGER_E ind_interp_npts_1d__
      INTEGER_E ind_interp_i__, ind_interp_j__, ind_interp_k__

c 
c ---- definition de la fonction
c 
        ind_interp( ind_interp_npts_1d__, 
     &              ind_interp_i__, ind_interp_j__, ind_interp_k__) = 
     &  ind_interp_i__ + 
     & (ind_interp_j__-1)*ind_interp_npts_1d__ + 
     & (ind_interp_k__-1)*ind_interp_npts_1d__*ind_interp_npts_1d__

C===================== IndInterpF.h ===============================
