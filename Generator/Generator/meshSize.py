import math
# Estimation de la taille de la premiere maille a imposer sur le maillage
# INPUT
# La densite non adimensionnee a l'infini
roInf = 1.2
# La densite a la paroi
roWall = roInf
# La vitesse non adimensionne a l'infini
uInf = 68
# Mu a l'infini
muInf = 0.00001711 # celui de l'air
# Longueur d'adimensionnement (longueur du profil)
LRef = 0.05

# Le e/c du profil. Rapport de l'epaisseur sur la corde quand on voit
# l'ecoulement lui arriver dessus
# On pourrait l'estimer par une rotation du corps pour l'avoir a incidence
# nulle + une BB
esurc = 0.012

# Le y+ voulu (2 c'est bien)
yplus = 2.

# le x/L ou on calcule y
xsurL = 0.5

# OUTPUT
# Reynolds inf
ReInf = roWall*uInf*LRef / muInf
print 'ReInf', ReInf

# Correction
Correction = math.exp(4.5*pow(esurc, 1.3))
print 'Correction', Correction

#==============================================================================
# turbulent dimensionne
#==============================================================================
# Frottement
Cf = 0.058*math.pow( ReInf*xsurL, -0.2) * Correction
print 'Cf', Cf

# Vitesse de frottement
utau = math.sqrt(roInf*Cf*uInf*uInf*0.5/roWall)
print 'utau', utau

# DeltaPP
deltapp = 0.16*xsurL*LRef/pow(ReInf*xsurL, 1./7.)
print 'deltapp', deltapp

# hauteur de la premier maille en 
y = muInf * yplus / roInf / utau
print 'y', y

print '------------------'
#==============================================================================
# Laminaire dimensionne
#==============================================================================
# Frottement
Cf = 0.664/math.sqrt(ReInf*xsurL)* Correction
print 'Cf', Cf

# Vitesse de frottement
utau = math.sqrt(roInf*Cf*uInf*uInf*0.5/roWall)
print 'utau', utau

# hauteur de la premier maille en 
y = muInf * yplus / roInf / utau
print 'y', y

print '------------------'
#==============================================================================
# Turbulent adimensionne
#==============================================================================
# Frottement
Cf = 0.058*math.pow( ReInf*xsurL, -0.2) * Correction
print 'Cf', Cf

# Vitesse de frottement
utau = math.sqrt(Cf*0.5)
print 'utau', utau

# DeltaPP
deltapp = 0.16*xsurL/pow(ReInf*xsurL, 1./7.)
print 'deltapp', deltapp

# hauteur de la premier maille en 
y = (1./ReInf) * yplus / utau
print 'y', y

print '------------------'

#==============================================================================
# Laminaire adimensionne
#==============================================================================
# Frottement
Cf = 0.664/math.sqrt(ReInf*xsurL)* Correction
print 'Cf', Cf

# Vitesse de frottement
utau = math.sqrt(Cf*0.5)
print 'utau', utau

# hauteur de la premier maille en 
y = (1./ReInf) * yplus / utau
print 'y', y

print '------------------'




