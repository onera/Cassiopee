# - getDistantIndex (array) -
import Geom as D

a = D.naca(12., 5001)
l = D.getLength(a)
print("Longueur totale du profil = %f."%l)
print("Indice du point eloigne du bord de %f = %d."%(l/10., D.getDistantIndex(a, 1, l/10.)))
print("Indice du point eloigne du bord de %f = %d."%(-l/10., D.getDistantIndex(a, 4701, -l/10.)))
