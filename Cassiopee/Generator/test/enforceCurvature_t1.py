# - enforceCurvature (array) -
import Geom as D
import Generator as G
import Transform as T
import KCore.test as test

# Naca profile with lines
a = D.naca(12., 500)
l1 = D.getLength(a)
a2 = D.line((1.,0.,0.),(2.,0.,0.), 500)
l2 = D.getLength(a2)
b = T.join(a, a2)
c = D.line((2.,0.,0.),(1.,0.,0.), 500)
res = T.join(c, b)

# Distribution on the profile
l = l1+2*l2
Ni = 100; Nj = 100
p1 = l2/l; p2 = (l2+l1)/l
h = (p2-p1)/(Ni-1)
distrib = G.cart((p1,0,0), (h, 0.25/Nj,1), (Ni,Nj,1))
distrib = G.enforceCurvature(distrib, res, 0.6)
test.testA([distrib], 1)
