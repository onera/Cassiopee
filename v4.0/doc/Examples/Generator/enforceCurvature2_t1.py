# test - enforceCurvature2 (array)
import Converter as C
import Geom as D
import Generator as G
import KCore.test as test

# Cas 1 : courbure constante
a = D.circle((0,0,0),0.1,N=51)
ni = 101; distrib = G.cart((0,0,0),(1./(ni-1),1,1),(ni,1,1))
distrib2 = G.enforceCurvature2(distrib,a)
test.testA([distrib2],1)

# Cas 2 : courbure variable
a = D.naca(12.)
ni = 101; distrib = G.cart((0,0,0),(1./(ni-1),1,1),(ni,1,1))
distrib2 = G.enforceCurvature2(distrib,a)
test.testA([distrib2],2)
