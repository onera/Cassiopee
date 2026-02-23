# - setHoleInterpolatedPts (array) -
import Converter as C
import Connector as X
import Generator as G
import KCore.test as test
def sphere(x,y,z):
    if x*x+y*y+z*z < 0.48**2 : return 0.
    else: return 1.

# Cas TETRA: champ cellN en noeud
a = G.cartHexa((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
a = C.initVars(a,'cellN', sphere, ['x','y','z'])
nod = 1
for d in [-2,-1,0,1,2,5]:
    celln = X.setHoleInterpolatedPoints(a,depth=d)
    test.testA([celln],nod); nod+=1
#
# Champ en centres
#
a = G.cartHexa((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
ac = C.node2Center(a)
ac = C.initVars(ac,'cellN', sphere, ['x','y','z'])
for d in [-2,-1,0,1,2,5]:
    celln = X.setHoleInterpolatedPoints(ac,depth=d)
    test.testA([celln],nod); nod+=1
