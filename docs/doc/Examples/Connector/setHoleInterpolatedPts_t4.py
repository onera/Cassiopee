# - setHoleInterpolatedPts (array) -
import Converter as C
import Connector as X
import Generator as G
import KCore.test as test

def sphere(x,y,z):
    return (x*x + y*y + z*z >= 0.48**2).astype(float)

# Cas PENTA: champ cellN en noeud
a = G.cartPenta((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
a = C.initVars(a,'cellN', sphere, ['x','y','z'], isVectorized=True)
nod = 1
for d in [-2,-1,0,1,2,5]:
    celln = X.setHoleInterpolatedPoints(a,depth=d)
    test.testA([celln],nod); nod+=1

# Champ en centres
a = G.cartPenta((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
ac = C.node2Center(a)
ac = C.initVars(ac,'cellN', sphere, ['x','y','z'], isVectorized=True)
for d in [-2,-1,0,1,2,5]:
    celln = X.setHoleInterpolatedPoints(ac,depth=d)
    test.testA([celln],nod); nod+=1
