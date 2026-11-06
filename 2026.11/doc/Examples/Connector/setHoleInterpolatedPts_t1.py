# - setHoleInterpolatedPts (array) -
import Converter as C
import Connector as X
import Generator as G
import KCore.test as test
import numpy

def sphere(x,y,z):
    return (x*x + y*y + z*z >= 0.48**2).astype(float)

# Cas structure: champ cellN en noeud
a = G.cart((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
a = C.initVars(a,'cellN', sphere, ['x','y','z'], isVectorized=True)
nod = 1
for d in [-2,-1,0,1,2,5]:
    celln = X.setHoleInterpolatedPoints(a, depth=d)
    test.testA([celln],nod); nod+=1

# Champ en centres
a = G.cart((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
ac = C.node2Center(a)
ac = C.initVars(ac,'cellN', sphere, ['x','y','z'], isVectorized=True)
for d in [-2,-1,0,1,2,5]:
    celln = X.setHoleInterpolatedPoints(ac, depth=d)
    test.testA([celln],nod); nod+=1

# Méthode "octahedron" - dir 3
# Cas structure: champ cellN en noeud
def cube(x, y, z):
    inside = (0. <= x) & (x < 0.05) & (0. <= y) & (y < 0.05) & (0. <= z) & (z < 0.05)
    return numpy.where(inside, 0., 1.)

depth = 5
a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
a = C.initVars(a,'cellN', cube, ['x','y','z'], isVectorized=True)
cellN = X.setHoleInterpolatedPoints(a, depth=depth, dir=3)
test.testA([cellN], nod)
