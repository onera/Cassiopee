# - setHoleInterpolatedPts (array) -
import Converter as C
import Connector as X
import Generator as G
import Post as P
import KCore.test as test
def sphere(x,y,z):
    if x*x+y*y+z*z < 0.48**2: return 0.
    else: return 1.

# Cas structure: champ cellN en noeud
a = G.cart((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
a = C.initVars(a,'cellN', sphere, ['x','y','z'])
nod = 1
for d in [-2,-1,0,1,2,5]:
    celln = X.setHoleInterpolatedPoints(a, depth=d)
    test.testA([celln],nod); nod+=1
#
# Champ en centres
#
a = G.cart((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
ac = C.node2Center(a)
ac = C.initVars(ac,'cellN', sphere, ['x','y','z'])
for d in [-2,-1,0,1,2,5]:
    celln = X.setHoleInterpolatedPoints(ac, depth=d)
    test.testA([celln],nod); nod+=1


# Méthode "octahedron" - dir 3
# Cas structure: champ cellN en noeud
def cube(x,y,z):
    if 0. <= x < 0.05 and 0. <= y < 0.05 and 0. <= z < 0.05: return 0.
    else: return 1.

depth = 5
a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
a = C.initVars(a,'cellN', cube, ['x','y','z'])
cellN = X.setHoleInterpolatedPoints(a, depth=depth, dir=3)
#t1 = P.selectCells(cellN, '{cellN}==2')
test.testA([cellN], nod)
