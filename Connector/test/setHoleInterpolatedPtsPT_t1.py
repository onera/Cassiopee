# - setHoleInterpolatedPts (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test
def sphere(x,y,z):
    if x*x+y*y+z*z < 0.5**2 : return 0.
    else: return 1.

# Cas structure
# Champ cellN en noeud
a = G.cart((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
b = T.translate(a,(2,0,0)); b[0] = 'cart2'
t = C.newPyTree(['Cart',a,b])
t = X.connectMatch(t)
t = C.fillEmptyBCWith(t,'nref','BCFarfield')
C._initVars(t,'Density',1.)
C._initVars(t,'cellN', sphere, ['CoordinateX','CoordinateY','CoordinateZ'])
nod = 1
for d in [-2,-1,0,1,2,5]:
    t2 = X.setHoleInterpolatedPoints(t,depth=d,loc='nodes')
    test.testT(t2,nod); nod+=1
#
# Champ en centres
#
a = G.cart((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
b = T.translate(a,(2,0,0)); b[0] = 'cart2'
t = C.newPyTree(['Cart',a,b])
t = X.connectMatch(t)
t = C.fillEmptyBCWith(t,'nref','BCFarfield')
C._initVars(t,'Density',1.)
C._initVars(t,'cellN', sphere, ['CoordinateX','CoordinateY','CoordinateZ'])
t = C.node2Center(t,'cellN'); t = C.rmVars(t,'cellN')
for d in [-2,-1,0,1,2,5]:
    t2 = X.setHoleInterpolatedPoints(t,depth=d,loc='centers')
    test.testT(t2,nod); nod+=1
