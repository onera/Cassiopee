# - setHoleInterpolatedPts (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import KCore.test as test
def sphere(x,y,z):
    if x*x+y*y+z*z < 0.48**2 : return 0.
    else: return 1.

# Cas structure
# Champ cellN en noeud
a = G.cart((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
b = T.translate(a,(2,0,0)); b[0] = 'cart2'
t = C.newPyTree(['Cart']); t[2][1][2]+=[a,b]
t = X.connectMatch(t)
t = C.fillEmptyBCWith(t, 'nref', 'BCFarfield')
t = C.initVars(t,'Density', 1.)
t = C.initVars(t,'cellN', sphere, ['CoordinateX','CoordinateY','CoordinateZ'])

tp = Internal.copyTree(t)
t2 = X.setHoleInterpolatedPoints(tp, depth=1, dir=2, loc='nodes')
test.testT(t2,1)

tp = Internal.copyTree(t)
t2 = X.setHoleInterpolatedPoints(tp, depth=2, dir=2, loc='nodes')
test.testT(t2,2)

tp = Internal.copyTree(t)
t2 = X.setHoleInterpolatedPoints(tp, depth=3, dir=2, loc='nodes')
test.testT(t2,3)

#
# Champ en centres
#
a = G.cart((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
b = T.translate(a,(2,0,0)); b[0] = 'cart2'
t = C.newPyTree(['Cart']); t[2][1][2]+=[a,b]
t = X.connectMatch(t)
t = C.fillEmptyBCWith(t,'nref','BCFarfield')
t = C.initVars(t,'Density', 1.)
t = C.initVars(t,'centers:cellN', sphere, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'])

tp = Internal.copyTree(t)
t2 = X.setHoleInterpolatedPoints(tp, depth=1, dir=2, loc='centers')
test.testT(t2,4)

tp = Internal.copyTree(t)
t2 = X.setHoleInterpolatedPoints(tp, depth=2, dir=2, loc='centers')
test.testT(t2,5)

tp = Internal.copyTree(t)
t2 = X.setHoleInterpolatedPoints(tp, depth=3, dir=2, loc='centers')
test.testT(t2,6)
