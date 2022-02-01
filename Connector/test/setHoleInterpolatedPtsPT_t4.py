# - setHoleInterpolatedPts (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test
import Converter.Internal as Internal
def sphere(x,y,z):
    if x*x+y*y+z*z < 0.48**2 : return 0.
    else: return 1.
#
# Champ en noeuds non structure PENTA
#
a = G.cartPenta((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
b = T.translate(a,(2,0,0)); b[0] = 'cart2'
t = C.newPyTree(['Cart']); t[2][1][2]+=[a,b]
t = C.initVars(t,'Density',1.)
t = C.initVars(t,'cellN', sphere, ['CoordinateX','CoordinateY','CoordinateZ'])
nod = 1
for d in [-2,-1,0,1,2,5]:
    tp = Internal.copyTree(t)
    t2 = X.setHoleInterpolatedPoints(t, depth=d)
    test.testT(t2,nod); nod+=1

# Champ cellN en centres
a = G.cartPenta((-2.,-1.,-1.),(0.1,0.1,0.1), (21,21,21))
b = T.translate(a,(2,0,0)); b[0] = 'cart2'
t = C.newPyTree(['Cart']); t[2][1][2]+=[a,b]
t = C.initVars(t,'Density',1.)
t = C.initVars(t,'centers:cellN', sphere, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'])
for d in [-2,-1,0,1,2,5]:
    tp = Internal.copyTree(t)
    t2 = X.setHoleInterpolatedPoints(tp,depth=d)
    test.testT(t2,nod); nod+=1
