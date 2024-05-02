# - merge (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Connector.PyTree as X
import Geom.PyTree as D
import KCore.test as test
def f(t,u):
    x = t+u
    y = t*t+1+u*u
    z = u
    return (x,y,z)

# surface grid
a = D.surface(f, N=50)
b = T.splitSize(a, 100)
b = X.connectMatch(b, dim=2)
t = C.newPyTree(['Surface', 2]); t[2][1][2] += b
C._initVars(t, 'F', 1.); C._initVars(t, 'centers:G', 2.)
C._addBC2Zone(t[2][1][2][0],'overlap','BCOverlap','imin')
C._fillEmptyBCWith(t,'wall','BCWall',dim=2)
b = T.merge(t)
t2 = C.newPyTree(['Surface', 2, b])
test.testT(t2,1)
b = T.merge(t, alphaRef=45.)
t2 = C.newPyTree(['Surface', 2, b])
test.testT(t2,2)
#
# cas d'angles vifs
#
a1 = G.cart((0,0,0), (1,1,1), (11,11,1))
a3 = G.cart((10,0,0), (1,1,1), (11,11,1))
a2 = T.rotate(a1, (0,0,0), (1,0,0), 90.)
a2 = T.reorder(a2, (-1,2,3)); a2[0] = C.getZoneName(a1[0])
t = C.newPyTree(['Base',2]); t[2][1][2] += [a1,a2,a3]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
C._addBC2Zone(t[2][1][2][0],'overlap','BCOverlap','imin')
t = X.connectMatch(t, dim=2)
C._fillEmptyBCWith(t, 'wall', 'BCWall', dim=2)
a = T.merge(t, alphaRef=45.)
t = C.newPyTree(['Base', 2]); t[2][1][2] += a
test.testT(t,3)

# volume grid
a = D.surface(f, N=50)
distrib = G.cart((0.,0.,0.),(0.1,1,1),(11,1,1))
a = G.addNormalLayers(a, distrib)
b = T.splitSize(a, 5000)
b = X.connectMatch(b, dim=3)
t = C.newPyTree(['Base']); t[2][1][2] += b
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
C._addBC2Zone(t[2][1][2][0],'overlap','BCOverlap','imin')
C._fillEmptyBCWith(t,'wall','BCWall',dim=3)
b = T.merge(t, dir=2)
t = C.newPyTree(['Base']); t[2][1][2] += b
test.testT(t,4)
