# - splitSize (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

t = C.newPyTree(['Base'])
a = G.cart((0,0,0),(1,1,1),(51,21,11)); a[0] = 'cart1'
a2 = T.translate(a,(46,0,0)); a2[0] = 'cart2'
a3 = T.translate(a,(-50,0,0)); a3[0] = 'cart3'
a = C.addBC2Zone(a,'overlap1','BCOverlap','imax')
a2 = C.addBC2Zone(a2,'overlap2','BCOverlap','imin')
t[2][1][2] += [a,a2,a3]
t = C.fillEmptyBCWith(t,'wall','BCWall')
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.splitSize(t, 700)
t = X.connectMatch(t)
test.testT(t)
