import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import Transform.PyTree as T
import KCore.test as test

NBBASES = 2
INTERSECT = 1

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,11,12))
b = G.cart((1.5,0.,0.), (0.1,0.1,0.2), (15,15,8))
a = C.addBC2Zone(a,'overlap','BCOverlap','imax')
b = C.addBC2Zone(b,'overlap','BCOverlap','imin')
b = C.addBC2Zone(b,'overlap','BCOverlap','kmax')

# 1 seule base, pas d'extrapolation
t = C.newPyTree(['Base']); t[2][1][2].append(a); t[2][1][2].append(b)
t = X.applyBCOverlaps(t, depth=2)
t = X.setInterpolations(t, loc='cell',storage='direct', sameBase=1)
test.testT(t,1)

# 1 seule base, extrapolations
a2 = T.rotate(a,(1.5,0.,0.), (0.,1.,0.), 15.)
b2 = T.rotate(b,(1.5,0.,0.), (0.,1.,0.), -15.)
t = C.newPyTree(['Base']); t[2][1][2].append(a2); t[2][1][2].append(b2)
t = X.applyBCOverlaps(t, depth=2)
t = X.setInterpolations(t, loc='cell', storage='direct', sameBase=1)
test.testT(t,2)

# 2 bases differentes, pas d'extrapolation
t = C.newPyTree(['Base1','Base2']); t[2][1][2].append(a); t[2][2][2].append(b)
t = X.applyBCOverlaps(t, depth=2)
t = X.setInterpolations(t, loc='cell', storage='direct', sameBase=1)
test.testT(t,3)

# 2 bases differentes, extrapolations
t = C.newPyTree(['Base1','Base2']); t[2][1][2].append(a2); t[2][2][2].append(b2)
t = X.applyBCOverlaps(t, depth=2)
t = X.setInterpolations(t, loc='cell', storage='direct', sameBase=1)
test.testT(t,4)
