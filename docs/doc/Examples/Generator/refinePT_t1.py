# - refine (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test

# cas 2D multibloc avec CL et raccords
a = G.cart( (0,0,0), (0.1,0.1,0.1), (11,11,1))
b = G.cart( (0,1,0), (0.1,0.1,0.1), (11,11,1))
C._addBC2Zone(a,'match1','BCMatch','jmax',b,'jmin')
C._addBC2Zone(a,"wall",'BCWall','imin')
C._addBC2Zone(a,"overlap",'BCOverlap','imax')
C._addBC2Zone(b,'match2','BCMatch','jmin',a,'jmax')
C._addBC2Zone(b,"wall",'BCWall','imin')
C._addBC2Zone(b,"overlap",'BCOverlap','imax')
t = C.newPyTree(['Base',2]); t[2][1][2] += [a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)

# Non integer refinement factor: GCs and BCs are removed
# topTree
factor = 1.5
t2 = G.refine(t, factor, dir=1)
test.testT(t2, 1)

t2 = G.refine(t, factor, dir=0)
test.testT(t2, 2)

# single zone
t2 = Internal.copyRef(t)
for noz in range(2): t2[2][1][2][noz] = G.refine(t2[2][1][2][noz], factor, dir=1)
test.testT(t2, 3)

# Integer refinement factor: GCs and BCs are removed
# topTree
factor = 2.
t2 = G.refine(t, factor, dir=1)
test.testT(t2, 4)

t2 = G.refine(t, factor, dir=0)
test.testT(t2, 5)

# single zone
t2 = Internal.copyRef(t)
for noz in range(2): t2[2][1][2][noz] = G.refine(t2[2][1][2][noz], factor, dir=1)
test.testT(t2, 6)

# Non integer refinement factor: GCs and BCs are removed
# topTree
factor = 0.3
t2 = G.refine(t, factor, dir=1)
test.testT(t2, 7)

t2 = G.refine(t, factor, dir=0)
test.testT(t2, 8)

# single zone
t2 = Internal.copyRef(t)
for noz in range(2): t2[2][1][2][noz] = G.refine(t2[2][1][2][noz], factor, dir=1)
test.testT(t2, 9)

# Integer refinement factor: GCs and BCs are removed
# topTree
factor = 0.5
t2 = G.refine(t, factor, dir=1)
test.testT(t2, 10)

t2 = G.refine(t, factor, dir=0)
test.testT(t2, 11)

# single zone
t2 = Internal.copyRef(t)
for noz in range(2): t2[2][1][2][noz] = G.refine(t2[2][1][2][noz], factor, dir=1)
test.testT(t2, 12)
