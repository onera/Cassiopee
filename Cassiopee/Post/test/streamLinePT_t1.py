# - streamLine (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import math as M
import KCore.test as test

ni = 30; nj = 40
def F(x): return M.cos(x)

# Maillage en noeuds
# 3D
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,2))
m = [m1,m2]
m = C.initVars(m, 'u', 1.)
m = C.initVars(m, 'v', F, ['CoordinateX'])
m = C.initVars(m, 'w', 0.)
m = C.initVars(m,'centers:G',1.)
m = C.addBC2Zone(m,'wall','BCWall','imin')
m = C.addBC2Zone(m,'BCOverlap','BCOverlap','imax')
# 3D struct
x0 = 0.1; y0 = 5.; z0 = 0.
t = C.newPyTree(['Base']); t[2][1][2] += m
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
p = P.streamLine(t, (x0,y0,z0),['u','v','w'])
t = C.addBase2PyTree(t,'Base2',1); t[2][2][2] = [p]
test.testT(t, 1)
