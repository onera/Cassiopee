# - streamSurf (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import Geom.PyTree as D
import math as M
import KCore.test as test

ni = 30; nj = 40; nk = 5
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk)); m1[0] = 'cart1'
m2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,nk)); m2[0] = 'cart2'
b = D.line( (0.1,5.,0.1), (0.1,5.,3.9), N=5 )
b = C.convertArray2Tetra(b)

def F(x): return M.cos(x)
m = [m1,m2]
m = C.initVars(m, 'u', 1.)
m = C.initVars(m, 'v', F, ['CoordinateX'])
m = C.initVars(m, 'w', 0.)
m = C.initVars(m,'centers:G',1.)
m = C.addBC2Zone(m,'wall','BCWall','imin')
m = C.addBC2Zone(m,'BCOverlap','BCOverlap','imax')
# 3D struct
t = C.newPyTree(['Base']); t[2][1][2] += m
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
p = P.streamSurf(t, b,['u','v','w'])
t = C.addBase2PyTree(t,'Base2',2); t[2][2][2] = [p]
test.testT(t, 1)
