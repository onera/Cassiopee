# - renameVars (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# a single zone
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
varsN = ['Density','MomentumX']
varsP = ['density','velocity_x']
for v in varsP: m = C.addVars(m, v)
m2 = P.renameVars(m, varsP, varsN)
test.testT([m2],1)

# 2 zones with different variables
ni = 30; nj = 40
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m2 = G.cart((10,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2)); m2[0] = 'cart2'
varsN = ['Density','MomentumX','centers:cellN']
varsP = ['density','velocity_x','centers:ichim']
for v in varsP: m1 = C.addVars(m, v)
C._addVars(m1,'ichim')
C._addVars(m2,'density')
C._initVars(m2,'centers:ichim', 1.)
t = C.newPyTree(['Base',m1,m2])
t = P.renameVars(t, varsP, varsN)
test.testT(t,2)

t = C.newPyTree(['Base',m1,m2])
P._renameVars(t, varsP, varsN)
test.testT(t,3)
