# - integ (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

# 2D
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m, 'vx', 1.)
m = C.initVars(m, 'centers:vy', 2.)
m = C.initVars(m, 'ratio', 1.)
res = P.integ(m,'vx') + P.integ(m,'centers:vy')
test.testO(res,1)

# 1D
ni = 30
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,1,1))
m = C.initVars(m, 'vx', 1.)
m = C.initVars(m, 'centers:vy', 2.)
m = C.initVars(m, 'ratio', 1.)
res = P.integ(m,'vx') + P.integ(m,'centers:vy')
test.testO(res,2)

# TRI
ni = 30; nj = 40
m = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m, 'vx', 1.)
m = C.initVars(m, 'centers:vy', 2.)
m = C.initVars(m, 'ratio', 1.)
res = P.integ(m,'vx') + P.integ(m,'centers:vy')
test.testO(res,3)

# BAR
ni = 30
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,1,1))
m = C.convertArray2Tetra(m)
m = C.initVars(m, 'vx', 1.)
m = C.initVars(m, 'centers:vy', 2.)
res = P.integ(m,'vx') + P.integ(m,'centers:vy')
test.testO(res,4)
