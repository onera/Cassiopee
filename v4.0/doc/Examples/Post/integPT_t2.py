# - integ (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

# 2D
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.addBC2Zone(m,'overlap','BCOverlap','imin')
m = C.initVars(m, 'vx', 1.); m = C.initVars(m, 'centers:vy', 1.)
m = C.initVars(m, 'ratio', 1.)
t = C.newPyTree(['Base',2]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m)
res = P.integ(t,'vx') + P.integ(t,'centers:vy')
test.testO(res,1)

# 1D
ni = 30
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,1,1))
m = C.initVars(m, 'vx', 1.)
m = C.initVars(m, 'centers:vy', 1.)
m = C.initVars(m, 'ratio', 1.)
m = C.addBC2Zone(m,'overlap','BCOverlap','imin')
m = C.initVars(m, 'vx', 1.); m = C.initVars(m, 'centers:vy', 1.)
m = C.initVars(m, 'ratio', 1.)
t = C.newPyTree(['Base',1]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m)
res = P.integ(t,'vx') + P.integ(t,'centers:vy')
test.testO(res,2)

# TRI
ni = 30; nj = 40
m = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m, 'vx', 1.)
m = C.initVars(m, 'centers:vy', 1.)
m = C.initVars(m, 'ratio', 1.)
t = C.newPyTree(['Base',2]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m)
res = P.integ(t,'vx') + P.integ(t,'centers:vy')
test.testO(res,3)

# BAR
ni = 30
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,1,1))
m = C.convertArray2Tetra(m)
m = C.initVars(m, 'vx', 1.)
m = C.initVars(m, 'centers:vy', 1.)
m = C.initVars(m, 'ratio', 1.)
t = C.newPyTree(['Base',1]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m)
res = P.integ(t,'vx') + P.integ(t,'centers:vy')
test.testO(res,4)
