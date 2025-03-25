# - streamLine (array) -
import Converter as C
import Post as P
import Generator as G
import math as M
import KCore.test as test

ni = 30; nj = 40
def F(x): return M.cos(x)

# Maillage en noeuds
# 3D
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,2))
m = [m1,m2]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)

# 2D
s1 = G.cart((0,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,1))
s2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,1))
s = [s1,s2]
s = C.initVars(s, 'rou', 1.)
s = C.initVars(s, 'rov', F, ['x'])
s = C.initVars(s, 'row', 0.)

# 3D struct
x0 = 0.1; y0 = 5.; z0 = 0.
p = P.streamLine(m, (x0,y0,z0),['rou','rov','row'])
test.testA([p], 1)

# 2D struct
x0 = 5.; y0 = 5.; z0 = 0.
p = P.streamLine(s, (x0,y0,z0),['rou','rov','row'], N=200)
test.testA([p], 2)

# 3D non struct
m2 = C.convertArray2Tetra(m)
p = P.streamLine(m2,(x0,y0,z0),['rou','rov','row'])
test.testA([p], 3)

# 2D non struct
s2 = C.convertArray2Tetra(s)
p = P.streamLine(s2, (x0,y0,z0),['rou','rov','row'])
test.testA([p], 4)

# 4e test : mixte
p = P.streamLine([m[0],m2[1]],(x0,y0,z0),['rou','rov','row'])
test.testA([p], 5)
