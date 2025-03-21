# - streamRibbon (array) -
import Converter as C
import Post as P
import Generator as G
import math as M
import KCore.test as test

ni = 30; nj = 40

def F(x): return M.cos(x)
n = (0.,0.2,0.)
pt = (0.55,0.5,0.)
#pt = (0.1,5,0) # ne marche pas
pt = (5,5,0) # ne marche pas
vect = ['rou','rov','row']
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
p = P.streamRibbon(m, pt, n, vect)
test.testA([p], 1)

# 3D non struct
m2 = C.convertArray2Tetra(m)
p = P.streamRibbon(m2, pt, n, vect)
test.testA([p], 2)

# Mixte
p = P.streamRibbon([m[0],m2[1]],pt, n, vect)
test.testA([p], 3)
