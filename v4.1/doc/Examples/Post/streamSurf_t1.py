# - streamSurf (array) -
import Converter as C
import Post as P
import Generator as G
import Geom as D
import math as M
import KCore.test as test

ni = 30; nj = 40; nk = 5

def F(x): return M.cos(x)

# Maillage en noeuds
# 3D
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,nk))
m = [m1,m2]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)

b = D.line( (0.1,5.,0.1), (0.1,5.,3.9), N=5 )
b = C.convertArray2Tetra(b)

# 3D struct
C._addVars(b,'rou')
C._addVars(b,'rov')
C._addVars(b,'row')

p = P.streamSurf(m, b,['rou','rov','row'])
test.testA([p], 1)

# 3D non struct TETRA
m2 = C.convertArray2Tetra(m)
p = P.streamSurf(m2,b,['rou','rov','row'])
test.testA([p], 2)

# 4e test : mixte
p = P.streamSurf([m[0],m2[1]], b,['rou','rov','row'])
test.testA([p], 3)


# 3D non struct HEXA
m2 = C.convertArray2Hexa(m)
p = P.streamSurf(m2,b,['rou','rov','row'])
test.testA([p], 2)

# 3D non struct NGon
m2 = C.convertArray2NGon(m)
p = P.streamSurf(m2,b,['rou','rov','row'])
test.testA([p], 4)
