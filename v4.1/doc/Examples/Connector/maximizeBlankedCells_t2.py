# - maximizeBlankedCells (array) 3D -
import Converter as C
import Connector as X
import Generator as G
import KCore.test as test

def F(x,y,z):
    if (x+y<1): return 1
    else: return 2

Ni = 50; Nj = 50; Nk = 50
a = G.cart((0,0,0),(1./(Ni-1),1./(Nj-1),1/(Nk-1)),(Ni,Nj,Nk))
a = C.initVars(a, 'cellN', F, ['x','y','z'])

a1 = X.maximizeBlankedCells(a, 2, dir=0)
test.testA([a1], 1)

a2 = X.maximizeBlankedCells(a, 1, dir=0)
test.testA([a2], 2)

a2 = X.maximizeBlankedCells(a, 1, dir=1)
test.testA([a2], 3)
