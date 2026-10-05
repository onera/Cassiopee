# - maximizeBlankedCells (array) -
import Converter as C
import Connector as X
import Generator as G
import KCore.test as test

def F(x,y):
    return (x+y>=1).astype(float) + 1

Ni = 50; Nj = 50
a = G.cart((0,0,0),(1./(Ni-1),1./(Nj-1),1),(Ni,Nj,1))
a = C.addVars(a, 'cellN')
a = C.initVars(a, 'cellN', F, ['x','y'], isVectorized=True)

a1 = X.maximizeBlankedCells(a, 2,dir=0)
test.testA([a1], 1)

a2 = X.maximizeBlankedCells(a, 1,dir=0)
test.testA([a2], 2)
