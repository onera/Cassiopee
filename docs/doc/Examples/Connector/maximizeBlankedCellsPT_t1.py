# - maximizeBlankedCells (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

def F(x,y):
    if (x+y<1): return 1
    else: return 2

# CAS 2D
Ni = 50; Nj = 50
a = G.cart((0,0,0),(1./(Ni-1),1./(Nj-1),1),(Ni,Nj,1))
C._initVars(a,'cellN', F,['CoordinateX','CoordinateY'])
a = C.node2Center(a, 'cellN')
C._rmVars(a, 'cellN')
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension',2)
C._addVars(t,'Density')
t2 = X.maximizeBlankedCells(t, depth=1,dir=0)
test.testT(t2,1)
t2 = X.maximizeBlankedCells(t, depth=2,dir=0)
test.testT(t2,2)
t2 = X.maximizeBlankedCells(t, depth=2, dir=1)
test.testT(t2,21)

# CAS 3D
Ni = 50; Nj = 50; Nk = 20
a = G.cart((0,0,0),(1./(Ni-1),1./(Nj-1),1./(Nk-1)),(Ni,Nj,Nk))
a = C.initVars(a,'cellN', F,['CoordinateX','CoordinateY'])
a = C.node2Center(a, 'cellN')
a = C.rmVars(a,'cellN')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension',3)
t = C.addVars(t,'Density')
t2 = X.maximizeBlankedCells(t, depth=1,dir=0)
test.testT(t2,3)
t2 = X.maximizeBlankedCells(t, depth=2,dir=0)
test.testT(t2,4)
