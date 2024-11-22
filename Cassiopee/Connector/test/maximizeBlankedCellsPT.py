# - maximizeBlankedCells (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

def F(x,y):
    if (x+y<1): return 1
    else: return 2

Ni = 50; Nj = 50
a = G.cart((0,0,0),(1./(Ni-1),1./(Nj-1),1),(Ni,Nj,1))
a = C.initVars(a,'cellN', F,
               ['CoordinateX','CoordinateY'])
a = C.node2Center(a, 'cellN')
a = C.rmVars(a,'cellN')
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t = X.maximizeBlankedCells(t, 2)
C.convertPyTree2File(t, 'out.cgns')
