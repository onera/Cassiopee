# - maximizeBlankedCells (array) -
import Converter as C
import Connector as X
import Generator as G

def F(x,y):
    return (x+y>=1).astype(float) + 1

Ni = 50; Nj = 50
a = G.cart((0,0,0),(1./(Ni-1),1./(Nj-1),1),(Ni,Nj,1))
a = C.initVars(a, 'cellN', F, ['x','y'], isVectorized=True)
a = X.maximizeBlankedCells(a, 2)
C.convertArrays2File([a], 'out.plt')
