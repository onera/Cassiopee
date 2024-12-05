# - setDoublyDefinedBC (array) 3D -
import Converter as C
import Connector as X
import Generator as G
import KCore.test as test

a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30))
celln = C.array('cellN',b[2]-1,b[3]-1,b[4]-1)
celln = C.initVars(celln, 'cellN', 2)

cellna = C.array('cellN',a[2]-1,a[3]-1,a[4]-1)
cellna = C.initVars(cellna, 'cellN', 1)
range = [1,a[2],1,a[3],1,1]
range = [1,1,1,a[3],1,a[4]]

celln = X.setDoublyDefinedBC(b, celln, [a], [cellna], range, depth=2)
bc = C.node2Center(b)
bc = C.addVars([bc,celln])

test.testA([bc])
