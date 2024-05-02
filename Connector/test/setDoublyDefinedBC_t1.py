# - setDoublyDefinedBC (array) 3D -
import Converter as C
import Connector as X
import Generator as G
import KCore.test as test

a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30))
celln = C.array('cellN',a[2]-1,a[3]-1,a[4]-1)
celln = C.initVars(celln, 'cellN', 1)
indmax = celln[2]*celln[3]
celln[1][0][0:indmax] = 2
cellnb = C.array('cellN',b[2]-1,b[3]-1,b[4]-1)
cellnb = C.initVars(cellnb, 'cellN', 1)
celln = X.setDoublyDefinedBC(a, celln, [b], [cellnb], [1,a[2],1,a[3],1,1])
ac = C.node2Center(a)
ac = C.addVars([ac,celln])
test.testA([ac],1)
