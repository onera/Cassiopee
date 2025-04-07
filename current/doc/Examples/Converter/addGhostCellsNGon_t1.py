# - addGhostCellsNGon (array) -
import Converter as C
import Generator as G
import Geom as D
import KCore.test as test

a = D.sphere((0,0,0),1)
a = C.convertArray2NGon(a)
a = G.close(a)
b = C.addGhostCellsNGon(a, depth=1)
test.testA([b],1)

a = G.cartNGon((0,0,0),(1,1,1),(21,21,1))
a = C.convertArray2NGon(a)
a = G.close(a)
b = C.addGhostCellsNGon(a, depth=1)
test.testA([b],2)

a = G.cartNGon((0,0,0),(1,1,1),(21,21,1))
a = C.convertArray2NGon(a)
a = G.close(a)
b = C.addGhostCellsNGon(a, depth=2)
test.testA([b],3)

a = G.cartNGon((0,0,0),(1,1,1),(21,21,21))
a = C.convertArray2NGon(a)
a = G.close(a)
b = C.addGhostCellsNGon(a, depth=1)
test.testA([b],4)
