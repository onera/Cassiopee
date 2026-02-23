# - rmGhostCellsNGon (array) -
import Converter as C
import Generator as G

a = G.cartNGon((0,0,0),(1,1,1),(21,21,1))
a = G.close(a)
b = C.addGhostCellsNGon(a, depth=2)
b = C.rmGhostCellsNGon(b, depth=2)
C.convertArrays2File(b, "out.plt")
