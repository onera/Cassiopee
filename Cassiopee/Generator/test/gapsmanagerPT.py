# - gapsmanager (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Generator.PyTree as G

a = D.sphere6((0,0,0), 1, N=10)
a = C.node2Center(a)
a = C.convertArray2Tetra(a)
b = G.gapsmanager(a, mode=2)
C.convertPyTree2File(b, 'out.cgns')
