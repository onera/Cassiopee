# - distrib1 (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C

a = D.line((0,0,0), (4,4,0), N=30)
b = D.distrib1(a, 0.01)
c = G.map(a, b)
C.convertPyTree2File(c, 'out.cgns')
