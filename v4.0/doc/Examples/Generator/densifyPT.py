# - densify (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

a = D.circle((0,0,0), 1., 10)
b = G.densify(a, 0.01)
C.convertPyTree2File(b, 'out.cgns')
