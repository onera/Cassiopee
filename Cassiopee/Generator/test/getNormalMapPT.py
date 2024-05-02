# - getNormalMap (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C

a = D.sphere((0,0,0), 1, 50)
a = G.getNormalMap(a)
C.convertPyTree2File(a, 'out.cgns')
