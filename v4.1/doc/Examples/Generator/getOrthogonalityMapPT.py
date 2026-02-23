# - getOrthogonalityMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,10))
a = G.getOrthogonalityMap(a)
C.convertPyTree2File(a, 'out.cgns')
