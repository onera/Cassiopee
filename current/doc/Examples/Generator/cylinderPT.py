# - cylinder (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
C.convertPyTree2File(a, 'out.cgns')
