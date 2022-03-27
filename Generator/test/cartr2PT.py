# - cartr2 (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cartr2((10,5,1), (1,1,1), (1.5,1.3,1.), (200.,100.,100.))
C.convertPyTree2File(a, 'out.cgns')
