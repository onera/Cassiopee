# - cartr1 (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cartr1((0,0,0), (1.,1.,1.), (1.1,1.2,1.), (10,10,10))
C.convertPyTree2File(a, 'out.cgns')
