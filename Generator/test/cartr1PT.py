import Generator.PyTree as G
import Converter.PyTree as C

a = G.cartr1((0,0,0), (1,1,1), (10,10,10),(1.1,1.2,1.3))
print(a)
C.convertPyTree2File(a, 'out.cgns')
