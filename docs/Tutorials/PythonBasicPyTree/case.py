# case.cgns for tutorial
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, 'Density=1.')
C.convertPyTree2File(a, 'case.cgns')
