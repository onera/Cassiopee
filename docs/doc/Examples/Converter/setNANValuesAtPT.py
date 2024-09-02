# - setNANValuesAt (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, 'F', 1.)
C._setNANValuesAt(a, 'F', 0.)
C.convertPyTree2File(a, 'out.cgns')
