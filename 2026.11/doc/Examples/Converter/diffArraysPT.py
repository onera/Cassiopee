# - diffArrays (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (11,11,11))
a1 = C.initVars(a, "F", 1.); a1 = C.initVars(a1, "centers:Q", 1.2)
a2 = C.initVars(a, "F", 3.); a2 = C.initVars(a2, "centers:Q", 2.)
ret = C.diffArrays(a1, a2)
C.convertPyTree2File(ret, 'out.cgns')
