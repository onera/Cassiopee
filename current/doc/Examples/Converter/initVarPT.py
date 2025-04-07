# - initVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))

# Init from a function
def F(x1, x2): return 3.*x1+2.*x2
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.initVars(a, 'centers:F', F, ['centers:CoordinateX','centers:CoordinateY'])
C.convertPyTree2File(a, 'out.cgns')
