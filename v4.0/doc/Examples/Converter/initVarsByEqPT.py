# - initVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, '{Density} = 3 * {CoordinateX} + sin({CoordinateY})')
a = C.initVars(a, '{centers:MomentumX} = 3 * {centers:CoordinateX} + sin({centers:CoordinateY})')
C.convertPyTree2File(a, 'out.cgns')
