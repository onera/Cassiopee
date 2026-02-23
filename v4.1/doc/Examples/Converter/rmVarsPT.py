# - rmVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = C.addVars(a, ['MomentumX','centers:Density'])
b = C.rmVars(b, ['MomentumX','centers:Density'])
C.convertPyTree2File(b, 'out.cgns')
