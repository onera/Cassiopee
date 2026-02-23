# - cellN2OversetHoles (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

a = G.cart((0,0,0),(1,1,1),(10,10,10))
C._initVars(a, 'centers:cellN', 0)
a = X.cellN2OversetHoles(a)
C.convertPyTree2File(a, 'out.cgns')
