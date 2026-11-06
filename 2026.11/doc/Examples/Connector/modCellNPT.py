# - modCellN (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X

a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._initVars(a, 'centers:cellN=1.')

C.setValue(a, 'centers:cellN', 22, 0.)
C.setValue(a, 'centers:cellN', 24, 2.)

X._modCellN1(a)
C.convertPyTree2File(a, 'out.cgns')
