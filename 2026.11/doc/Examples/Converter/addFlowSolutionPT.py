# - addFlowSolution (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a])
t = elsAProfile.addFlowSolution(t, governingEquations='Euler')
import Converter.Internal as Internal
Internal.printTree(t)
C.convertPyTree2File(t, 'out.cgns')
