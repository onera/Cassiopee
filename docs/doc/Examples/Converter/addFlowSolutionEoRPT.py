# - addFlowSolutionEoR (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a])
t = elsAProfile.addFlowSolutionEoR(t, governingEquations='Euler')
C.convertPyTree2File(t, 'out.cgns')
