# - addTurbulentDistanceIndex (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.elsAProfile as elsAProfile

Internal.__FlowSolutionCenters__='FlowSolution#Init'
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
t = C.newPyTree(['Base',a])
C._initVars(t, "centers:TurbulentDistance", 1.)
tp = elsAProfile.addTurbulentDistanceIndex(t)
C.convertPyTree2File(tp, 'out.cgns')
