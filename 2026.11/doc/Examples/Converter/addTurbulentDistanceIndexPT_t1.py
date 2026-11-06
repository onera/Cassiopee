# - addTurbulentDistanceIndex (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.elsAProfile as CE
import KCore.test as test
import Converter.Internal as Internal

Internal.__FlowSolutionCenters__='FlowSolution#Init'
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
t = C.newPyTree(['Base',a])
C._initVars(t,'centers:TurbulentDistance', 1.)
tp = CE.addTurbulentDistanceIndex(t)
test.testT(tp,1)
#in place
CE._addTurbulentDistanceIndex(t)
test.testT(t,2)
