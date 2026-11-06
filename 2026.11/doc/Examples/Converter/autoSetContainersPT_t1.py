# - autoSetContainers (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
Internal.__FlowSolutionCenters__ = 'FlowSolution#EndOfRun'
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
Internal.__FlowSolutionCenters__ = 'FlowSolution#Centers'
t = C.newPyTree(['Base',a])
Internal.autoSetContainers(t)
test.testO([Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__])
