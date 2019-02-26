# - autoSetContainers (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
Internal.__FlowSolutionCenters__ = 'FlowSolution#EndOfRun'
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={CoordinateY}')
Internal.__FlowSolutionCenters__ = 'FlowSolution#Centers'
t = C.newPyTree(['Base',a])
Internal.autoSetContainers(t)
print(Internal.__FlowSolutionNodes__)
print(Internal.__FlowSolutionCenters__)
