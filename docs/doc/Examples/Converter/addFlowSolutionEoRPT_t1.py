# - addFlowSolution (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10))
C._addBC2Zone(a,"wall","BCWall","imin")
C._initVars(a,'centers:Density',1.)
t = C.newPyTree(['Base', a])
C.addState2Node__(t[2][1],'GoverningEquations','Euler')
notest = 1
for gvEq in ["Euler","NSLaminar","NSTurbulent",None]:
    tp = elsAProfile.addFlowSolutionEoR(t,governingEquations=gvEq)
    test.testT(tp,notest)
    notest+=1

for extractBC in [False,True]:
    tp = elsAProfile.addFlowSolutionEoR(t,addBCExtract=extractBC)
    test.testT(tp,notest)
    notest+=1

for frame in ['relative','absolute']:
    tp = elsAProfile.addFlowSolutionEoR(t,writingFrame=frame)
    test.testT(tp,notest)
    notest+=1

for protocol in ['end','iteration','after']:
    tp = elsAProfile.addFlowSolutionEoR(t,protocol=protocol)
    test.testT(tp,notest)
    notest+=1
