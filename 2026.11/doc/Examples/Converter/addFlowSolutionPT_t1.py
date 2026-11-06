# - addFlowSolution (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
C._addBC2Zone(a,"wall","BCWall","imin")
C._initVars(a,'centers:Density',1.)

t = C.newPyTree(['Base', a])
C.addState2Node__(t[2][1],'GoverningEquations','Euler')
VARS_CONS=["Density","MomentumX","MomentumY","MomentumZ","EnergyStagnationDensity"]
VARS_SA = ["TurbulentSANuTildeDensity"]
outputDict={}
outputDict["tata"]=0
outputDict["toto"]='machaine'
notest = 1

for name in ["","#OUT"]:
    tp = elsAProfile.addFlowSolution(t, name=name, variables=[])
    test.testT(tp, notest)
    notest += 1
for loc in ["Vertex","CellCenter","cellfict"]:
    tp = elsAProfile.addFlowSolution(t, loc=loc, variables=[])
    test.testT(tp, notest)
    notest += 1
for variables in [None, [],VARS_CONS,VARS_CONS+VARS_SA]:
    tp = elsAProfile.addFlowSolution(t, variables=variables)
    test.testT(tp, notest)
    notest += 1
for gvEq in [None,"Euler","NSLaminar","NSTurbulent"]:
    tp = elsAProfile.addFlowSolution(t, governingEquations=gvEq)
    test.testT(tp, notest)
    notest += 1
for wmode in [None,0,1]:
    tp = elsAProfile.addFlowSolution(t,writingMode=wmode,variables=[])
    test.testT(tp, notest)
    notest += 1

for frame in ['relative','absolute']:
    tp = elsAProfile.addFlowSolution(t,writingFrame=frame,variables=[])
    test.testT(tp, notest)
    notest += 1

for period in [None,10]:
    tp = elsAProfile.addFlowSolution(t,period=period,variables=[])
    test.testT(tp, notest)
    notest += 1

for output in [None,outputDict]:
    tp = elsAProfile.addFlowSolution(t,output=output,variables=[])
    test.testT(tp, notest)
    notest += 1

for extractBC in [False,True]:
    tp = elsAProfile.addFlowSolution(t,addBCExtract=extractBC,variables=[])
    test.testT(tp, notest)
    notest += 1

for protocol in ['end','iteration','after']:
    tp = elsAProfile.addFlowSolution(t,protocol=protocol,variables=[])
    test.testT(tp, notest)
    notest += 1
