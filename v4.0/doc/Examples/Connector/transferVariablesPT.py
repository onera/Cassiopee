# - transferVariables (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

a = G.cart((0.,0.,0.),(1,1,1),(11,11,11))
b = G.cart((9.,0.,0.),(0.5,0.5,1),(22,10,11))
c = G.cart((9.,3.,0.),(1,1,1),(11,10,11))
C._initVars(a,"centers:cellN=2*({centers:CoordinateX}>9)+({centers:CoordinateX}<9.)")
indicesI,XI,YI,ZI = X.getInterpolatedPoints(a,loc='centers')
C._initVars(a,"centers:Density=1.")
C._initVars(b,"centers:Density=20.")
C._initVars(c,"centers:Density=100.")
C._initVars(b,"centers:cellN=1.")
C._initVars(c,"centers:cellN=1.")
for zdnr in [b,c]:
    zdnrc = C.node2Center(zdnr)
    adt = C.createHook(zdnrc,'adt')
    fields = X.transferFields(zdnrc, XI, YI, ZI, hook=adt, variables=['Density'])
    print(fields[1][-1,:])
    C.freeHook(adt)
