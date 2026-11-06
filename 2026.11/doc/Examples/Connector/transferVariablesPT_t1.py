import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test
import Geom.PyTree as D

# transfers at centers
a = D.sphere((5.,5.,5.),1)
C._initVars(a,"centers:cellN=2.")
C._initVars(a,"centers:Density=0.")
indicesI,XI,YI,ZI = X.getInterpolatedPoints(a,loc='centers')

zdnr = G.cart((0.,0.,0.),(1,1,1),(11,11,11))
C._initVars(zdnr,"centers:Density=1.")
C._initVars(zdnr,"centers:cellN=1.")
zdnrc = C.node2Center(zdnr)
adt = C.createHook(zdnrc,'adt')
fields = X.transferFields(zdnrc, XI, YI, ZI, hook=adt, variables=['Density'])
C.freeHook(adt)
test.testA(fields,1)

# transfers at nodes
a = D.sphere((5.,5.,5.),1)
C._initVars(a,"cellN=2.")
C._initVars(a,"Density=0.")
indicesI,XI,YI,ZI = X.getInterpolatedPoints(a,loc='nodes')

zdnr = G.cart((0.,0.,0.),(1,1,1),(11,11,11))
C._initVars(zdnr,"Density=1.")
C._initVars(zdnr,"cellN=1.")
adt = C.createHook(zdnr,'adt')
fields = X.transferFields(zdnr, XI, YI, ZI, hook=adt, variables=['Density'])
C.freeHook(adt)
test.testA(fields,2)

# transfers with orphan or extrapolated points

a = D.sphere((0.5,0.5,0.5),0.1)
C._initVars(a,"cellN=2.")
C._initVars(a,"Density=0.")
indicesI,XI,YI,ZI = X.getInterpolatedPoints(a,loc='nodes')

Np = 101; hp = 1./(Np-1)
zdnr = G.cart((0.,0.,0.),(hp,hp,hp),(Np,Np,Np))
C._initVars(zdnr,"Density=1.")
C._initVars(zdnr,"cellN=2*({CoordinateX}>0.4999)+({CoordinateX}>0.4999)")
adt = C.createHook(zdnr,'adt')
fields = X.transferFields(zdnr, XI, YI, ZI, hook=adt, variables=['Density'])
C.freeHook(adt)
#print Converter.getMinValue(fields,'donorVol'), Converter.getMaxValue(fields,'donorVol')
test.testA(fields,3)
