import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0),(1,1,1),(11,11,1))
C._initVars(a,"cellN=(2.*({CoordinateX}<0.5)+({CoordinateX}>0.5))")
res =  X.getInterpolatedPoints(a,loc='nodes')
test.testO(res,1)
#

a = G.cart((0,0,0),(1,1,1),(11,11,1))
C._initVars(a,"centers:cellN=(2.*({centers:CoordinateX}<1.)+({centers:CoordinateX}>1.))")
res =  X.getInterpolatedPoints(a,loc='centers')
test.testO(res,2)
