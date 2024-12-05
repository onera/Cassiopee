# - getInterpolatedPoints (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

a = G.cart((0,0,0),(1,1,1),(11,11,1))
C._initVars(a,"cellN=(2.*({CoordinateX}<0.5)+({CoordinateX}>0.5))")
res =  X.getInterpolatedPoints(a,loc='nodes')
print(res)
