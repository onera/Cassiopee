# - setInterpData2 (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test
aD = G.cart((0,0,0),(1,1,1), (11,11,11))
aR = G.cart((0,0,0),(0.5,0.5,0.5), (21,21,21))
X._setInterpData2(aR, aD, loc='centers', cartesian=False)
test.testT(aD,1)
