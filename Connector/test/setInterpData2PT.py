# - setInterpData2 (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X

aD = G.cart((0,0,0),(1,1,1), (11,11,11))
aR = G.cart((0,0,0),(0.5,0.5,0.5), (21,21,21))
X._setInterpData2(aR, aD, loc='centers', cartesian=False)
C.convertPyTree2File(aD, "out.cgns")
