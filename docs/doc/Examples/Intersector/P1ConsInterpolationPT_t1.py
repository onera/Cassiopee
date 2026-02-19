# - P1ConservativeInterpolation (pyTree) -
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Generator.PyTree as G
import KCore.test as test

N = 10
aD = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (N,N,N))
aD = C.convertArray2NGon(aD)
aD = C.initVars(aD, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')
#C.convertPyTree2File(aD, 'aD.cgns')

aR = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (N+1,N+1,N+1))
aR = C.convertArray2NGon(aR)

aR = XOR.P1ConservativeInterpolation(aR, aD)
#C.convertPyTree2File(aR, 'aR.cgns')

test.testT(aR,1)
