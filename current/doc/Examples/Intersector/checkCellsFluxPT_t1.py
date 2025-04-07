# - selectCells (pyTree) -
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Converter.Internal as I
import KCore.test as test


t = C.convertFile2PyTree('boolNG_M1.tp')
t = C.convertArray2NGon(t)
t = XOR.closeCells(t)

I._createElsaHybrid(t, method=1, methodPE=1)

XOR.checkCellsFlux(t)

test.testT(t,1)
