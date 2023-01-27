# - boolean checkCellsClosure (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Intersector.PyTree as XOR

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

err = XOR.checkCellsClosure(M1)