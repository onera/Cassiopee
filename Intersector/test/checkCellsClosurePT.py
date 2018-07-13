# - boolean checkCellsClosure (array) -
import Generator.PyTree as G
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Geom.PyTree as D

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

err = XOR.checkCellsClosure(M1)



