# - closeCells (array) -
import Intersector as XOR
import Converter as C

m = C.convertFile2Arrays('boolNG_M1.tp')
m = C.convertArray2NGon(m[0])

m = XOR.closeCells(m)
C.convertArrays2File([m], 'out.plt')