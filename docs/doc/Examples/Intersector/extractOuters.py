# - extractOuterLayers (array) -
import Converter as C
import Intersector as XOR
import Generator as G

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])
M1 = C.conformizeNGon(M1); M1 = G.close(M1)
m = XOR.extractOuterLayers(M1, 1, discard_external=0)

C.convertArrays2File(m, "out.plt")
