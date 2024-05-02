# - convertArray2NGon(pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

def convertArray2NGon(a):
    a = C.convertArray2NGon(a); a = G.close(a, tol=1e-12)
    return a

test.stdTestT(convertArray2NGon)
