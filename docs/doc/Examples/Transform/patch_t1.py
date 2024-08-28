# - patch (array) - ijk indices
import Transform as T
import Generator as G
import Converter as C
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (5,1,1))
a = C.initVars(a, '{F}=1')
test.stdTestA(T.patch,a,(3,1,1))
