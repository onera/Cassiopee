# - deform (array) -
import Converter as C
import Transform as T
import KCore.test as test

def F(a):
    a = C.initVars(a, 'dx', 1)
    a = C.initVars(a, 'dy', 0.5)
    a = C.initVars(a, 'dz', 0.5)
    a = T.deform(a, vector=['dx', 'dy', 'dz'])
    a = C.rmVars(a,['dx','dy','dz'])
    return a
test.stdTestA(F)
