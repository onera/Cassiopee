# - initVars (aray) -
import Converter as C
import KCore.test as test

# Create a function
def F(x1, x2): return 3.*x1+2.*x2

test.stdTestA(C.initVars, 'cellN', F, ['x', 'y'])
