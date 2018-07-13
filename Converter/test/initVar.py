# - initVars (array) -
import Converter as C
import Generator as G

# Create a function
def F(x1, x2): return 3.*x1+2.*x2

a = G.cart((0,0,0), (1,1,1), (11,11,1))
a = C.initVars(a, 'F', F, ['x','y'])
C.convertArrays2File([a], "out.plt")
