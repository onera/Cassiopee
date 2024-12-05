# - extractVars (array) -
import Generator as G
import Converter as C

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'F', 2.)
# Var defined by a string
r = C.extractVars(a, 'F')
# Vars defined by a list
r = C.extractVars(a, ['x','y'])
C.convertArrays2File(r, 'out.plt')
