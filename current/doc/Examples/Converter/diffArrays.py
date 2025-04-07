# - diffArrays (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (11,11,11))
a = C.initVars(a, "F", 1.)
a = C.initVars(a, "Q", 1.2)

b = G.cart((0,0,0), (1,1,1), (11,11,11))
b = C.initVars(b, "Q", 2.)
b = C.initVars(b, "F", 3.)

ret = C.diffArrays([a], [b])
C.convertArrays2File(ret, 'out.plt')
