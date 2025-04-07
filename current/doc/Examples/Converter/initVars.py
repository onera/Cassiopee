# - initVars (array) -
import Converter as C
a = C.array('x,y,z', 10, 10, 10)
a = C.initVars(a, 'celln', 2.)
a = C.initVars(a, ['F','G'], 0.)
C.convertArrays2File(a, 'out.plt')
