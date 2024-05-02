# - usurp (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P
import KCore.test as test

a = G.cart((-5,-5,0.), (1,1,1), (10,10,1))
C._initVars(a, 'centers:cellN=1')
b = G.cylinder((0,0,0.), 1, 2, 360, 0, 1, (30,30,1))
C._initVars(b, 'centers:cellN=1')
t = C.newPyTree(['Base',2,a,b])
try:
    P._usurp(t)
    f = P.integ(t, 'centers:cellN')
    test.testT(t, 1)
except: pass
