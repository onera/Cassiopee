# - stack (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

a = G.cylinder((0,0,0), 1, 1.3, 360, 0, 1., (50,10,1))
b = T.rotate(a, (0,0,0), (1,0,0), 5.)
b = T.translate(b, (0,0,0.5))
c = G.stack(a, b)
test.testT(c, 1)
