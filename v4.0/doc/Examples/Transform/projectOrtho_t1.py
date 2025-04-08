# - projectOrtho (array) -
import Geom as D
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

#a = D.sphere((5,0,0), 5., 200)
#test.stdTestA(T.projectOrtho, [a])

a = D.sphere((0,0,0), 1., 200)
a = C.initVars(a,'F',1)
b = G.cartNGon((-0.5,-0.5,-0.5),(0.05,0.05,0.1), (10,10,1))
b = C.initVars(b,'F',1)
c = T.projectOrtho([b], [a])
test.testA([a,b]+c, 1)
