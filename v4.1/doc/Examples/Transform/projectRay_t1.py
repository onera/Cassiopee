# - projectRay (array) -
import Geom as D
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

a = D.sphere((0,0,0), 1., 20); a = C.initVars(a,'F',1)
b = G.cart((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
b = C.initVars(b,'F',1)
c = T.projectRay([b],[a],(0,0,0))
test.testA([a,b]+c, 1)
