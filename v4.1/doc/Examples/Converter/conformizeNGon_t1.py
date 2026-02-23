# - conformizeNGon (array) -
import Generator as G
import Converter as C
import Transform as T
import KCore.test as test
a = G.cartNGon((0,0,0),(0.1,0.1,1),(11,11,1))
b = G.cartNGon((1.,0,0),(0.1,0.2,1),(11,6,1))
res = T.join(a,b)
res = C.conformizeNGon(res)
test.testA([res],1)
