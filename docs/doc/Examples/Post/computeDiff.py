# - computeDiff (array) -
import Converter as C
import Post as P
import Generator as G

def F(x):
    return (x > 5.)

ni = 30; nj = 40; nk = 1
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1.), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x'], isVectorized=True)
p = P.computeDiff(m,'ro')
p = C.addVars([m, p])
C.convertArrays2File([p], 'out.plt')
