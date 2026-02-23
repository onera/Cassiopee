# - integNormProduct (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m,'MomentumX',1.)
m = C.initVars(m,'MomentumY',1.)
m = C.initVars(m,'MomentumZ',1.)
res = P.integNormProduct(m,['MomentumX','MomentumY','MomentumZ']); print(res)
#>> 100.0
