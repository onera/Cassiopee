# - integ (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
C._initVars(m, 'vx', 1.); C._initVars(m, 'ratio', 1.)
resn = P.integ(m, 'vx'); print(resn)
#>> [99.99999999999989]
