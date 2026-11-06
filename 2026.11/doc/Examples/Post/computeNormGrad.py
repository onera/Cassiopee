# - computeNormGrad (array) -
import Converter as C
import Post as P
import Generator as G

ni = 11; nj = 11; nk = 1
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, '{ro}= 2*{x}+{x}*{y}')
p = P.computeNormGrad(m,'ro') # p is defined on centers
p = C.center2Node(p) # back on initial mesh
p = C.addVars([m, p])
C.convertArrays2File([p], 'out.plt')
