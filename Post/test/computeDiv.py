# - computeDiv (array) -
import Converter as C
import Post as P
import Generator as G

ni = 1001; nj = 1001; nk = 1
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, '{veloX}= 2*{x}+{x}*{y}')
m = C.initVars(m, '{veloY}= 4.*{y}')
m = C.initVars(m, '{veloZ}= {x}+{z}*{z}')
p = P.computeDiv(m, ['veloX','veloY','veloZ']) # p is defined on centers
p = C.center2Node(p) # back to initial mesh
p = C.addVars([m, p])
C.convertArrays2File([p], 'out.plt')
