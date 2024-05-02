# - computeGrad2 (array) -
import Converter as C
import Post as P
import Generator as G

m = G.cartNGon((0,0,0), (1,1,1), (4,4,4))

mc = C.node2Center(m)
mc = C.initVars(mc, '{ro}= 2*{x}+{x}*{y}')
#mc = C.initVars(mc, '{ro}=1.')
mv = C.extractVars(mc, ['ro'])

p = P.computeGrad2(m, mv) # p is defined on centers
p = C.center2Node(p) # back on initial mesh
p = C.addVars([m, p])
C.convertArrays2File([p], 'out.plt')
