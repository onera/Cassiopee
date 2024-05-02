# - join (array) -
import Transform as T
import Generator as G
import Converter as C
import KCore.test as test

# test 3D
a1 = G.cart((0.,0.,0.), (1.,1.,1.), (11,11,11) )
a2 = G.cart((0.,10.,0.), (1.,1.,1.), (11,11,11) )
a2 = T.reorder(a2,(3,-1,2))
a2 = T.reorder(a2,(2,-3,-1))
a1 = C.initVars(a1, 'F', 2)
a2 = C.initVars(a2, 'F', 3)
a = T.join (a1, a2)
test.testA([a],1)
# avec champs en centres
a1 = C.initVars(a1,'{F}={x}+{y}**2+{z}**3'); a2 = C.initVars(a2,'{F}={x}+{y}**2+{z}**3')
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join(a1, a2, ac1, ac2)
test.testA(res,2)
