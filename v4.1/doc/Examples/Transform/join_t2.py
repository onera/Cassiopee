# - join (array) -
import Converter as C
import Transform as T
import Generator as G
import KCore.test as test

# test 2D
# raccord (im,1)
a1 = G.cart((0.,0.,0.), (1.,1.,1.), (11,11,1) )
a2 = G.cart((0.,10.,0.), (1.,1.,1.), (11,11,1) )
a2 = T.reorder(a2,(1,3,2))
a = T.join (a1, a2)
test.testA([a],1)

# avec champs en centres
a1 = C.initVars(a1,'{F}={x}+{y}**2'); a2 = C.initVars(a2,'{F}={x}+{y}**2')
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join (a1, a2, ac1, ac2)
test.testA(res,21)

# raccord (1,im)
a1 = G.cart((0.,0.,0.), (1.,1.,1.), (11,11,1) )
a1 = T.reorder(a1,(-1,2,3))
a2 = G.cart((10.,0.,0.), (1.,1.,1.), (11,11,1) )
a2 = T.reorder(a2,(-1,2,3))
a = T.join (a1, a2)
test.testA([a],2)

# avec champs en centres
a1 = C.initVars(a1,'{F}={x}+{y}**2'); a2 = C.initVars(a2,'{F}={x}+{y}**2')
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join (a1, a2, ac1, ac2)
test.testA(res,22)


# raccord (1,1)
a1 = G.cart((0.,0.,0.), (1.,1.,1.), (11,11,1) )
a1 = T.reorder(a1,(-1,2,3))
a2 = G.cart((10.,0.,0.), (1.,1.,1.), (11,11,1) )
a = T.join (a1, a2)
test.testA([a],3)

# avec champs en centres
a1 = C.initVars(a1,'{F}={x}+{y}**2'); a2 = C.initVars(a2,'{F}={x}+{y}**2')
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join (a1, a2, ac1, ac2)
test.testA(res,32)

# raccord (im,im)
a1 = G.cart((0.,0.,0.), (1.,1.,1.), (11,11,1) )
a2 = G.cart((10.,0.,0.), (1.,1.,1.), (11,11,1) )
a2 = T.reorder(a2,(-1,2,3))
a = T.join (a1, a2)
test.testA([a],4)

# avec champs en centres
a1 = C.initVars(a1,'{F}={x}+{y}**2'); a2 = C.initVars(a2,'{F}={x}+{y}**2')
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join (a1, a2, ac1, ac2)
test.testA(res,42)

# raccord (im,jm)
a1 = G.cart((0.,0.,0.), (1.,1.,1.), (11,11,1) )
a2 = G.cart((0.,10.,0.), (1.,1.,1.), (11,1,11) )
a2 = T.reorder(a2,(-2,-1,3))
a = T.join (a1, a2)
test.testA([a],5)
# avec champs en centres
a1 = C.initVars(a1,'{F}={x}+{y}**2'); a2 = C.initVars(a2,'{F}={x}+{y}**2')
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join (a1, a2, ac1, ac2)
test.testA(res,52)
