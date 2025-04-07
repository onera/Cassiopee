# - join (array) -
import Geom as D
import Transform as T
import Converter as C
import KCore.test as test

# test1 lignes
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a2 = T.reorder(a2,(-1,2,3))
a1 = C.initVars(a1, 'F', 2)
a2 = C.initVars(a2, 'F', 3)
a = T.join(a1, a2)
test.testA([a],1)
# avec champs en centres
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join(a1, a2, ac1, ac2)
test.testA(res,21)

# test2 lignes
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a2 = T.reorder(a2,(2,-1,3))
a1 = C.initVars(a1, 'F', 2)
a2 = C.initVars(a2, 'F', 3)
a = T.join (a1, a2)
test.testA([a],2)
# avec champs en centres
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join (a1, a2, ac1, ac2)
test.testA(res,22)

# test3 lignes
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a1 = C.initVars(a1, 'F', 2)
a2 = C.initVars(a2, 'F', 3)
a2 = T.reorder(a2,(2,3,-1))
a = T.join (a1, a2)
test.testA([a],3)
# avec champs en centres
ac1 = C.node2Center(a1); ac2 = C.node2Center(a2)
res = T.join (a1, a2, ac1, ac2)
test.testA(res,23)
