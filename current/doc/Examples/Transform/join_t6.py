# - join (array) -
import Transform as T
import Converter as C
import Generator as G
import KCore.test as test

# Join 2 NGON avec TETRA
a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.convertArray2NGon(a1)
a1 = G.close(a1)
a1 = C.initVars(a1, 'F', 2.)
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.convertArray2NGon(a2)
a2 = G.close(a2)
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 1)

# Join 2 NGON avec TETRA + champs en centres
a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.convertArray2NGon(a1)
ac1 = C.node2Center(a1)
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.convertArray2NGon(a2)
ac2 = C.node2Center(a2)
ac1 = C.initVars(ac1,'G',2.); ac1 = C.extractVars(ac1,['G'])
ac2 = C.initVars(ac2,'G',1.); ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res, 21)

# Join 2 NGON issu d un TETRA et d un HEXA
a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.convertArray2NGon(a1)
a1 = G.close(a1)
a1 = C.initVars(a1, 'F', 2.)
a2 = G.cartHexa((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.convertArray2NGon(a2)
a2 = G.close(a2)
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 3)
