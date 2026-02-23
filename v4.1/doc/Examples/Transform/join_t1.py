# - join (array) -
import Geom as D
import Transform as T
import Converter as C
import Generator as G
import KCore.test as test

# Join 2 NON-STRUCT TETRA
a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.initVars(a1, 'F', 2)
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 1)

# Join 2 NON-STRUCT TETRA avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',2.)
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',2.)
ac1 = C.extractVars(ac1,['G'])
ac2 = C.extractVars(ac2,['G'])

res = T.join(a1, a2, ac1, ac2)
test.testA(res, 21)

# mix struct 3D + HEXA
a1 = G.cart((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.initVars(a1, 'F', 2)
a2 = G.cartHexa((10.,0.,0.), (1.,1.,1), (11,11,10))
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 2)

# mix struct 3D + HEXA avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',4.);ac1 = C.extractVars(ac1,['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.);ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res, 22)

# Join 2 NON-STRUCT PENTA
a1 = G.cartPenta((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.initVars(a1, 'F', 2)
a2 = G.cartPenta((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 9)

# Join 2 NON-STRUCT PENTA avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',4.);ac1 = C.extractVars(ac1,['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.);ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res, 92)

# Join 2 NON-STRUCT PYRA
a1 = G.cartPyra((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.initVars(a1, 'F', 2)
a2 = G.cartPyra((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 10)

# Join 2 NON-STRUCT PYRA avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',4.);ac1 = C.extractVars(ac1,['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.);ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res, 102)

# Join 2 NON-STRUCT TRI
a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,1))
a1 = C.initVars(a1, 'F', 2)
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,1))
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 3)

# Join 2 NON-STRUCT TRI avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',4.);ac1 = C.extractVars(ac1,['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.);ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res, 32)

# Join 1 STRUCT et 1 NON-STRUCT QUAD
a1 = G.cart((13.,0.,0.), (1.,1.,1), (11,11,1))
a1 = C.initVars(a1, 'F', 2)
a2 = G.cartHexa((9.,0.,0.), (1.,1.,1), (10,10,1))
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 4)

# Join 1 STRUCT et 1 NON-STRUCT QUAD avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',4.);ac1 = C.extractVars(ac1,['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.);ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res,42)

# Join 2 SRUCT-1D
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a1 = C.initVars(a1, 'F', 2)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a1, a2)
test.testA([a], 5)

# Join 2 STRUCT-1D avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1, 'G', 4.)
ac1 = C.extractVars(ac1, ['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.)
ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA([res[1]], 52)

# Join 2 BARS
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a1 = C.convertArray2Tetra(a1); a1 = C.initVars(a1, 'F', 2)
a2 = C.convertArray2Tetra(a2); a2 = C.initVars(a2, 'F', 3.)
a = T.join (a1, a2)
test.testA([a], 6)

# Join 2 BARS avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',4.)
ac1 = C.extractVars(ac1,['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.)
ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res,62)

# Join 1 BAR et 1 i-array
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a2 = C.convertArray2Tetra(a2)
a1 = C.initVars(a1, 'F', 2)
a2 = C.initVars(a2, 'F', 3.)
a = T.join (a1, a2)
test.testA([a], 7)

# Join 1 BAR et 1 i-array avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',4.);ac1 = C.extractVars(ac1,['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.);ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res,72)

# Join 2 arrays nodes
a1 = G.cart((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.initVars(a1, 'F', 2)
a1 = C.convertArray2Node(a1)
a2 = G.cart((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.initVars(a2, 'F', 3.)
a2 = C.convertArray2Node(a2)
a = T.join(a1, a2)
test.testA([a], 8)

# Join 2 arrays nodes avec champs en centres
ac1 = C.node2Center(a1); ac1 = C.initVars(ac1,'G',4.);ac1 = C.extractVars(ac1,['G'])
ac2 = C.node2Center(a2); ac2 = C.initVars(ac2,'G',5.);ac2 = C.extractVars(ac2,['G'])
res = T.join(a1, a2, ac1, ac2)
test.testA(res,82)
