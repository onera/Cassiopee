# - usurp (array) -
import Post as P
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

cyln = []
# Creation des cylindres
a1 = G.cylinder((0,0,0), 0, 2, 360, 0, 1., (100,2,10))
a1 = T.subzone(a1, (1,2,1), (a1[2],2,a1[4]))
cyln.append(a1)

a2 = G.cylinder((0,0,0), 0, 2, 90, 0, 0.5, (10,2,10))
a2 = T.translate(a2, (0,0,0.2))
a2 = T.subzone(a2, (1,2,1), (a2[2],2,a2[4]))
cyln.append(a2)

c1 = cyln[0]
ib1 = C.array('cellN', c1[2]-1, c1[3], c1[4]-1)
ib1 = C.initVars(ib1,'cellN', 1)
ib1[1][0,586] = 0.

c2 = cyln[1]
ib2 = C.array('cellN', c2[2]-1, c2[3], c2[4]-1)
ib2 = C.initVars(ib2, 'cellN', 1)

ibc = [ib1, ib2]

r = P.usurp(cyln, ibc)
if r is not None: test.testA(r,1)
