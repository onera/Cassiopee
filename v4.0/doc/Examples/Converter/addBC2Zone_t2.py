# - addBC2Zone (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# - Structure Near Match -
# cas ou le transform est (1,2,3)
a = G.cart((0,0,0), (0.1,0.1,1.), (11,11,11))
b = G.cart((1,0,0), (0.1,0.05,1.), (11,21,11)); b[0] = 'cart2'
a = C.initVars(a,'F1',1.); a = C.initVars(a,'centers:G1',2.)
b = C.initVars(b,'F2',3.); b = C.initVars(b,'centers:G2',4.)
a = C.addBC2Zone(a,'nearmatch', 'BCNearMatch', 'imax',b,'imin',[1,2,3])
b = C.addBC2Zone(b,'nearmatch', 'BCNearMatch', 'imin',a,'imax',[1,2,3])
t = C.newPyTree(['Base', a, b])
test.testT(t)

# transform different de (1,2,3)
a = C.addBC2Zone(a,'nearmatch1', 'BCNearMatch', 'imax',b,'jmax',[-2,1,3])
b = C.addBC2Zone(b,'nearmatch2', 'BCNearMatch', 'jmax',a,'imax',[2,-1,3])
t = C.newPyTree(['Base', a, b])
test.testT(t,2)

# ajout par une subzone
a = G.cart((2,0,0), (0.1,0.1,1), (10,10,10))
a = C.addBC2Zone(a, 'stage', 'BCStage', 'imin')
test.testT(a,3)

# BCStage by a family name
a = G.cart((2,0,0), (0.1,0.1,1), (10,10,10))
a = C.addBC2Zone(a, 'stage', 'FamilySpecified:StageMxpl','imin')
test.testT(a,4)
