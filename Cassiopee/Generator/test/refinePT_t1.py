# - refine (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# cas 2D multibloc avec CL et raccords
a = G.cart( (0,0,0), (0.1,0.1,0.1), (11,11,1))
b = G.cart( (0,1,0), (0.1,0.1,0.1), (11,11,1))
a = C.addBC2Zone(a,'match1','BCMatch','jmax',b,'jmin')
a = C.addBC2Zone(a,"wall",'BCWall','imin')
a = C.addBC2Zone(a,"overlap",'BCOverlap','imax')
b = C.addBC2Zone(b,'match2','BCMatch','jmin',a,'jmax')
b = C.addBC2Zone(b,"wall",'BCWall','imin')
b = C.addBC2Zone(b,"overlap",'BCOverlap','imax')
t = C.newPyTree(['Base',2]); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
# facteur de raffinement non entier
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],1.5,dir=1)
test.testT(t,1)
# facteur de raffinement entier 1 direction
t = C.newPyTree(['Base']); t[2][1][2]+=[a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],2.,dir=1)
test.testT(t,2)
#
# facteur de raffinement entier 3 directions
#
t = C.newPyTree(['Base',2]); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],2,dir=0)
test.testT(t,3)
