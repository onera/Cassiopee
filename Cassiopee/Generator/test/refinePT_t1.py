# - refine (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# cas 2D multibloc avec CL et raccords
a = G.cart( (0,0,0), (0.1,0.1,0.1), (11,11,1))
b = G.cart( (0,1,0), (0.1,0.1,0.1), (11,11,1))
C._addBC2Zone(a,'match1','BCMatch','jmax',b,'jmin')
C._addBC2Zone(a,"wall",'BCWall','imin')
C._addBC2Zone(a,"overlap",'BCOverlap','imax')
C._addBC2Zone(b,'match2','BCMatch','jmin',a,'jmax')
C._addBC2Zone(b,"wall",'BCWall','imin')
C._addBC2Zone(b,"overlap",'BCOverlap','imax')
t = C.newPyTree(['Base',2]); t[2][1][2] += [a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
# facteur de raffinement non entier > 1
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],1.5,dir=1)
test.testT(t,1)

# facteur de raffinement entier 1 direction
t = C.newPyTree(['Base']); t[2][1][2]+=[a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],2.,dir=1)
test.testT(t,2)
#
# facteur de raffinement entier 3 directions
#
t = C.newPyTree(['Base',2]); t[2][1][2] += [a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],2,dir=0)
test.testT(t,3)

# facteur de raffinement (1/factor entier) 1 direction
t = C.newPyTree(['Base']); t[2][1][2]+=[a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],0.5,dir=1)
test.testT(t,4)

# facteur de raffinement (1/factor non entier) 1 direction
t = C.newPyTree(['Base']); t[2][1][2]+=[a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],0.3,dir=1)
test.testT(t,5)

# facteur de raffinement (1/factor entier) 3 directions
t = C.newPyTree(['Base']); t[2][1][2]+=[a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],0.5,dir=0)
test.testT(t,6)
C.convertPyTree2File(t,"out.cgns")

# facteur de raffinement (1/factor non entier) 3 directions
t = C.newPyTree(['Base']); t[2][1][2]+=[a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
for noz in range(2): t[2][1][2][noz] = G.refine(t[2][1][2][noz],0.3,dir=0)
test.testT(t,7)

