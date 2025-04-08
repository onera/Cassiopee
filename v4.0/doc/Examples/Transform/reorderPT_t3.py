# - reorder (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
a = C.addBC2Zone(a,'wall','BCWall','jmin')
a = C.addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
b = C.addBC2Zone(b,'wall','BCWall','jmin')
b = C.addBC2Zone(b,'overlap','BCOverlap','jmax')
a = C.addBC2Zone(a, 'match', 'BCMatch', [11,11,1,11,1,11], zoneDonor=b, rangeDonor=[1,1,1,11,1,11], trirac=[1,2,3])
b = C.addBC2Zone(b, 'match', 'BCMatch', [1,1,1,11,1,11], zoneDonor=a, rangeDonor=[11,11,1,11,1,11], trirac=[1,2,3])
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
#
# reorder l arbre
t2 = T.reorder(t,(-1,2,3))
test.testT(t2,1)
#
# reorder une zone et ne rend pas coherent le reste de l arbre
t[2][1][2][0] = T.reorder(t[2][1][2][0],(-1,-2,3))
test.testT(t,2)

# reorder une zone et rend coherent le reste de l'arbre
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
a = C.addBC2Zone(a,'wall','BCWall','jmin')
a = C.addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
b = C.addBC2Zone(b,'wall','BCWall','jmin')
b = C.addBC2Zone(b,'overlap','BCOverlap','jmax')
a = C.addBC2Zone(a, 'match', 'BCMatch', [11,11,1,11,1,11], zoneDonor=b, rangeDonor=[1,1,1,11,1,11], trirac=[1,2,3])
b = C.addBC2Zone(b, 'match', 'BCMatch', [1,1,1,11,1,11], zoneDonor=a, rangeDonor=[11,11,1,11,1,11], trirac=[1,2,3])
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t[2][1][2][0] = T.reorder(t[2][1][2][0],(-1,-2,3),t)
test.testT(t,3)

# reorder sur la base  et rend coherent le reste de l'arbre
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
a = C.addBC2Zone(a,'wall','BCWall','jmin')
a = C.addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
b = C.addBC2Zone(b,'wall','BCWall','jmin')
b = C.addBC2Zone(b,'overlap','BCOverlap','jmax')
a = C.addBC2Zone(a, 'match', 'BCMatch', [11,11,1,11,1,11], zoneDonor=b, rangeDonor=[1,1,1,11,1,11], trirac=[1,2,3])
b = C.addBC2Zone(b, 'match', 'BCMatch', [1,1,1,11,1,11], zoneDonor=a, rangeDonor=[11,11,1,11,1,11], trirac=[1,2,3])
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t[2][1] = T.reorder(t[2][1],(-1,-2,3),t)
test.testT(t,4)

# reorder sur la base  et rend coherent le reste de l arbre
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
a = C.addBC2Zone(a,'wall','BCWall','jmin')
a = C.addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
b = C.addBC2Zone(b,'wall','BCWall','jmin')
b = C.addBC2Zone(b,'overlap','BCOverlap','jmax')
a = C.addBC2Zone(a, 'match', 'BCMatch', [11,11,1,11,1,11], zoneDonor=b, rangeDonor=[1,1,1,11,1,11], trirac=[1,2,3])
b = C.addBC2Zone(b, 'match', 'BCMatch', [1,1,1,11,1,11], zoneDonor=a, rangeDonor=[11,11,1,11,1,11], trirac=[1,2,3])
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t[2][1][2][0] = T.reorder(t[2][1][2][0],(-2,1,3),t)
test.testT(t,5)
