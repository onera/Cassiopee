# - projectAllDirs (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# Structure
a = D.sphere((0,0,0), 1., 20)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin',a,'jmax',[1,2])
C._addVars(a, 'F'); C._addVars(a,'centers:G')
b = G.cart((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
b = G.getNormalMap(b); b = C.center2Node(b,['centers:sx','centers:sy','centers:sz'])
c = T.projectAllDirs(b, a, ['sx','sy','sz']); c[0] = 'projection'
t = C.newPyTree(['Base',2,c])
test.testT(t,1)

# Non structure sur structure
a = D.sphere((0,0,0), 1., 20)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin',a,'jmax',[1,2])
C._addVars(a, 'F'); C._addVars(a, 'centers:G')
b = G.cartTetra((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
b = G.getNormalMap(b); b = C.center2Node(b,['centers:sx','centers:sy','centers:sz'])
c = T.projectAllDirs(b, a, ['sx','sy','sz']); c[0] = 'projection'
t = C.newPyTree(['Base',2,c])
test.testT(t,2)

# Structure sur NS
a = D.sphere((0,0,0), 1., 20)
a = C.convertArray2Tetra(a)
b = G.cart((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
C._initVars(b, 'F',2.); C._initVars(b, 'centers:G',1.)
b = G.getNormalMap(b); b = C.center2Node(b,['centers:sx','centers:sy','centers:sz'])
c = T.projectAllDirs(b, a, ['sx','sy','sz']); c[0] = 'projection'
t = C.newPyTree(['Base',2]); t[2][1][2] += [c, b]
test.testT(t,3)

# NS/NS
a = D.sphere((0,0,0), 1., 20)
a = C.convertArray2Tetra(a)
b = G.cartHexa((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
C._initVars(b, 'F',2.); C._initVars(b, 'centers:G',1.)
b = G.getNormalMap(b); b = C.center2Node(b,['centers:sx','centers:sy','centers:sz'])
c = T.projectAllDirs(b, a, ['sx','sy','sz']); c[0] = 'projection'
t = C.newPyTree(['Base',2]); t[2][1][2] += [c, b]
test.testT(t,4)

# sur un arbre
t = C.newPyTree(['Base',2,a])
tp = C.newPyTree(['Proj',2,b])
t = G.getNormalMap(t)
t = C.center2Node(t,['centers:sx','centers:sy','centers:sz'])
t2 = T.projectAllDirs(t, tp, ['sx','sy','sz']); c[0] = 'projection'
test.testT(t2,5)
