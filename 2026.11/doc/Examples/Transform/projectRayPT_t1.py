# - projectRay (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# structure
a = D.sphere((0,0,0), 1., 20)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin',a, 'jmax',[1,2])
C._initVars(a, 'F',1.); C._initVars(a, 'centers:G',2.)
b = G.cart((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
c = T.projectRay(b, a, (0,0,0))
test.testT(c,1)

# sur une zone non structuree
a = D.sphere((0,0,0), 1., 20)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin',a, 'jmax',[1,2])
C._initVars(a, 'F',1.); C._initVars(a, 'centers:G',2.)
b = G.cartTetra((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
c = T.projectRay(b, a, (0,0,0))
test.testT(c,2)

# NS sur une zone non structuree
a = D.sphere((0,0,0), 1., 20); a = C.convertArray2Tetra(a)
C._initVars(a, 'F',1.); C._initVars(a, 'centers:G',2.)
b = G.cartTetra((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
c = T.projectRay(b, a, (0,0,0))
test.testT(c,3)

# sur un arbre
t = C.newPyTree(['Base', 2, b])
tp = C.newPyTree(['Proj', 2, a])
t2 = T.projectRay(t, tp, (0,0,0))
test.testT(t2,4)
