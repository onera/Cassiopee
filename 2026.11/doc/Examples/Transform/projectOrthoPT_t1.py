# - projectOrtho (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# sur une zone structuree
a = D.sphere((0,0,0), 1., 20)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin',a, 'jmax',[1,2])
b = G.cart((-0.5,-0.5,-0.5),(0.05,0.05,0.1), (20,20,1))
C._addVars(b, 'F'); C._addVars(b, 'centers:G')
c = T.projectOrtho(b, a)
test.testT(c,1)

# sur une zone non structuree
a = D.sphere((0,0,0), 1., 20)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin',a, 'jmax',[1,2])
b = G.cartTetra((-0.5,-0.5,-0.5),(0.05,0.05,0.1), (20,20,1))
C._initVars(b, 'F', 2.); C._addVars(b, 'centers:G')
c = T.projectOrtho(b, a)
test.testT(c,2)

# NS sur une zone non structuree
a = D.sphere((0,0,0), 1., 20); a = C.convertArray2Tetra(a)
b = G.cartTetra((-0.5,-0.5,-0.5),(0.05,0.05,0.1), (20,20,1))
C._initVars(b, 'F',2.); C._addVars(b, 'centers:G')
c = T.projectOrtho(b, a)
t  = C.newPyTree(['Base',a,b,c])
test.testT(c,3)

# test sur un arbre
t = C.newPyTree(['Base',2]); t[2][1][2]+=[b]
tp = C.newPyTree(['Proj',2]); tp[2][1][2]+=[a]
t2 = T.projectOrtho(t,tp)
test.testT(t2,4)
