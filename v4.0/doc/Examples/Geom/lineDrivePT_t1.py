# - lineDrive (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test
import Converter.PyTree as C

# 1D structure
a = D.naca(12.)
b = D.line((0,0,0), (0,0.,1.))
c = D.lineDrive(a, b)
t = C.newPyTree(['Base',2,c])
test.testT(t, 1)

# 1D structure + champ noeud
a = D.circle( (0,0,0), 1)
a = C.addVars(a, 'var')
b = D.line((0,0,0), (0,0,1))
c = D.lineDrive(a, b)
t = C.newPyTree(['Base',2,c])
test.testT(t, 2)

# 2D structure + champ en noeuds + champ en centres
a = G.cylinder((0,0,0), 1, 2, 360, 0, 1, (50,21,1))
a = C.addVars(a, 'var')
a = C.addVars(a, 'centers:var2')
b = D.line((0,0,0), (0,0.,1.))
c = D.lineDrive(a, b)
t = C.newPyTree(['Base',3,c])
test.testT(t, 3)
