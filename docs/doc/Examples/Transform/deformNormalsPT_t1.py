# - deformNormals (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test

# Zone structuree
a = G.cart((0,0,0), (0.1,0.1,1.), (11,11,1))
a = C.fillEmptyBCWith(a, 'wall', 'BCWall', dim=2)
a = C.initVars(a, '{F}={CoordinateX}+{CoordinateY}**2'); a = C.initVars(a,'{G}=1.')
a = C.initVars(a,'{alpha}=0.5')
a = T.deformNormals(a,'alpha',niter=1)
test.testT(a, 1)

# Liste de zones structuree
a = D.sphere6( (0,0,0), 1.,20)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = X.connectMatch(a, dim=2)
a = C.initVars(a, '{F}={CoordinateX}+{CoordinateY}**2'); a = C.initVars(a,'{G}=1.')
a = C.initVars(a, '{alpha}=0.5')
a = T.deformNormals(a, 'alpha', niter=1)
test.testT(a, 2)

# Zone non structuree
a = D.sphere( (0,0,0), 1.,50)
a = C.initVars(a, '{F}={CoordinateX}+{CoordinateY}**2'); a = C.initVars(a,'{G}=1.')
a = C.initVars(a,'{alpha}=0.5')
a = T.deformNormals(a,'alpha',niter=1)
test.testT(a, 3)

# Arbre
a = D.sphere6( (0,0,0), 1.,20)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = X.connectMatch(a, dim=2)
a = C.initVars(a, '{F}={CoordinateX}+{CoordinateY}**2'); a = C.initVars(a, '{G}=1.')
a = C.initVars(a,'{alpha}=0.5')
t = C.newPyTree(['Base',2]); t[2][1][2]+=a
t = T.deformNormals(t,'alpha',niter=1)
test.testT(t, 4)
