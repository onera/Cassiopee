# - extractMesh  sur maillage surfacique -
import Converter.PyTree as C
import Post.PyTree as P
import Geom.PyTree as D
import KCore.test as test

# Fonction de surface
def FS(t,u):
    x = t+u
    y = t*t+1+u*u
    z = u
    return (x,y,z)

# Create an init function
def FI(x,y,z):
    return x*x + y*y + z*z

# Creation de la surface portant la solution
a = D.surface(FS, 50)
a = C.initVars(a, 'sol', FI, ['CoordinateX','CoordinateY','CoordinateZ'])

# Creation de la surface d'extraction
e = D.surface(FS, 100)
e = P.extractMesh([a], e, order=2, tol=1.e-3)
test.testT(e, 1)

# unstructured surface
a = C.convertArray2Tetra(a);
a = C.initVars(a, 'sol', FI, ['CoordinateX','CoordinateY','CoordinateZ'])
e = C.convertArray2Tetra(e)
e = P.extractMesh([a], e, order=2, tol=1.e-3)
test.testT(e, 2)

# arbre
t= C.newPyTree(['Base',2]); t[2][1][2] += [a]
a = C.initVars(a, 'sol', FI, ['CoordinateX','CoordinateY','CoordinateZ'])
e = P.extractMesh(t, e, order=2, tol=1.e-3)
test.testT(e, 3)
