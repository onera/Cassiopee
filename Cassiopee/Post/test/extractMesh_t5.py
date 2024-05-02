# - extractMesh  sur maillage surfacique -
import Converter as C
import Post as P
import Geom as D
import KCore.test as test

# Fonction de surface
def FS(t,u):
    x = t+u
    y = t*t+1+u*u
    z = u
    return (x,y,z)

# Creation de fonction d'initialisation
def FI(x,y,z):
    return x*x + y*y + z*z

# Creation de la surface portant la solution
a = D.surface(FS, 50)
a = C.initVars(a, 'sol', FI, ['x','y','z'])

# Creation de la surface d'extraction
e = D.surface(FS, 100)
e2 = P.extractMesh([a], e, order=2, tol=1.e-3)
test.testA([a,e2], 1)
e2 = P.extractMesh([a], [e], order=2, tol=1.e-3)
test.testA(e2, 2)
