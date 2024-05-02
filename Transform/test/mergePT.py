# - merge (pyTree) -
import Converter.PyTree as C
import Transform.PyTree as T
import Connector.PyTree as X
import Geom.PyTree as D

def f(t,u):
    x = t+u; y = t*t+1+u*u; z = u
    return (x,y,z)

a = D.surface(f)
b = T.splitSize(a, 100)
b = X.connectMatch(b, dim=2)
t = C.newPyTree(['Surface', b])
b = T.merge(t)
t[2][1][2] = b
C.convertPyTree2File(t, "out.cgns")
