# - merge (array) -
import Converter as C
import Transform as T
import Geom as D
def f(t,u):
    x = t+u; y = t*t+1+u*u; z = u
    return (x,y,z)

a = D.surface(f)
b = T.splitSize(a, 100)
b = T.merge(b)
C.convertArrays2File(b, "out.plt")
