# - surface (PyTree) -
import Converter.PyTree as C
import Geom.PyTree as D

# User definition of parametric surface
def f(t,u):
    x = t+u; y = t*t+1+u*u; z = u
    return (x,y,z)

a = D.surface(f)
C.convertPyTree2File(a, 'out.cgns')
