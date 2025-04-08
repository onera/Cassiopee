# - curve (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D

# User definition of parametric curve
def f(t):
    x = t; y = t*t+1; z = 0.
    return (x,y,z)

a = D.curve(f)
C.convertPyTree2File(a, 'out.cgns')
