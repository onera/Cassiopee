# - surface (array) -
import Converter as C
import Geom as D

# User definition of parametric surface by a function
def f(t,u):
    x = t+u; y = t*t+1+u*u; z = u
    return (x,y,z)

a = D.surface(f)

# Definition by formula
b = D.surface('{x} = cos(pi*{t}); {y} = sin(pi*{u}); {z} = {t}*{u}')
C.convertArrays2File([a, b], 'out.plt')
