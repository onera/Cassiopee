# - curve (array) -
import Converter as C
import Geom as D

# Definition of parametric curve by a function
def f(t):
    x = t; y = t*t+1; z = 0.
    return (x,y,z)
a = D.curve(f)

# Definition by equation
b = D.curve('{x}=cos(2*pi*{t}); {y}=sin(2*pi*{t}); {z} = 0.')

# Definition from data base
from Geom.Parametrics import base
c = D.curve(base['circle'])
C.convertArrays2File([a,b], "out.plt")
