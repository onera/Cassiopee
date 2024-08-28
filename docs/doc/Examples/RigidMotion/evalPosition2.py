# - evalPosition (array) -
import RigidMotion as R
import Transform as T
import Generator as G
import Converter as C
from math import *

# Coordinates of rotation center in absolute frame
def centerAbs(t): return [t, 0, 0]

# Coordinates of rotation center in relative frame
def centerRel(t): return [5, 5, 0]

# Rotation matrix
def rot(t):
    omega = 0.1
    m = [[cos(omega*t), -sin(omega*t), 0],
         [sin(omega*t), cos(omega*t),  0],
         [0,            0,             1]]
    return m

# Complete motion
def F(t): return (centerAbs(t), centerRel(t), rot(t))

# Create a structured Cartesian grid
a = G.cart((0,0,0), (1,1,1), (11,11,1))
# Move the mesh
t = 3.
b = R.evalPosition(a, t, F)
C.convertArrays2File([b], "out1.plt")

# Equivalent to:
c = T.rotate(a, (5,5,0), (0,0,1), 0.1*t*180/3.14)
c = T.translate(c, (t-5,-5,0))
C.convertArrays2File([c], "out2.plt")
