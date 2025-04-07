# - evalPosition (PyTree) -
import RigidMotion.PyTree as R
import Generator.PyTree as G
import Converter.PyTree as C
from math import *

# Coordonnees du centre de rotation dans le repere absolu
def centerAbs(t): return [t, 0, 0]

# Coordonnees du centre de la rotation dans le repere entraine
def centerRel(t): return [5, 5, 0]

# Matrice de rotation
def rot(t):
    omega = 0.1
    m = [[cos(omega*t), -sin(omega*t), 0],
         [sin(omega*t), cos(omega*t),  0],
         [0,            0,             1]]
    return m

# Mouvement complet
def F(t): return (centerAbs(t), centerRel(t), rot(t))

a = G.cart((0,0,0), (1,1,1), (11,11,2))

# Move the mesh
time = 3.
b = R.evalPosition(a, time, F); b[0]='moved'
C.convertPyTree2File([a,b], "out.cgns")
