# - getCEBBIntersectingDomains 3D (pyTree) -
import Connector.PyTree as X
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test
from math import cos, sin

# Coordonnees du centre de rotation dans le repere absolu
def centerAbs(t): return [t, 0, 0]

# Coordonnees du centre de la rotation dans le repere entraine
def centerRel(t): return [5, 5, 0]

# Matrice de rotation
def rot(t):
    omega = 30.
    m = [[cos(omega*t), -sin(omega*t), 0],
         [sin(omega*t), cos(omega*t),  0],
         [0,            0,             1]]
    return m

# Mouvement complet
def F(t): return (centerAbs(t), centerRel(t), rot(t))

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 2., (50,50,3))
# --- CL
a = C.addBC2Zone(a,'wall','BCWall','jmin')
# --- champ aux noeuds
t = C.newPyTree(['Cylindre', a])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
b = G.cylinder((1.5,0.,0.), 0.5, 1., 360., 0., 4., (50,50,3))
# --- champ aux centres
# --- CL
b = C.addBC2Zone(b,'wall','BCWall','jmin')
t[2][1][2].append(b); b[0]='cylinder2'
#
dt = 1.; Funcs = [F,[]]
b1 = G.cart((-2.,-2.,0.),(0.4,0.4,1),(11,11,4))
b2 = G.cart((-4.,-2.,0.),(0.4,0.4,1),(11,11,4))
b3 = G.cart((2.,-2.,0.),(0.4,0.4,1),(11,11,4))
b4 = G.cart((-2.,-2.,3.),(0.4,0.4,1),(11,11,4))
b5 = G.cart((-4.,-2.,3.),(0.4,0.4,1),(11,11,4))
b6 = G.cart((2.,-2.,3.),(0.4,0.4,1),(11,11,4))
#
t = C.addBase2PyTree(t,'Cart');t[2][2][2] = [b1,b2,b3,b4,b5,b6]
C._initVars(t, 'centers:cellN', 1.)
C._initVars(t, 'Density', 2.)
bases = Internal.getNodesFromType(t,'CGNSBase_t'); base = bases[0]
doms = X.getCEBBTimeIntersectingDomains(base, F, bases, Funcs, 0, 6, dt, sameBase=1)
test.testO(doms)
doms = X.getCEBBTimeIntersectingDomains(base, F, bases, Funcs, 0, 6, dt, sameBase=0)
test.testO(doms,2)
