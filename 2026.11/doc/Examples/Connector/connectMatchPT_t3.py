# - connectMatch sur meme bloc (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
t = C.newPyTree(['Base', a])
t = X.connectMatch(t)
test.testT(t,1)

# 3D raccord i = 1 partiel profil NACA
a = D.naca(12., 5001)
a2 = D.line((1.,0.,0.),(2.,0.,0.),5001); a = T.join(a, a2)
a2 = D.line((2.,0.,0.),(1.,0.,0.),5001); a = T.join(a2, a)
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
naca = G.hyper2D(a, distrib, "C")
a = T.addkplane(naca)
# --- champ aux centres
C._initVars(a, 'centers:Density', 1.)
# --- champ aux noeuds
C._initVars(a, 'cellN', 2.)
t = C.newPyTree(['Base', a])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = X.connectMatch(t)
test.testT(t,2)

# 2D raccord i=1 avec i = imax
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,1))
t = C.newPyTree(['Base',2,a])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = X.connectMatch(t, dim=2)
test.testT(t,3)

# 2D raccord i = 1 partiel
t = C.newPyTree(['Base',2, naca])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = X.connectMatch(t, dim=2)
test.testT(t,4)
