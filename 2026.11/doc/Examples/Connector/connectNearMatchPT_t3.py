# - connectNearMatch (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,51,30))
a1 = T.subzone(a,(1,1,1),(25,51,30))
a2 = T.subzone(a,(25,1,1),(50,51,30))
a2 = T.oneovern(a2,(1,2,1))
t = C.newPyTree(['Base']); t[2][1][2]+=[a1,a2]
t = X.connectNearMatch(t)
test.testT(t,1)

# 3D raccord i = 1 partiel profil NACA
a = D.naca(12., 5001)
a2 = D.line((1.,0.,0.),(2.,0.,0.),5001); a = T.join(a, a2)
a2 = D.line((2.,0.,0.),(1.,0.,0.),5001); a = T.join(a2, a)
Ni = 301; Nj = 51
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
naca = G.hyper2D(a, distrib, "C")
a = T.addkplane(naca)
a1 = T.subzone(a,(1,1,1),(151,Nj,2))
a2 = T.subzone(a,(151,1,1),(Ni,Nj,2))
a2 = T.oneovern(a2,(2,2,1))
t = C.newPyTree(['Base']); t[2][1][2]+=[a1,a2]
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
# --- champ aux centres
t = C.initVars(t, 'centers:Density', 1.)
# --- champ aux noeuds
t = C.initVars(t, 'cellN', 2.)
t = X.connectNearMatch(t,dim=2)
test.testT(t,2)
