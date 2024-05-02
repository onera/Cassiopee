# - hyper2D2 (pyTree) -
import KCore.test as test
import Geom.PyTree as D
import Generator.PyTree as G

msh = D.naca(12,5001)

# Distribution
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
a = G.hyper2D2(msh, distrib, "C", 110.)
test.testT(a,1)
