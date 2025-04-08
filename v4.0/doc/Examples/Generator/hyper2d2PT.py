# - hyper2D2 (PyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G

msh = D.naca(12.,5001)
# Distribution
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
a = G.hyper2D2(msh, distrib, "C", 110.)
t = C.newPyTree(['Base',2,a])
C.convertPyTree2File(t, 'out.cgns')
