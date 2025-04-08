# - hyper2D (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C

line = D.naca(12., 5001)
# Distribution
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))

a = G.hyper2D(line, distrib, "C")
C.convertPyTree2File(a, 'out.cgns')
