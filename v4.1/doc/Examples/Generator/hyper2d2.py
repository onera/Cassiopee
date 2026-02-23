# - hyper2D2 (array) -
import Converter as C
import Geom as D
import Generator as G

msh = D.naca(12.,5001)

# Distribution
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
a = G.hyper2D2(msh, distrib, "C", 110.)
C.convertArrays2File([a],"out.plt")
