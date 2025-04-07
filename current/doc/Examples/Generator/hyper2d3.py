# - hyper2D3 (array) -
import Geom as D
import Generator as G
import Converter as C

a = D.naca(12., 5001)

# distribution
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 20./(Nj-1),1), (Ni,Nj,1))
#
msh = G.hyper2D3(a, distrib, "C", 10., 350.)
C.convertArrays2File([msh], 'out.plt')
