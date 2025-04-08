# - hyper2D4 (array) -
import Geom as D
import Generator as G
import Converter as C

# Put a naca profile in a
a = D.naca(12., 5001)

# distribution
Ni = 100; Nj = 100
distrib = G.cart((0,0,0), (1./(Ni-1), 20./(Nj-1),1), (Ni,Nj,1))
msh = G.hyper2D4(a, distrib, "C")
C.convertArrays2File([msh], 'out.plt')
