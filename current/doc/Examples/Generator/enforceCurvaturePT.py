# - enforceCurvature (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C

a = D.naca(12., 501)

# Distribution on the profile
Ni = 20; Nj = 20; Nk = 1; h = 1./(Ni-1)
b = G.cart((0,0,0), (h, 0.25/Nj,1), (Ni,Nj,Nk))
b = G.enforceCurvature(b, a, 0.6)
C.convertPyTree2File(b, 'out.cgns')
