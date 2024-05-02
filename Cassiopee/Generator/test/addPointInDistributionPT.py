# - addPointInDistribution (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C

# Distribution
Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,2))
b = G.addPointInDistribution(a, Ni)
C.convertPyTree2File(b, 'out.cgns')
