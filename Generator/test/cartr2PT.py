import Generator.PyTree as G
import Converter.PyTree as C
# cartr2(X0,H,R,Xf)
a = G.cartr2((10,2,3), (3,2,5), (1.5,1.3,1.),(200,179.5,160))
# print(a)
C.convertPyTree2File(a, 'out.cgns')
