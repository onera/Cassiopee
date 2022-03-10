import Generator.PyTree as G
import Converter.PyTree as C
# cartr2(X0,H,R,Xf)
a = G.cartr2((10,5,1), (1,1,1), (1.5,1.3,1.),(200,100,100))
# print(a)
C.convertPyTree2File(a, 'out.cgns')
