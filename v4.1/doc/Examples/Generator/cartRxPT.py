# - cartRx (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cartRx((0,0,0), (1,1,1), (30,30,30), (5,5,5), depth=0, addCellN=False)
C.convertPyTree2File(a, 'out.cgns')
