# - copyGridInit2Grid (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import RigidMotion.PyTree as R

a = G.cart((0,0,0), (1,1,1), (10,10,10))
R._copyGrid2GridInit(a, mode=1)
R._copyGridInit2Grid(a)
C.convertPyTree2File(a, 'out.cgns')
