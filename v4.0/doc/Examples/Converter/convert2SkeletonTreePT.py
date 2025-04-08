# - convert2SkeletonTree (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = Cmpi.convert2SkeletonTree(a)
C.convertPyTree2File(a, 'out.cgns')
