# - cartRx2 (pyTree) -
import Generator.PyTree as G
import Converter.Mpi as Cmpi

a = G.cartRx2((0,0,0), (10,10,10), (0.1,0.1,0.1), (-10,-10,-10), (30,20,20), (1.1,1.1,1.1), rank=Cmpi.rank, size=Cmpi.size)

Cmpi.convertPyTree2File(a, 'out.cgns')
