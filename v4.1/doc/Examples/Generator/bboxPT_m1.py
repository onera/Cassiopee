# - bbox (pyTree) -
import Generator.PyTree as G
import Generator.Mpi as Gmpi
import Converter.Mpi as Cmpi
import KCore.test as test

# bbox de deux zones en parallele
a = G.cart((Cmpi.rank*20,0.,0.),(0.1,0.1,0.1),(11,11,11))
bb = Gmpi.bbox(a)
if Cmpi.rank == 0: test.testO(bb, 1)
