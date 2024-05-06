# - cartRx (pyTree) -
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import KCore.test as test

# cartRx avec BCMatch
a = G.cartRx((0,0,0), (1,1,1), (10,10,10), (5,5,5), addBCMatch=True, rank=Cmpi.rank, size=Cmpi.size)

if Cmpi.rank == 0: test.testT(a, 1)
