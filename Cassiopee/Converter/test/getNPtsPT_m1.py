# - getNPts (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import KCore.test as test

if Cmpi.rank == 0: a = G.cart((0,0,0), (1,1,1), (10,10,10)) # 1000
else: a = G.cart((0,0,0), (1,1,1), (20,10,10)) # 2000

npts = Cmpi.getNPts(a) # 3000

if Cmpi.rank == 0: test.testO(npts, 1)
