# - getSizeOf (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import numpy
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
s = Cmpi.getSizeOf(a, 'sum')
if Cmpi.rank == 0: test.testO(s, 1)

a = [numpy.zeros((100), dtype=numpy.float64)]
s = Cmpi.getSizeOf(a, 'sum')
if Cmpi.rank == 0: test.testO(s, 2)
