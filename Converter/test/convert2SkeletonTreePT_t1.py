# - convert2SkeletonTree (pyTree) -
import Converter.Mpi as Cmpi
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = Cmpi.convert2SkeletonTree(a)
test.testT(a, 1)
