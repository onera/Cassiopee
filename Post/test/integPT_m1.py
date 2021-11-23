# - integ (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Converter.Mpi as Cmpi
import Post.Mpi as Pmpi
import KCore.test as test

ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m, 'vx', Cmpi.rank)
m = C.initVars(m, 'vz', 1.)
m = C.initVars(m, 'centers:vy', 2.)
res1 = Pmpi.integ(m, 'vx')
res2 = Pmpi.integ(m,'centers:vy')

if Cmpi.rank == 0: test.testO(res1+res2, 1)
