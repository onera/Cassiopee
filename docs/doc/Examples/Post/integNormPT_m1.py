# - integNorm (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Post.PyTree as P
import Post.Mpi as Pmpi
import Converter.Mpi as Cmpi
import KCore.test as test

def f1(x,y): return 2*x + y
def f2(x,y): return 3*x*y + 4

ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m,'vx', f1, ['CoordinateX','CoordinateY'])
m = C.initVars(m,'vy', f2, ['CoordinateX','CoordinateY'])
res = Pmpi.integNorm(m, Internal.__FlowSolutionNodes__)
test.testO(res[0],1)

# Test avec des []
t = C.newPyTree(['Base'])
if Cmpi.rank == 0:
    m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
    m = C.initVars(m, 'vx', 1.)
    t[2][1][2] += [m]

res1 = Pmpi.integNorm(t, 'vx')
if Cmpi.rank == 0: test.testO(res1[0], 2)
