# - streamLine2 (pyTree) -
import Converter.PyTree as C
import Post.Mpi as Pmpi
import Post.PyTree as P
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import KCore.test as test

a = G.cartRx3((0,0,0), (10,10,10), (0.1,0.1,0.1), (-10,-10,-10), (30,20,20), (1.1,1.1,1.1), rank=Cmpi.rank, size=Cmpi.size)

C._initVars(a, '{VelocityX}=1.')
C._initVars(a, '{VelocityY}=cos(0.2*{CoordinateX})')
C._initVars(a, '{VelocityZ}=0.')
#Cmpi.convertPyTree2File(a, 'out.cgns')

x0=0.1; y0=2; z0=5
points = []
for i in range(100): points.append( (x0,y0+0.1*i,z0) )
s = Pmpi.streamLine2(a, points, ['VelocityX','VelocityY','VelocityZ'])
#Cmpi.convertPyTree2File(s, 'out.cgns')
if Cmpi.rank == 0: test.testT(s, 1)