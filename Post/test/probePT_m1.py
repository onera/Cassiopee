# - Probe (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Mpi as Cmpi

# test case
a = G.cartRx((0,0,0), (1,1,1), (30,30,30), (5,5,5), depth=0, addCellN=False, rank=Cmpi.rank, size=Cmpi.size)
C._initVars(a, '{centers:F} = {centers:CoordinateX}')
t = C.newPyTree(['Base',a])
#Cmpi.convertPyTree2File(t, 'out.cgns')

# create a probe
p1 = Probe.Probe('probe1.cgns', t, (10.,10.,10.), fields=['centers:F'], append=False)
p1.print()
for i in range(110):
    p1.extract(time=0.1*i)
p1.flush()

# test append
p1 = Probe.Probe('probe1.cgns', t, (10.,10.,10.), fields=['centers:F'], append=True)
for i in range(50):
    p1.extract(time=11.+0.1*i)
p1.flush()
