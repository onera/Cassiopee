# - Probe (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import KCore.test as test

# test case
a = G.cartRx((0,0,0), (1,1,1), (20,20,20), (3,3,3), depth=0, addCellN=False, rank=Cmpi.rank, size=Cmpi.size)
C._initVars(a, '{centers:F} = {centers:CoordinateX}')
t = C.newPyTree(['Base', a])

# create a probe
p1 = Probe.Probe('probe1.cgns', t, X=(10.,10.,10.), fields=['centers:F'], append=False)
for i in range(20):
    time = 0.1*i
    C._initVars(t, f'{{centers:F}} = {{centers:CoordinateX}}+10.*sin({time})')
    p1.extract(t, time=time)
p1.flush()

# test append
p1 = Probe.Probe('probe1.cgns', t, X=(10.,10.,10.), fields=['centers:F'], append=True)
for i in range(20):
    time = 2.+0.1*i
    C._initVars(t, f'{{centers:F}} = {{centers:CoordinateX}}+10.*sin({time})')
    p1.extract(t, time=time)
p1.flush()

if Cmpi.rank == 0: test.testT(p1._probeZones, 1)

# create a probe
p2 = Probe.Probe('probe2.cgns', t, ind=(2,2,2), blockName='cart0-0-0', fields=['centers:F'], append=False)
for i in range(20):
    time = 0.1*i
    C._initVars(t, f'{{centers:F}} = {{centers:CoordinateX}}+10.*sin({time})')
    p2.extract(t, time=time)
p2.flush()

if Cmpi.rank == 0: test.testT(p2._probeZones, 2)

