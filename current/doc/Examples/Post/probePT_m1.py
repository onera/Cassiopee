# - Probe (pyTree) -
# probe tous les modes
import Post.Probe as Probe
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import KCore.test as test

LOCAL = test.getLocal()

# test case
a = G.cartRx((0,0,0), (1,1,1), (20,20,20), (3,3,3), depth=0, addCellN=False, rank=Cmpi.rank, size=Cmpi.size)
C._initVars(a, '{centers:F} = {centers:CoordinateX}')
t = C.newPyTree(['Base', a])

# create a probe from X
p1 = Probe.Probe(LOCAL+'/probe1.cgns', t, X=(10.,10.,10.), fields=['centers:F'], bufferSize=15, append=False)
for i in range(20):
    time = 0.1*i
    C._initVars(t, '{centers:F} = {centers:CoordinateX}+10.*sin(%20.16g)'%time)
    p1.extract(t, time=time)
p1.flush()

# test append
p1 = Probe.Probe(LOCAL+'/probe1.cgns', t, X=(10.,10.,10.), fields=['centers:F'], bufferSize=15, append=True)
for i in range(20):
    time = 2.+0.1*i
    C._initVars(t, '{centers:F} = {centers:CoordinateX}+10.*sin(%20.16g)'%time)
    p1.extract(t, time=time)
p1.flush()

if Cmpi.rank == 0: test.testT(p1._probeZones, 1)

# create a probe from ind
p2 = Probe.Probe(LOCAL+'/probe2.cgns', t, ind=(2,2,2), blockName='cart0-0-0', fields=['centers:F'], bufferSize=15, append=False)
for i in range(20):
    time = 0.1*i
    C._initVars(t, '{centers:F} = {centers:CoordinateX}+10.*sin(%20.16g)'%time)
    p2.extract(t, time=time)
p2.flush()

if Cmpi.rank == 0: test.testT(p2._probeZones, 2)
Cmpi.barrier()

# create a probe from zones
p3 = Probe.Probe(LOCAL+'/probe3.cgns', fields=['centers:F'], append=False, bufferSize=15)
for i in range(20):
    time = 0.1*i
    a = G.cart((Cmpi.rank,0,0), (0.1,0.1,0.1), (11,11,1))
    a[0] = 'cart%d'%Cmpi.rank
    C._initVars(a, '{centers:F} = %20.16g'%time)
    p3.extract(a, time=time)
p3.flush()

p3 = Probe.Probe(LOCAL+'/probe3.cgns', fields=['centers:F'], append=True, bufferSize=15)
for i in range(20):
    time = 2+0.1*i
    a = G.cart((Cmpi.rank,0,0), (0.1,0.1,0.1), (11,11,1))
    a[0] = 'cart%d'%Cmpi.rank
    C._initVars(a, '{centers:F} = %20.16g'%time)
    p3.extract(a, time=time)
p3.flush()

if Cmpi.rank == 0: test.testT(p3._probeZones, 3)
Cmpi.barrier()
