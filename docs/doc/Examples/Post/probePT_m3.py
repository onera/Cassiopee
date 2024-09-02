# - Probe (pyTree) -
# Exemple de probe mode=0
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

#out = p1.read(ind=0)
#if Cmpi.rank == 0: test.testT(out, 2)

#out = p1.read(cont=0)
#if Cmpi.rank == 0: test.testT(out, 3)
#Cmpi.convertPyTree2File(out, 'out.cgns')
