# - Probe (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

LOCAL = test.getLocal()

# test case
a = G.cartRx((0,0,0), (1,1,1), (20,20,20), (3,3,3), depth=0, addCellN=False)
C._initVars(a, '{centers:F} = {centers:CoordinateX}')
t = C.newPyTree(['Base',a])

# create a probe with X
p1 = Probe.Probe(LOCAL+'/probe1.cgns', t, X=(10.,10.,10.), fields=['centers:F'], append=False, bufferSize=15)
for i in range(20):
    time = 0.1*i
    C._initVars(t, '{centers:F} = {centers:CoordinateX}+10.*sin(%20.16g)'%time)
    p1.extract(t, time=time)
p1.flush()

# test append
p1 = Probe.Probe(LOCAL+'/probe1.cgns', t, X=(10.,10.,10.), fields=['centers:F'], append=True, bufferSize=15)
for i in range(20):
    time = 2.+0.1*i
    C._initVars(t, '{centers:F} = {centers:CoordinateX}+10.*sin(%20.16g)'%time)
    p1.extract(t, time=time)
p1.flush()
test.testT(p1._probeZones, 1)

# create a probe with index
p1 = Probe.Probe(LOCAL+'/probe2.cgns', t, ind=(10,10,10), blockName='cart1-1-1', fields=['centers:F'], append=False, bufferSize=15)
for i in range(20):
    time = 0.1*i
    C._initVars(t, '{centers:F} = {centers:CoordinateX}+10.*sin(%20.16g)'%time)
    p1.extract(t, time=time)
p1.flush()
test.testT(p1._probeZones, 2)

# create probe from zones
p1 = Probe.Probe(LOCAL+'/probe3.cgns', fields=['centers:F'], append=False, bufferSize=15)
for i in range(20):
    time = 0.1*i
    a = G.cart((0,0,0), (1,1,1), (5,5,1))
    C._initVars(a, '{centers:F} = %20.16g'%time)
    p1.extract(a, time=time)
p1.flush()

p1 = Probe.Probe(LOCAL+'/probe3.cgns', fields=['centers:F'], append=True, bufferSize=15)
for i in range(20):
    time = 2+0.1*i
    a = G.cart((0,0,0), (1,1,1), (5,5,1))
    C._initVars(a, '{centers:F} = %20.16g'%time)
    p1.extract(a, time=time)
p1.flush()

# reread probe
p1 = Probe.Probe(LOCAL+'/probe1.cgns')
out = p1.read()
#C.convertPyTree2File(out, LOCAL+'/out1.cgns')
test.testT(out, 3)

p1 = Probe.Probe(LOCAL+'/probe2.cgns')
out = p1.read()
#C.convertPyTree2File(out, LOCAL+'/out2.cgns')
test.testT(out, 4)

p1 = Probe.Probe(LOCAL+'/probe3.cgns')
out = p1.read()
#C.convertPyTree2File(out, LOCAL+'/out3.cgns')
test.testT(out, 5)
