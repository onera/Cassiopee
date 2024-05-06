# - Probe (pyTree) -
# Exemple de probe mode=2 (empilement de surfaces avec champs en noeuds)
import Post.Probe as Probe
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import KCore.test as test

LOCAL = test.getLocal()

# create a probe from zones
p3 = Probe.Probe(LOCAL+'/probe3.cgns', fields=['F'], append=False, bufferSize=15)
for i in range(20):
    time = 0.1*i
    a = G.cart((Cmpi.rank,0,0), (0.1,0.1,0.1), (11,11,1))
    a[0] = 'cart%d'%Cmpi.rank
    C._initVars(a, '{F} = %20.16g'%time)
    p3.extract(a, time=time)
p3.flush()

p3 = Probe.Probe(LOCAL+'/probe3.cgns', fields=['F'], append=True, bufferSize=15)
for i in range(20):
    time = 2+0.1*i
    a = G.cart((Cmpi.rank,0,0), (0.1,0.1,0.1), (11,11,1))
    a[0] = 'cart%d'%Cmpi.rank
    C._initVars(a, '{F} = %20.16g'%time)
    p3.extract(a, time=time)
p3.flush()

if Cmpi.rank == 0: test.testT(p3._probeZones, 1)

out = p3.read(ind=[(0,0),(1,1)], probeName='cart0')
if Cmpi.rank == 0: test.testT(out, 2)
out = p3.read(cont=0)
if Cmpi.rank == 0: test.testT(out, 3)
