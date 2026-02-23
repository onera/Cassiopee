# - probe (pyTree) -
# Exemple de probe mode=3 (surface permeable)
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Geom.PyTree as D
import Post.Probe as Probe
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import KCore.test as test

LOCAL = test.getLocal()

if Cmpi.rank == 0:
    # Domaine donneur
    a = G.cart((0,0,0),(0.1,0.1,0.1),(51,51,51))
    C._initVars(a, "{F}={CoordinateX}*{CoordinateY}")
    C._initVars(a, 'cellN', 1.)
    C.convertPyTree2File(a, LOCAL+'/tc.cgns')

    # Surface de probe permeable
    s = D.sphere((2.5,2.5,2.5), 1.)
    C.convertPyTree2File(s, LOCAL+'/ts.cgns')
Cmpi.barrier()

# Prepare la probe
h = Filter.Handle(LOCAL+'/ts.cgns')
ts = h.loadAndSplit()
h = Filter.Handle(LOCAL+'/tc.cgns')
tc = h.loadAndSplit()

probe = Probe.Probe(LOCAL+'/probe.cgns', tPermeable=ts, fields=['F'], bufferSize=15)
tcs = probe.prepare(tc)

# Extraction
for i in range(20):
    time = 0.1*i
    C._initVars(tcs, "{F}={CoordinateX}*{CoordinateY}+%g"%time)
    probe.extract(tcs, time)
probe.flush()

if Cmpi.rank == 0:
    test.testT(probe._probeZones, 1)
    test.testT(probe._ts, 4)
Cmpi.barrier()

# Reread
#out = probe.read(ind=1)
#if Cmpi.rank == 0: test.testT(out, 2)

#out = probe.read(cont=0)
#if Cmpi.rank == 0: test.testT(out, 3)
#Cmpi.convertPyTree2File(out, LOCAL+'/out.cgns')
