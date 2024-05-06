# - probe (pyTree) -
# Relecture probe mode 2 pour donner un numpy qui concatene toute la solution
# a tous les instants
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import Post.Probe as Probe
import KCore.test as test

LOCAL = test.getLocal()

# Surface probe
s = D.sphere((2.5,2.5,2.5), 1.)
s = T.splitNParts(s, 2)

# Create probe to store surface
probe = Probe.Probe(LOCAL+'/probe.cgns', t=None, fields=['centers:F'], bufferSize=15)

# Extractions
for i in range(40):
    time = 0.1*i
    C._initVars(s, "{centers:F}={centers:CoordinateX}*{centers:CoordinateY}+%g"%time)
    probe.extract(s, time)
probe.flush()

# Relecture par containers
p2 = Probe.Probe(LOCAL+'/probe.cgns', fields=['centers:F'])

out = p2.readCont2(0, 'centers:F', [])
out = p2.readCont2(1, 'centers:F', out)

# out[0] is zone 0, out[0][0] is index 0 (all time)

test.testO(out, 1)
