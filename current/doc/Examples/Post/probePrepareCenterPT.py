# - Probe prepare - centers (pyTree) -
import Post.Probe as Probe
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

# test case
a = G.cart((0,0,0), (0.1,0.1,0.1), (11,11,11))
t = C.newPyTree(['Base', a])
C._initVars(t, '{centers:F} = {centers:CoordinateX}')

b = D.box((0.2,0.2,0.2), (0.7,0.7,0.7), 6, ntype='STRUCT')
tR = C.newPyTree(['Base', b]) # receiver surface tree storing the interpolated flow variable data

tD = Internal.copyRef(t) # donor tree storing the interpolation data (interpolation coefficients)

# create a probe using mode 3
p3 = Probe.Probe('probe3.cgns', tPermeable=tR, fields=['centers:F'], append=False, bufferSize=100)

# compute the interpolation data for probe extraction (mode 3 only)
C._initVars(tD, '{cellN} = 1') # donor solution is located at the nodes
C._initVars(tR, '{centers:cellN} = 2')
tD = p3.prepare(tD, loc='centers')