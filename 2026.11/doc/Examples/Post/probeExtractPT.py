# - Probe extract (pyTree) -
import Post.Probe as Probe
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

# test case
a = G.cart((0,0,0), (0.1,0.1,0.1), (51,11,11))
t = C.newPyTree(['Base', a])
C._initVars(t, '{Density} = 1.')
C._initVars(t, '{VelocityX} = {CoordinateX}-1')
C._initVars(t, '{VelocityY} = {CoordinateY}-1')
C._initVars(t, '{VelocityZ} = {CoordinateZ}-1')

b1 = G.cart((0,0,0), (0.,0.1,0.1), (1,11,11)); b1[0] = 'upstream'
b2 = G.cart((5,0,0), (0.,0.1,0.1), (1,11,11)); b2[0] = 'downstream'
tR = C.newPyTree(['Base', [b1,b2]])
G._getSmoothNormalMap(tR)

tD = Internal.copyRef(t)

# create a first probe using mode 3
fields = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ']
p3 = Probe.Probe('probe3.cgns', tPermeable=tR, fields=fields, append=False, bufferSize=100)
tD = p3.prepare(tD, loc='nodes')

# create a second probe using mode 4
p4 = Probe.Probe('probe4.cgns', fields=['massflow_upstream', 'massflow_downstream'], bufferSize=100, append=False)

# extract flow information on the receiver surface tree
for varname in fields: C._cpVars(t, varname , tD, varname)
p3.extract(tD, onlyTransfer=True)

# tR interpolated data are now updated: integral quantities can be computed on each surface
C._initVars(tR, '{massflow}={Density}*({sx}*{VelocityX}+{sy}*{VelocityY}+{sz}*{VelocityZ})')

z_upstream = Internal.getNodeFromName(tR, 'upstream')
massflow_upstream = P.integ(z_upstream, 'massflow')[0]

z_downstream = Internal.getNodeFromName(tR, 'downstream')
massflow_downstream = P.integ(z_downstream, 'massflow')[0]

# store data information in the p4 probe using single integrated values
p4.extract(time=1, value=[massflow_upstream, massflow_downstream])

p4.flush()