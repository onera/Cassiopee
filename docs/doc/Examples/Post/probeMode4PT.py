# - Probe mode 4 (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

import numpy

# test case
a = G.cartRx((0,0,0), (1,1,1), (30,30,30), (5,5,5))
t = C.newPyTree(['Base', a])
C._initVars(t, '{centers:Fx} = {centers:CoordinateX}')
C._initVars(t, '{centers:Fy} = {centers:CoordinateY}')

# create a probe using mode 4
p4 = Probe.Probe('probe4.cgns', fields=['minFx', 'maxFx', 'minFy', 'maxFy'], append=False, bufferSize=100)
p4.printInfo()

# extract probe information over time
for i in range(110):
    time = 0.1*i
    maxFx, minFx = 0., 1.e6
    maxFy, minFy = 0., 1.e6
    C._initVars(t, '{centers:Fx} = {centers:Fx} + cos(%f)'%time, isVectorized=True)
    C._initVars(t, '{centers:Fy} = {centers:Fx} + sin(%f)'%time, isVectorized=True)
    for z in Internal.getZones(t):
        Fx = Internal.getNodeFromName(z, 'Fx')[1]
        Fy = Internal.getNodeFromName(z, 'Fy')[1]

        minFx = min(minFx, numpy.min(Fx))
        maxFx = max(maxFx, numpy.max(Fx))
        
        minFy = min(minFy, numpy.min(Fy))
        maxFy = max(maxFy, numpy.max(Fy))
        
    p4.extract(t, time=time, value=[minFx, maxFx, minFy, maxFy])

# force probe flush
p4.flush()