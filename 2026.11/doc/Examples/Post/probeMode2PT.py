# - Probe mode 2 (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D

# test case
b = D.cylinder((50,50,50), 10., 40., N=50, ntype='STRUCT')
tR = C.newPyTree(['Base', b])
C._initVars(tR, '{centers:Fx} = {centers:CoordinateX}')
C._initVars(tR, '{centers:Fy} = {centers:CoordinateY}')

# create a probe using mode 2
fields = ['centers:Fx', 'centers:Fy']
p2 = Probe.Probe('probe2.cgns', tR, fields=fields, append=False, bufferSize=100)
p2.printInfo()

# extract probe information over time
for i in range(110):
    time = 0.1*i
    T._rotate(tR, (50.,50.,50.), (0.,0.,1.), 0.1)
    C._initVars(tR, '{centers:Fx} = {centers:CoordinateX}')
    C._initVars(tR, '{centers:Fy} = {centers:CoordinateY}')
    p2.extract(tR, time=time)

# force probe flush
p2.flush()