# - Probe mode 1 (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Generator.PyTree as G

# test case
a = G.cartRx((0,0,0), (1,1,1), (30,30,30), (5,5,5))
t = C.newPyTree(['Base', a])
C._initVars(t, '{centers:Fx} = {centers:CoordinateX}')
C._initVars(t, '{centers:Fy} = {centers:CoordinateY}')

# create a probe using mode 1
fields = ['centers:Fx', 'centers:Fy']
p1 = Probe.Probe('probe1.cgns', t, ind=(10,10,10), blockName='cart1-1-1', fields=fields, append=False, bufferSize=100)
p1.printInfo()

# extract probe information over time
for i in range(110):
    time = 0.1*i
    C._initVars(t, '{centers:Fx} = {centers:Fx} + cos(%f)'%time, isVectorized=True)
    C._initVars(t, '{centers:Fy} = {centers:Fx} + sin(%f)'%time, isVectorized=True)
    p1.extract(t, time=time)

# force probe flush
p1.flush()