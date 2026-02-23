# - Probe read (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Converter.Internal as Internal
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
for i in range(210):
    time = 0.1*i
    C._initVars(t, '{centers:Fx} = {centers:Fx} + cos(%f)'%time, isVectorized=True)
    C._initVars(t, '{centers:Fy} = {centers:Fx} + sin(%f)'%time, isVectorized=True)
    p1.extract(t, time=time)

# flush current buffered data to disk to save last extractions
p1.flush()

# reread probe from file (not necessary here as the probe object p1 is already loaded)
p1 = Probe.Probe('probe1.cgns')

# extract all points
out = p1.read(cont=0) # get all the information located in the first container
Internal.printTree(out) # list of 100 zones named 'probe@it' where 'it' goes from 0 to 99 (bufferSize-1)

# extract all times
out = p1.read(ind=0, probeName=0) # get all the information over time of the first point of the first probe zone
Internal.printTree(out) # mono zone with numpy arrays of size 210
