# - Probe mode 3 (pyTree) -
import Post.Probe as Probe
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

# test case
a = G.cartRx((0,0,0), (1,1,1), (30,30,30), (5,5,5))
t = C.newPyTree(['Base', a])
C._initVars(t, '{centers:Fx} = {centers:CoordinateX}')
C._initVars(t, '{centers:Fy} = {centers:CoordinateY}')

b = D.cylinder((50,50,50), 10., 40., N=50, ntype='STRUCT')
tR = C.newPyTree(['Base', b]) # receiver surface tree storing the interpolated flow variable data

tD = Internal.copyRef(t) # donor tree storing the interpolation data (interpolation coefficients)

# create a probe using mode 3
fields = ['centers:Fx', 'centers:Fy']
p3 = Probe.Probe('probe3.cgns', tPermeable=tR, fields=fields, append=False, bufferSize=100)

# compute the interpolation data for probe extraction (mode 3 only)
tD = p3.prepare(tD, loc='centers')

# extract probe information over time
for i in range(110):
    time = 0.1*i
    C._initVars(t, '{centers:Fx} = {centers:Fx} + cos(%f)'%time, isVectorized=True)
    C._initVars(t, '{centers:Fy} = {centers:Fx} + sin(%f)'%time, isVectorized=True)
    for varname in fields: C._cpVars(t, varname, tD, varname)
    tD = C.center2Node(tD, fields) # donor solution is located at the nodes
    p3.extract(tD, time=time)

# force probe flush
p3.flush()