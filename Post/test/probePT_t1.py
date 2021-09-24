# - Probe (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G

# test case
a = G.cartRx((0,0,0), (1,1,1), (30,30,30), (5,5,5), depth=0, addCellN=False)
C._initVars(a, '{centers:F} = {centers:CoordinateX}')
t = C.newPyTree(['Base',a])
#C.convertPyTree2File(a, 'out.cgns')

# create a probe
p1 = Probe.Probe(t, (10.,10.,10.), 'probe1.cgns', fields=['centers:F'], append=False)
p1.print()
for i in range(110):
    p1.extract(time=0.1*i)
p1.flush()

# test append
p1 = Probe.Probe(t, (10.,10.,10.), 'probe1.cgns', fields=['centers:F'], append=True)
for i in range(50):
    p1.extract(time=11.+0.1*i)
p1.flush()
