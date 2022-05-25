# - Probe (pyTree) -
import Post.Probe as Probe
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Internal as Internal

# test case
a = G.cartRx((0,0,0), (1,1,1), (30,30,30), (5,5,5), depth=0, addCellN=False)
C._initVars(a, '{centers:F} = {centers:CoordinateX}')
t = C.newPyTree(['Base', a])

# create a probe
p1 = Probe.Probe('probe1.cgns', t, X=(10.,10.,10.), fields=['centers:F'], append=False)
p1.printInfo()
for i in range(110):
    time = 0.1*i
    C._initVars(t, f'{{centers:F}} = {{centers:CoordinateX}}+10.*sin({time})')
    p1.extract(time=time)
p1.flush()

# reread probe from file
p1 = Probe.Probe('probe1.cgns')
out = p1.read()
C.convertPyTree2File(out, 'out.cgns')
