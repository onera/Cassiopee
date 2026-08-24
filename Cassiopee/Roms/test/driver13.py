# driver: parametric boolean + dmesh
import Roms.Driver as D
import Converter.PyTree as C

T1 = D.Part("Part1")

# Create parameter
length = T1.Scalar('length')
length.range = [3., 5., 0.1]

# Create sketch 1
circle1 = T1.Circle('circle1', (0,0,0), 1.)
sketch1 = T1.Sketch('sketch1', [circle1])

# Create sketch 2
circle2 = T1.Circle('circle2', (0,0,5), 1.)
sketch2 = T1.Sketch('sketch2', [circle2])

# surface1
surface1 = T1.Loft('surface1', [sketch1, sketch2])

# Create sketch 3
circle3 = T1.Circle('circle3', (0,0,0), 0.5)
sketch3 = T1.Sketch('sketch3', [circle3])

# Create sketch 4
circle4 = T1.Circle('circle4', (0,0,1), 0.5)
sketch4 = T1.Sketch('sketch4', [circle4])
T1.Eq(circle4.P[0].z, length)

# surface2
surface2 = T1.Loft('surface2', [sketch3, sketch4])
surface2.rotAxis.x.v = 0
surface2.rotAxis.y.v = 1
surface2.rotAxis.z.v = 0
surface2.rotAngle.v = 90.
surface2.position.z.v = 2.

# surface finale
#surface = T1.Merge('surface', listSurfaces=[surface1,surface2])
surface = T1.Union('surface', listSurfaces1=[surface1], listSurfaces2=[surface2], h=[0.1,0.1,0.1])

# test
T1.solve()
T1.instantiate({'length': 4.})
#surface.writeCAD('out.step')
m = surface.MeshAsReference()
#C.convertPyTree2File(m, 'out.cgns')

import CPlot.PyTree as CPlot, time
point = T1.walkDOE()
i = 0
while point is not None:
    T1.instantiate(point)
    #mesh = surface.Mesh()
    m = surface.Dmesh()
    C.convertPyTree2File(m, 'out%02d.cgns'%i)
    CPlot.display(m)
    point = T1.walkDOE()
    time.sleep(0.5)
    i += 1
print("done", flush=True)