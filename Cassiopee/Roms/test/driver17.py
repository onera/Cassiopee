# driver: parametric sphere (with derivatives)
import Roms.Driver as D
import Converter.PyTree as C

# Create a parameter
radius = D.Scalar('radius')
radius.range = [0.1, 10, 0.3]

# Create parametric sphere
sphere1 = D.Sphere('sphere1', (0,0,0), radius, h=(0.1,0.1,0.01))

# solve for free parameters
D.DRIVER.solve()
D.DRIVER.instantiate({'radius': 1.5})

sphere1.writeCAD('out.step')
sphere1.print()

m1 = sphere1.MeshAsReference()
m1[0][0] = 'sphere1'
D.DRIVER.instantiate({'radius': 1.6})
m2 = sphere1.Dmesh()
m2[0][0] = 'sphere2'

#D.DRIVER._dXdmu(sphere1, Mesh=M, freeParams=['radius'], deps=0.1)
C.convertPyTree2File(m1+m2, 'out.cgns')
#D.DRIVER.plot(sphere1)
