# driver: parametric sphere (with derivatives)
import Roms.Driver as D
import Converter.PyTree as C

# Create a parameter
radius = D.Scalar('radius')
radius.range = [0.1, 10, 0.3]

# Create parametric sphere
sphere1 = D.Sphere('sphere1', (0,0,0), radius, h=(0.1,0.1,0.01))

# solve and instantiate
D.DRIVER.solve()
D.DRIVER.instantiate({'radius': 1.5})

sphere1.writeCAD('out.step')
sphere1.print()

# compute dXdmu
m = sphere1.MeshAsReference()
D.DRIVER._dXdmu(sphere1, Mesh=m, freeParams=['radius'], deps=1.e-6)
C.convertPyTree2File(m, 'out.cgns')

D.DRIVER.plot(sphere1)
