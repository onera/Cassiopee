# driver: parametric circle (with derivatives)
import Roms.Driver as D
import Converter, Geom

# Create a parameter
radius = D.Scalar('radius')
radius.range = [0.1, 10, 0.3]

# Create parametric circle
circle1 = D.Circle('circle1', (0,0,0), radius)

# Create parametric sketch
sketch1 = D.Sketch('sketch1', [circle1], h=[0.01,0.01,0.01])

# solve for free parameters
D.DRIVER.solve()

# compute dL/dR by FD
deps = 1.e-10
D.DRIVER.instantiate({'radius': 1.5})
m1 = sketch1.mesh()
L1 = Geom.getLength(m1)

D.DRIVER.instantiate({'radius': 1.5+deps})
m2 = sketch1.mesh()
L2 = Geom.getLength(m2)

dLdR = (L2-L1)/1.e-10
print("dLdR by FD:", dLdR)

# compute dL/dR by dL/dX * dX/dR by FD and AD
D.DRIVER.instantiate({'radius': 1.5})
D.DRIVER._dXdmu(sketch1, m1, freeParams=['radius'])
print(m1)