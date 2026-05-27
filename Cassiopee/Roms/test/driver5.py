# driver: parametric circle (with derivatives)
import Roms.Driver as D
import Converter, Geom, numpy

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
print("dLdR by FD:", dLdR, 2*numpy.pi)

# compute dX/dmu by FD
D.DRIVER.instantiate({'radius': 1.5})
D.DRIVER._dXdmu(sketch1, mesh=m1, freeParams=['radius'])

# compute  dD/dmu by FD
dDdmu = D.DRIVER.dDdmu(sketch1, mesh=m1, freeParams=['radius'], deps=1.e-6)
print("dDdR by FD:", dDdmu[0], 2*numpy.pi*1.5)

Converter.convertArrays2File(m1, 'out.plt')
