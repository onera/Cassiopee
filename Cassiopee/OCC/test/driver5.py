# driver: parametric circle (with derivatives)
import OCC.Driver as D
import Converter

# Create a parameter
radius = D.Scalar(1., name='radius')
radius.range = [0.1, 10]

# Create parametric circle
circle1 = D.Circle( (0,0,0), radius, name='circle1' )

# Create parametric sketch
sketch1 = D.Sketch([circle1], name='sketch1')

# solve for free parameters
D.DRIVER.solve2()

# instantiate a CAD
D.DRIVER.instantiate({'radius': 1.5})
sketch1.writeCAD('out.step')
mesh = sketch1.mesh(0.01, 0.01, 0.01)
D.DRIVER._diff(sketch1, mesh)
Converter.convertArrays2File(mesh, 'out.plt')

# Build DOE
D.DRIVER.setDOE({'radius': 0.3})
D.DRIVER.createDOE('doe.hdf')
D.DRIVER.walkDOE(sketch1, 0.01, 0.01, 0.01)

# reread one snaphsot from DOE file
m = D.DRIVER.readSnaphot(0)
Converter.convertArrays2File(m, 'reread.plt')

# read snapshots as matrix
F = D.DRIVER.readAllSnapshots()
D.DRIVER.createROM(F, K=-1)
D.DRIVER.writeROM('rom.hdf')
# add to file the coordinates of snapshots on POD vectors
D.DRIVER.addAllCoefs()

# reread and build a snapshot from ROM
coords = D.DRIVER.readCoefs(0)
m = D.DRIVER.evalROM(coords)
Converter.convertArrays2File(m, 'reread2.plt')

# instantiate CADs, mesh and display
import CPlot, time
for i in range(50):
    D.DRIVER.instantiate({'radius': 10-i/10.})
    mesh = sketch1.mesh(0.01, 0.01, 0.01)
    CPlot.display(mesh)
    time.sleep(0.5)
