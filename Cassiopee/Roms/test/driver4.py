# driver: parametric profile
import Roms.Driver as D
import Geom
import Generator
import Converter

# Create a parameter
epaisseur = D.Scalar('epaisseur')
epaisseur.range = [0, 5, 0.1]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create parameter grid
grid1 = D.Grid('grid1', bbox[0:3], bbox[3:], N=(3,3,1))
D.Eq(epaisseur, grid1.P[1][2][0].y)

# Create parametric profile
spline1 = D.Spline3('spline1', grid1, mesh=naca)

# Create parametric sketch
sketch1 = D.Sketch('sketch1', [spline1], h=[0.01,0.01,0.01])

# solve for free parameters
D.DRIVER.solve()
#grid1.P[1][2][0].y.print()

# instantiate a CAD from free parameters
# then mesh and get sensibilities
D.DRIVER.instantiate({'epaisseur': 0.8})
sketch1.writeCAD('out.step')
mesh = sketch1.mesh()
D.DRIVER._dXdmu(sketch1, mesh)
Converter.convertArrays2File(mesh, 'dout.plt')

# Build DOE
D.DRIVER.createDOE('doe.hdf')
D.DRIVER.walkDOE3(sketch1)

# reread one snapshot from DOE file
m = D.DRIVER.readSnaphot(0)
Converter.convertArrays2File(m, 'reread.plt')

# read snapshots as matrix
F = D.DRIVER.readAllSnapshots()
# Create a ROM limited to K modes
Phi,S,Vt = D.DRIVER.createROM(F, K=20)
D.DRIVER.writeROM('rom.hdf')
# add to file the coordinates of snapshots on POD vectors
D.DRIVER.addAllCoefs()

# reread and build a snapshot from ROM
coords = D.DRIVER.readCoefs(0)
m = D.DRIVER.evalROM(coords)
Converter.convertArrays2File(m, 'reread2.plt')

# instantiate CADs, mesh and display
import CPlot, time
for i in range(20):
    D.DRIVER.instantiate({'epaisseur': 0.3+i/50.})
    mesh = sketch1.mesh()
    CPlot.display(mesh)
    time.sleep(0.5)