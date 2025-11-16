# driver: parametric profile
import OCC.Driver as D
import Geom
import Generator
import Converter

# Create a parameter
epaisseur = D.Scalar(12., name='epaisseur')
epaisseur.range = [0, 5, 0.1]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create parameter grid
grid1 = D.Grid(bbox[0:3], bbox[3:], N=(3,3,1), name="grid1")
D.Eq(epaisseur.s, grid1.P[1][2][0].y.s)

# Create parametric profile
spline1 = D.Spline3( grid1, mesh=naca, name='spline1' )

# Create parametric sketch
sketch1 = D.Sketch([spline1], name='sketch1')

# solve for free parameters
D.DRIVER.solve2()
#grid1.P[1][2][0].y.print()

# instantiate a CAD from free parameters
# then mesh and get sensibilities
D.DRIVER.instantiate({'epaisseur': 0.8})
sketch1.writeCAD('out.step')
mesh = sketch1.mesh(0.01, 0.01, 0.01)
D.DRIVER._diff(sketch1, mesh)
Converter.convertArrays2File(mesh, 'dout.plt')

# Build DOE
#D.DRIVER.setDOE({'epaisseur': 0.1})
D.DRIVER.setDOE()
D.DRIVER.createDOE('doe.hdf')
D.DRIVER.walkDOE(sketch1, 0.01, 0.01, 0.01)

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
    mesh = sketch1.mesh(0.01, 0.01, 0.01)
    CPlot.display(mesh)
    time.sleep(0.5)

# build dmesh
mesh1 = D.DRIVER.dmesh(sketch1, mesh, ['epaisseur'], 0.1)
mesh2 = D.DRIVER.dmesh(sketch1, mesh1, ['epaisseur'], 0.1)
Converter.convertArrays2File(mesh+mesh1+mesh2, 'out.plt')
