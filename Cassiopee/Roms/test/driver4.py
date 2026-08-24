# driver: parametric profile
import Roms.Driver as D
import Geom
import Generator
import Converter

T1 = D.Part("Part1")

# Create a parameter
epaisseur = T1.Scalar('epaisseur')
epaisseur.range = [0, 5, 0.1]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create parameter grid
grid1 = T1.Grid('grid1', bbox[0:3], bbox[3:], N=(3,3,1))
T1.Eq(epaisseur, grid1.P[1][2][0].y)

# Create parametric profile
spline1 = T1.Spline3('spline1', grid1, mesh=naca)

# Create parametric sketch
sketch1 = T1.Sketch('sketch1', [spline1], h=[0.01,0.01,0.01])

# solve for free parameters
T1.solve()
#grid1.P[1][2][0].y.print()

# instantiate a CAD from free parameters
# then mesh and get sensibilities
T1.instantiate({'epaisseur': 0.8})
sketch1.writeCAD('out.step')
mesh = sketch1.mesh()
T1._dXdmu(sketch1, mesh=mesh, deps=1.e-3)
Converter.convertArrays2File(mesh, 'dout.plt')

# Build DOE
T1.createDOE('doe.hdf')
T1.walkDOE3(sketch1)

# reread one snapshot from DOE file
m = T1.readSnaphot(0)
Converter.convertArrays2File(m, 'reread.plt')

# read snapshots as matrix
F = T1.readAllSnapshots()
# Create a ROM limited to K modes
Phi,S,Vt = T1.createROM(F, K=20)
T1.writeROM('rom.hdf')
# add to file the coordinates of snapshots on POD vectors
T1.addAllCoefs()

# reread and build a snapshot from ROM
coords = T1.readCoefs(0)
m = T1.evalROM(coords)
Converter.convertArrays2File(m, 'reread2.plt')

# instantiate CADs, mesh and display
import CPlot, time
for i in range(20):
    T1.instantiate({'epaisseur': 0.3+i/50.})
    mesh = sketch1.mesh()
    CPlot.display(mesh)
    time.sleep(0.5)
