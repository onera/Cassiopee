# driver: parametric profile in parallel
import Roms.Driver as D
import Geom
import Generator

# Create a parameter
epaisseur = D.Scalar('epaisseur', 12.)
epaisseur.range = [10,15]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create parameter grid
grid1 = D.Grid('grid1', bbox[0:3], bbox[3:], N=(3,3,1))
grid1.P[1][2][0].y.range = [0, 5, 0.1]
D.Eq(epaisseur.s, grid1.P[1][2][0].x.s)

# Create parametric profile
spline1 = D.Spline3('spline1', grid1, mesh=naca)

# Create parametric sketch
sketch1 = D.Sketch('sketch1', [spline1])

# solve for free parameters
D.DRIVER.solve2()

# Build DOE
D.DRIVER.createDOE('doe.hdf')
D.DRIVER.walkDOE3(sketch1, 0.01, 0.01, 0.01)

# read snapshots as matrix
#F = D.DRIVER.readAllSnapshots()
# Create a ROM limited to K modes
#Phi,S,Vt = D.DRIVER.createROM(F, K=20)
#D.DRIVER.writeROM('rom.hdf')
# add to file the coordinates of snapshots on POD vectors
#D.DRIVER.addAllCoefs()