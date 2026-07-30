# driver: parametric profile in parallel
import Roms.Driver as D
import Geom
import Generator

T1 = D.Part("Part1")

# Create a parameter
epaisseur = T1.Scalar('epaisseur')
epaisseur.range = [10, 15]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create parameter grid
grid1 = T1.Grid('grid1', bbox[0:3], bbox[3:], N=(3,3,1))
grid1.P[1][2][0].y.range = [0, 5, 0.1]
T1.Eq(epaisseur, grid1.P[1][2][0].y)

# Create parametric profile
spline1 = T1.Spline3('spline1', grid1, mesh=naca)

# Create parametric sketch
sketch1 = T1.Sketch('sketch1', [spline1], h=[0.01,0.01,0.01])

# solve for free parameters
T1.solve()

# Build DOE
T1.createDOE('doe.hdf')
T1.walkDOE3(sketch1)

# read snapshots as matrix
#F = T1.readAllSnapshots()
# Create a ROM limited to K modes
#Phi,S,Vt = T1.createROM(F, K=20)
#T1.writeROM('rom.hdf')
# add to file the coordinates of snapshots on POD vectors
#T1.addAllCoefs()