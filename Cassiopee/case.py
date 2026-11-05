# create a cube test case in h5 format, variable size
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import numpy

# size in M points
size = 100 # 1 10 30 50 100

npts = (size*1.e6) ** (1./3.)
npts = int(npts)
npts = max(npts, 2)
print("Generating %s M points (%d)."%(size, npts))

# create grid
a = G.cartHexa((0,0,0), (10./(npts-1),10./(npts-1),10./(npts-1)), (npts,npts,npts))

# add fields
def initDensity(x,y,z): return 1.0 + 0.03*x + 0.01*y - 0.02*z
def initMomentumX(x,y,z): return 0.1 + 0.05*x - 0.05*y + 0.02*z
def initMomentumY(x,y,z): return     - 0.03*x + 0.05*y - 0.02*z
def initMomentumZ(x,y,z): return 0.1 + 0.01*x + 0.03*y + 0.02*z
def initEnergy(x,y,z): return 15. + 0.02*x + 0.03*y + 0.01*z
def initNut(x,y,z): return 500.*numpy.exp(-5*x**2)
C._initVars(a, 'centers:Density', initDensity, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)
C._initVars(a, 'centers:MomentumX', initMomentumX, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)
C._initVars(a, 'centers:MomentumY', initMomentumY, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)
C._initVars(a, 'centers:MomentumZ', initMomentumZ, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)
C._initVars(a, 'centers:EnergyStagnationDensity', initEnergy, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)
C._initVars(a, 'centers:Viscosity_EddyMolecularRatio', initNut, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)

# add bc
p = P.exteriorFaces(a)
C._addBC2Zone(a, 'far', 'BCFarfield', subzone=p)

# save h5
C.convertPyTree2File(a, 'case%sM.h5'%str(size))
