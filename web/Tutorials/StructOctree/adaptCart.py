# Adaptation of an structured octree mesh
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Geom.PyTree as D
import Converter.Internal as Internal

# Creation of an unstructured octree mesh around a NACA0012 profile
naca = D.naca(12)
o = G.octree([naca], [0.01], dfar=10., balancing=1)

# Creation of the structured Cartesian set of grids 
res0 = G.octree2Struct(o, vmin=5, ext=0)

# Fill the undefined BCs: external borders are defined by farfield BCs here
res0 = C.fillEmptyBCWith(res0, 'nref', 'BCFarfield', dim=2)

res0 = C.initVars(res0, 'Density={CoordinateX}')

# Define a model sensor field
import math
K = 1./math.sqrt(0.001*2*math.pi); K2 = 1./(2*0.001)
L = 1./math.sqrt(0.01*2*math.pi); L2 = 1./(2*0.01)
def F(x,y):
    if abs(y) < 0.7 and abs(y) > 0.05: return K*math.exp(-K2*(x-0.5)**2)
    else: return 0.
res0 = C.initVars(res0, 'F', F, ['CoordinateX','CoordinateY'])   

# Project the sensor field F on the unstructured octree:
o2 = P.computeIndicatorValue(o, res0, 'F')

# Computation of the indicator onto the unstructured octree mesh
npts = Internal.getZoneDim(o)[1]
o2, valInf, valSup = P.computeIndicatorField(o2, 'F', nbTargetPts=1.1*npts,
                                             refineFinestLevel=0, bodies=[naca])

# Adaptation of the octree wrt the indicator
o2 = G.adaptOctree(o2)

# Creation of the adapted structured Cartesian grids
res = G.octree2Struct(o2, vmin=5, ext=0)

# Undefined BC: external borders are defined by farfield BCs
res = C.fillEmptyBCWith(res, 'nref', 'BCFarfield', dim=2)

# Interpolate solution on adapted set of Cartesian grids
res0 = C.rmVars(res0, ['centers:indicator', 'F'])
res = P.extractMesh(res0, res)

# Save file
t = C.newPyTree(['CARTESIAN0','CARTESIAN']); t[2][1][2] += res0; t[2][2][2] += res
C.convertPyTree2File(t, 'out.cgns')
