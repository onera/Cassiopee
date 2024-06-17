# Adaptation of an unstructured octree mesh
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Geom.PyTree as D
import Converter.Internal as Internal

# Creation of the octree mesh around a NACA0012 profile
naca = D.naca(12)
o = G.octree([naca], [0.01], dfar=10., balancing=1)
o = C.initVars(o, 'Density = {CoordinateX}') 

# Define a model sensor field
import math
K = 1./math.sqrt(0.001*2*math.pi)
K2 = 1./(2*0.001)
def F(x,y):
    if abs(y) < 0.7 and abs(y) > 0.05: return K*math.exp(-(x-0.5)**2*K2)
    else: return 0.
o = C.initVars(o, 'F', F, ['CoordinateX','CoordinateY'])   

# Computation of the indicator onto the octree mesh
npts = Internal.getZoneDim(o)[1]
o2, valInf, valSup = P.computeIndicatorField(o, 'F', nbTargetPts=1.2*npts,
                                             refineFinestLevel=0, bodies=[naca])

# Adaptation of the octree with respect to the indicator
o2 = G.adaptOctree(o2)

# Interpolate the solution on adapted octree
o = C.rmVars(o, ['centers:indicator', 'F'])
o2 = P.extractMesh([o], o2)

# Save file
t = C.newPyTree(['Octree', 'Octree2']); t[2][1][2] += [o]; t[2][2][2] = [o2]
C.convertPyTree2File(t, 'out.cgns')
