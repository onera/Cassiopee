import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as I
import XCore.PyTree as X
import XCore.xcore
import numpy as np

# Tag cells between two ellipses
def F(x, y):
    a = 0.3
    b = 0.25
    dx2 = ((x-0.5)/a)**2
    dy2 = ((y-0.5)/b)**2
    sr = 1.
    br = 1.2
    return float(dx2 + dy2 > sr and dx2 + dy2 < br)

# Make initial mesh
nc = 50
t = G.cartHexa((0,0,0),(1./nc,1./nc,1./nc),(nc+1,nc+1,2))
t = C.convertArray2NGon(t)
t = G.close(t)
print("Initial cells:", C.getNCells(t))

# Make adaptMesh
h = X.createAdaptMesh(t)

# Freeze adaptation in this direction (normal vector)
fv = np.array([0., 0., 1.])

# Adaptation loop
niter = 3
for i in range(niter):
    # Extract refined mesh
    l = X.extractLeafMesh(h)

    # Compute custom sensor on the refined mesh
    C._initVars(l, 'centers:alpha', F, ['centers:CoordinateX', 'centers:CoordinateY'])
    fsolc = I.getNodeFromName(l, I.__FlowSolutionCenters__)
    fld = I.getNodeFromName(l, 'alpha')[1]
    
    # Adapt
    X._adaptMeshSeq(h, [fld], fv)

# Output
l = X.extractLeafMesh(h)
print('Leaves:', C.getNCells(l))
C.convertPyTree2File(l, 'refined.cgns')
