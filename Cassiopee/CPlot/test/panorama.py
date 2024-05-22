# - panorama (pyTree) --
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import KCore.Vector as Vector
import Converter.Internal as Internal

offscreen = 2

# Create model
a = G.cart((-4,-4,-4), (1,1,1), (9,9,9))
#a = D.sphere((0,0,0), 1., 50)
C._initVars(a, '{f}={CoordinateZ}')

# Initial view of model
posCam = (0,0,0); posEye = (1,0,0); dirCam = (0,0,1)

# Create cube map views

# Compute all view vectors
vz = Vector.normalize(dirCam)
v1 = Vector.sub(posEye, posCam)
v2 = Vector.cross(vz, v1)
n = Vector.norm(v1)
v3 = Vector.mul(n, vz)

# right
posEye0 = Vector.sub(posCam, v2); dirCam0 = dirCam
CPlot.display(a, mode='Mesh', isoEdges=1,
              posCam=posCam, posEye=posEye0, dirCam=dirCam0, viewAngle=90.,
              export='cube_right.png', offscreen=offscreen, exportResolution='800x800')
CPlot.finalizeExport(offscreen)

# left
posEye0 = Vector.add(posCam, v2); dirCam0 = dirCam
CPlot.display(a, isoEdges=1,
              posCam=posCam, posEye=posEye0, dirCam=dirCam0, viewAngle=90.,
              export='cube_left.png', offscreen=offscreen)
CPlot.finalizeExport(offscreen)

# front
posEye0 = posEye; dirCam0 = dirCam
CPlot.display(a, isoEdges=1,
              posCam=posCam, posEye=posEye0, dirCam=dirCam0, viewAngle=90.,
              export='cube_front.png', offscreen=offscreen)
CPlot.finalizeExport(offscreen)

# back
posEye0 = Vector.sub(posCam, v1); dirCam0 = dirCam
CPlot.display(a, isoEdges=1,
              posCam=posCam, posEye=posEye0, dirCam=dirCam0, viewAngle=90.,
              export='cube_back.png', offscreen=offscreen)
CPlot.finalizeExport(offscreen)

# top
posEye0 = Vector.add(posCam, v3); dirCam0 = Vector.mul(-1, v1)
CPlot.display(a, isoEdges=1,
              posCam=posCam, posEye=posEye0, dirCam=dirCam0, viewAngle=90.,
              export='cube_top.png', offscreen=offscreen)
CPlot.finalizeExport(offscreen)

# bot
posEye0 = Vector.sub(posCam, v3); dirCam0 = Vector.mul(+1, v1)
CPlot.display(a, isoEdges=1,
              posCam=posCam, posEye=posEye0, dirCam=dirCam0, viewAngle=90.,
              export='cube_bottom.png', offscreen=offscreen)
CPlot.finalizeExport(offscreen)

# Create 360 image
t = C.newPyTree(['Base'])
CPlot.display(t, panorama=1,
              offscreen=offscreen, export='image360.png', exportResolution="3200x1600")
CPlot.finalizeExport(offscreen)

# Create a sphere for skybox rendering (optional)
a = D.sphere((0,0,0), 1., N=30)

# build _u_ _v_
C._initVars(a, '_u_=0')
C._initVars(a, '_v_=0')
pu = Internal.getNodeFromName(a, '_u_')[1]
pv = Internal.getNodeFromName(a, '_v_')[1]
(ni,nj) = pu.shape
for j in range(nj):
    for i in range(ni):
        pu[i,j] = j*1./(nj-1)
        pv[i,j] = i*1./(ni-1)

# add render tags
CPlot._addRender2Zone(a, material='Texmat')

t = C.newPyTree(['Base', a])
CPlot._addRender2PyTree(t, materials=['image360.png'] )

C.convertPyTree2File(t, 'out.cgns')