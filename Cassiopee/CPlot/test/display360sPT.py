# - display360 (pyTree) -
# stereo version
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import CPlot.PyTree as CPlot

offscreen = 1

# Create model
#a = G.cart((-4,-4,-4), (1,1,1), (9,9,9))
a = D.sphere((0,0,0), 1., 50)
C._initVars(a, '{f}={CoordinateZ}')

# Initial view of model
posCam = (0,0,0); posEye = (1,0,0); dirCam = (0,0,1)

CPlot.display360(a, mode="Mesh", posCam=posCam, posEye=posEye, dirCam=dirCam,
                 offscreen=offscreen, type360=0,
                 stereo=1, stereoDist=0.1,
                 export='image360.png', exportResolution='4096x2048')
