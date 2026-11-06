# - display360 (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import CPlot.PyTree as CPlot

offscreen = 1

# Create model
a = D.sphere((0,0,0), 1., 50)
C._initVars(a, '{f}={CoordinateZ}')

# Initial view of model
posCam = (0,0,0); posEye = (1,0,0); dirCam = (0,0,1)

CPlot.display360(a, mode="Mesh", posCam=posCam, posEye=posEye, dirCam=dirCam,
                 offscreen=offscreen,
                 export='image360.png', exportResolution='4096x2048')
