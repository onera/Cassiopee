# - display (array) -
# display composite using GL
import Geom      as D
import Generator as G
import CPlot
import os

sphere = D.sphere((0.,0.,0.),1.0)
cube   = G.cart((0.2,0.2,0.2),(1.,1.,1.),(2,2,2))
cube2  = G.cart((-1.5,-0.2,-1.2),(1.,1.,1.),(2,2,2))

# Full display
CPlot.display([sphere,cube2], posCam=(-0.449,1.63,-2.9), posEye=(0.35,0.35,0.35), dirCam=(0.94,-0.16,-0.3), export="ref.png", offscreen=2)
CPlot.finalizeExport()

# Accumulate in images
CPlot.display([sphere], posCam=(-0.449,1.63,-2.9), posEye=(0.35,0.35,0.35), dirCam=(0.94,-0.16,-0.3), export="composited.png", offscreen=3)
CPlot.finalizeExport(3)

CPlot.display([cube2], posCam=(-0.449,1.63,-2.9), posEye=(0.35,0.35,0.35),dirCam=(0.94,-0.16,-0.3), export="composited.png", offscreen=4)
CPlot.finalizeExport(4)

os._exit(0)
