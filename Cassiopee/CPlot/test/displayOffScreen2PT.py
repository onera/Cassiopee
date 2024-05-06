# - display (pyTree) -
# display offscreen using GL
import CPlot.PyTree as CPlot
import Geom.PyTree as D

a = D.sphere((0,0,0), 1)
b = D.sphere((1,1,1), 1.2)

# One image, all scene
CPlot.display([a,b], mode='mesh',
              posCam=(0,6,0), posEye=(0,0,0), dirCam=(0,0,1),
              export="one.png", offscreen=2)
CPlot.finalizeExport()

# Compositing in one image
CPlot.display(a, mode='mesh',
              posCam=(0,6,0), posEye=(0,0,0), dirCam=(0,0,1),
              export="composite.png", offscreen=3)
CPlot.finalizeExport(3)
CPlot.display(b, mode='mesh',
              posCam=(0,6,0), posEye=(0,0,0), dirCam=(0,0,1),
              export="composite.png", offscreen=4)
CPlot.finalizeExport(4)

import os; os._exit(0)
