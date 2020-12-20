# - display (pyTree) -
# Offscreen display using osmesa
import CPlot.PyTree
import Geom as D

a = D.sphere((0,0,0), 1)
b = D.sphere((1,1,1), 1.2)

# One image, all scene
CPlot.display([a,b], mode='mesh',
              posCam=(0,6,0), posEye=(0,0,0), dirCam=(0,0,1),
              export="one.png", offscreen=1)

# Compositing in one image
CPlot.display(a, mode='mesh',
              posCam=(0,6,0), posEye=(0,0,0), dirCam=(0,0,1),
              export="composite.png", offscreen=5)
CPlot.display(b, mode='mesh',
              posCam=(0,6,0), posEye=(0,0,0), dirCam=(0,0,1),
              export="composite.png", offscreen=6)
