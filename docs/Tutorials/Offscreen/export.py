# Offscreen rendering

import Generator as G
import Converter as C
import CPlot
import time
import Transform as T

a = G.cart( (0,0,0), (1,1,1), (5,5,5) )

# dump a PNG 
CPlot.display([a], mode=0, posCam=(10,2,2), posEye=(2,2,2),
              export='export.png',
              exportResolution='600x600', offscreen=2) ; CPlot.finalizeExport()

# dump a MPEG movie
for i in range(200):
    a = T.rotate(a, (2,2,2), (0,0,1), 1) ; time.sleep(0.1)
    CPlot.display([a], mode=1, export='export.mpeg',
                  exportResolution='600x600', offscreen=2)
time.sleep(2) ; CPlot.finalizeExport(action=1)
import os ; os._exit(1)
