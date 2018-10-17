# - display (array) -
# Offscreen using mesa
import CPlot
import Transform as T
import Converter as C
import Geom as D
import time

a = D.sphere( (0,0,0), 1)

# One image
CPlot.display([a], offscreen=1, bgColor=1, mode=1, meshStyle=2,
              solidStyle=1, posCam=(0,6,0), export="one.png")

# Movie
for i in xrange(50):
    a = T.rotate(a, (0,0,0), (0,0,1), 1.)
    CPlot.display([a], offscreen=1, bgColor=1, mode=1, meshStyle=2,
                  solidStyle=1, posCam=(0,6,0), 
                  exportResolution='680x600', export="export.mpeg")
time.sleep(0.1)
time.sleep(1); CPlot.finalizeExport()
import os; os._exit(0)
