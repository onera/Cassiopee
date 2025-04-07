# - display (export) -
import Generator as G
import CPlot
import time
import Transform as T

a = G.cart( (0,0,0), (1,1,1), (5,5,5) )

# Bitmap PPM format
CPlot.display(a, mode=0, export='export.ppm'); time.sleep(1)
# Bitmap postscript
CPlot.display(a, mode=0, export='export.ps'); time.sleep(1)
# PNG
CPlot.display(a, mode=0, export='export.png'); time.sleep(1)
# MPEG
for i in range(200):
    a = T.rotate(a, (2,2,2), (0,0,1), 1); time.sleep(0.1)
    CPlot.display(a, mode=0, meshStyle=2, export='export.mpeg')
time.sleep(2); CPlot.finalizeExport()
