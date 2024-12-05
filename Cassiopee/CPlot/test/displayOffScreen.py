# - display (array) -
# display offscreen using osmesa
import CPlot
import Transform as T
import Geom as D

a = D.sphere((0,0,0), 1)

# One image
CPlot.display(a, offscreen=1, bgColor=1, mode=1, meshStyle=2,
              solidStyle=1, posCam=(0,6,0), export="one.png")

# Movie
for i in range(50):
    a = T.rotate(a, (0,0,0), (0,0,1), 1.)
    CPlot.display(a, offscreen=1, bgColor=1, mode=1, meshStyle=2,
                  solidStyle=1, posCam=(0,6,0),
                  exportResolution='680x600', export="export.mpeg")

# free osmesa context (not compulsary)
CPlot.finalizeExport(1)