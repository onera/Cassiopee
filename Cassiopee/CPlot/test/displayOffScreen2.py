# - display (array) -
# display offscreen using OpenGL
import CPlot
import Transform as T
import Geom as D

a = D.sphere((0,0,0), 1, N=200)

# Multi images
CPlot.display(a, offscreen=2, bgColor=1, mode=0, meshStyle=2,
              solidStyle=1, posCam=(0,6,0), export='one.png')
CPlot.finalizeExport(2) # wait for end of file write
for i in range(5):
    a = T.rotate(a, (0,0,0), (0,0,1), 1.)
    CPlot.display(a, offscreen=2, bgColor=1, mode=0, meshStyle=2,
                  solidStyle=1, posCam=(0,6,0), export='one%d.png'%i)
    CPlot.finalizeExport(2) # wait for end of file write
import os; os._exit(0)

# Mpeg Movie
for i in range(50):
    a = T.rotate(a, (0,0,0), (0,0,1), 1.)
    CPlot.display(a, bgColor=1, mode=0,
                  solidStyle=1, posCam=(0,6,0), export='export.mpeg',
                  exportResolution='680x600', offscreen=2)
    CPlot.finalizeExport(2) # wait for end of file write
CPlot.finalizeExport(1) # close mpeg file
import os; os._exit(0)
