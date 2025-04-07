# - display (array) -
# display offscreen using GL
import CPlot
import Transform as T
import Geom as D

a = D.sphere((0,0,0), 1, N=200)

# Multi images
#CPlot.displayFBO(a, bgColor=1, mode=0, meshStyle=2,
#              solidStyle=1, posCam=(0,6,0), export='one.png')

CPlot.displayFBO(a, mode=0, offscreen=2)
