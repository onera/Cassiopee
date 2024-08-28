# - render (array) -
import Generator as G
import CPlot

a = G.cart((0,0,0), (1,1,1), (10,10,10))
CPlot.display(a)
CPlot.render()
