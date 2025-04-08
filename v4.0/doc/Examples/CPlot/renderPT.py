# - render (pyTree) -
import Generator.PyTree as G
import CPlot.PyTree

a = G.cart((0,0,0), (1,1,1), (10,10,10))
CPlot.PyTree.display([a])
CPlot.PyTree.render()
