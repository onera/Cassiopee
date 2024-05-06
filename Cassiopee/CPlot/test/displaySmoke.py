# - display (pyTree) -
# Affichage du shader smoke (mode solid)
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Converter.PyTree as C

a = D.sphere((0,0,0),1.)
a = C.convertArray2Hexa(a)
a = CPlot.addRender2Zone(a, material='Smoke', color='White')
t = C.newPyTree(['Base',a])

CPlot.display(t, mode=2)
