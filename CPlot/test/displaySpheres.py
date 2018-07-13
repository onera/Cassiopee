# - display (pyTree) -
# Affichage du shader spheres (mode solid)
import Generator.PyTree as G
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Converter.PyTree as C

a = D.sphere((0,0,0),1.)
a = C.convertArray2Hexa(a)
a = CPlot.addRender2Zone(a, material='Spheres', color='White')
t = C.newPyTree(['Base',a])

CPlot.display(t, displayBB=0, mode=2)
