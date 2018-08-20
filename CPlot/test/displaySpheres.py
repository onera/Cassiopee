# - display (pyTree) -
# Affichage du shader sphere
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import time

a = D.sphere((0,0,0),1.,N=5)
a = C.convertArray2Hexa(a)

for i in xrange(21):
    print i*1./10.
    CPlot._addRender2Zone(a, material='Sphere', color='White',
                          shaderParameters=[1.5,i*1./10.])

    CPlot.display(a, displayBB=0, mode='render', bgColor=1)
    time.sleep(1.)
