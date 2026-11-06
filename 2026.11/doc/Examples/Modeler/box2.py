# - cobbleStone (array) -
import Converter as C
import Modeler.Models as Models
import Modeler.WorldZ as World
import CPlot
import random

zmap = World.createZMap(50, 50)
zmap = World.initGround(zmap, type=0)
ground = C.copy(zmap)
allObjects = [ground]

for i in range(30):
    hx = random.random()*5.
    hy = random.random()*5.
    hz = random.random()*1.
    chamfer = -0.2+0.2*random.random()
    a = Models.box((0,0,0),(hx,hy,hz),chamfer=chamfer)
    a = World.placeRandomObject(zmap,a,impact=1,rotateZ=True)
    allObjects.append(a)

pos = World.placeViewerAtCenter(zmap)
CPlot.display(allObjects, displayInfo=0, bgColor=1, shadow=0, posCam=pos[0],
              posEye=pos[1], meshStyle=3)
World.simpleLoop2(zmap, pos)
#C.convertArrays2File(allObjects, 'out.plt')
