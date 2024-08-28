# - cobbleStone (array) -
import Converter as C
import Modeler.Models as Models
import Modeler.WorldZ as World
import CPlot

zmap = World.createZMap(50, 50)
zmap = World.initGround(zmap, type=0)
ground = C.copy(zmap)
allObjects = [ground]

for i in range(800):
    a = Models.cobbleStone(hx=1., hy=0.3, hz=0.2)
    a = World.placeRandomObject(zmap,a,impact=0,rotateZ=True)
    allObjects.append(a)

pos = World.placeViewerAtCenter(zmap)
CPlot.display(allObjects, displayInfo=0, bgColor=1, shadow=0, posCam=pos[0],
              posEye=pos[1], meshStyle=3)
World.simpleLoop2(zmap, pos)
#C.convertArrays2File(allObjects, 'out.plt')
