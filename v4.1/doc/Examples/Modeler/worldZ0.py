# - worldZ test 0 -
import CPlot
import Modeler.WorldZ as World
import Modeler.Models as Models
import Converter as C

# Objects
allObjects = []

# Create Ground
zmap = World.createZMap(200, 200)
zmap = World.initGround(zmap, type=0)
ground = C.copy(zmap)
allObjects += [ground]

# Create static zmap objects
box = Models.box((0,0,0), (40,40,1), chamfer=0.08)
box = World.placeObject(zmap,box,120,100,impact=1)
allObjects.append(box)
box = Models.box((0,0,0), (32,32,1), chamfer=0.07)
box = World.placeObject(zmap,box,124,104,impact=1)
allObjects.append(box)
box = Models.box((0,0,0), (24,24,1), chamfer=0.06)
box = World.placeObject(zmap,box,128,108,impact=1)
allObjects.append(box)
box = Models.box((0,0,0), (16,16,1), chamfer=0.05)
box = World.placeObject(zmap,box,132,112,impact=1)
allObjects.append(box)
box = Models.box((0,0,0), (8,8,1), chamfer=0.04)
box = World.placeObject(zmap,box,136,116,impact=1)
allObjects.append(box)
box = Models.box((0,0,0), (3,3,1), chamfer=0.03)
box = World.placeObject(zmap,box,138.5,118.5,impact=1)
allObjects.append(box)
box = Models.box((0,0,0), (1,1,1), chamfer=0.02)
box = World.placeObject(zmap,box,139.5,119.5,impact=1)
allObjects.append(box)

# position du viewer (posCam, posEye, Height of viewer, deltaZ)
pos = World.placeViewerAtCenter(zmap)
CPlot.display(allObjects, displayInfo=0, bgColor=1, shadow=0, posCam=pos[0],
              posEye=pos[1], meshStyle=3) #, stereo=2, stereoDist=0.3)

World.simpleLoop(zmap, pos)
