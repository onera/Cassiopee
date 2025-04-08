import CPlot
import Modeler.WorldZ as World
import Modeler.Models as Models
import Converter as C
import time

# Objects
staticNeutral = []; staticZMap = []; dynamicNeutral = []; dynamicZMap = []
allObjects = []

# Create Ground
zmap = World.createZMap(200, 200)
zmap = World.initGround(zmap, type=0)
ground = C.copy(zmap)
allObjects += [ground]

# Create static zmap objects
# Pyramide
box = Models.box((0,0,0), (40,40,1), chamfer=0.08)
box = World.placeObject(zmap,box,120,100,impact=1)
staticZMap.append(box)
box = Models.box((0,0,0), (32,32,1), chamfer=0.07)
box = World.placeObject(zmap,box,124,104,impact=1)
staticZMap.append(box)
box = Models.box((0,0,0), (24,24,1), chamfer=0.06)
box = World.placeObject(zmap,box,128,108,impact=1)
staticZMap.append(box)
box = Models.box((0,0,0), (16,16,1), chamfer=0.05)
box = World.placeObject(zmap,box,132,112,impact=1)
staticZMap.append(box)
box = Models.box((0,0,0), (8,8,1), chamfer=0.04)
box = World.placeObject(zmap,box,136,116,impact=1)
staticZMap.append(box)
box = Models.box((0,0,0), (3,3,1), chamfer=0.03)
box = World.placeObject(zmap,box,138.5,118.5,impact=1)
staticZMap.append(box)
box = Models.box((0,0,0), (1,1,3), chamfer=0.02)
box = World.placeObject(zmap,box,139.5,119.5,impact=1)
staticZMap.append(box)
# bloc isole
box = Models.box((0,0,0), (1,1,3), chamfer=0.02)
box = World.placeObject(zmap,box,80.,80.,impact=1)
staticZMap.append(box)
# Jumping blocks
box = Models.box((0,0,0), (3,0.5,2), chamfer=0.02)
box = World.placeObject(zmap,box,110.,50.,impact=1)
staticZMap.append(box)
box = Models.box((0,0,0), (3,0.5,2), chamfer=0.02)
box = World.placeObject(zmap,box,110.,40.,impact=1)
staticZMap.append(box)

#column = Models.column(R=1., N=1, h=1)
#for i in range(67,56):
#    staticZMap.append(World.placeRandomObject(zmap,column,impact=1,zshift=-0.5))
allObjects += staticZMap

# position du viewer (posCam, posEye, Height of viewer, deltaZ, speed)
pos = [(zmap[2]/2,zmap[3]/2,0), (zmap[2]/2+1,zmap[2]/2,0), 1., 1.5, 0.]
pos = World.placeViewer(zmap, pos)
CPlot.display(allObjects, displayInfo=0, bgColor=1, shadow=0, posCam=pos[0],
              posEye=pos[1], meshStyle=3) #, stereo=2, stereoDist=0.3)

#CPlot.setState(message='u: forward, j: backward, h: turn left, k: turn right')
CPlot.setState(activateShortCuts=0)
kstate = World.createKeyState()

t = 0
while 1 != 2:
    World.getKeys(kstate)
    if kstate['forward'] == 1: # forward
        pos = World.moveForward(zmap, pos, 0.3)
    if kstate['backward'] == 1: #backward
        pos = World.moveBackward(zmap, pos, 0.3)
    if kstate['left'] == 1: # left
        pos = World.turnLeft(zmap, pos, 5.)
    if kstate['right'] == 1: # right
        pos = World.turnRight(zmap, pos, 5.)
    if kstate['jump'] == 1: # jump
        pos = World.jump(zmap, pos, 2.)
        kstate['jump']=0
    pos = World.completeMotion(zmap, pos)
    time.sleep(0.05)
    t += 1
