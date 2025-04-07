import CPlot
import Modeler.WorldZ as World
import Modeler.Models as Models
import Geom as D
import Transform as T
import Converter as C
import Generator as G
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
#box = Models.box((0,0,0), (2,2,10))
box = Models.panel(text="Trouve les\netoiles!", h=1.5)
box = World.placeObject(zmap,box,zmap[2]/2+10,zmap[3]/2+10,impact=1)
staticZMap.append(box)
wall = Models.wall(D.line((0,0,0), (200,0,0)), 2., 1., 1, nlayers=20, chamfer=0.1)
wall = World.placeObject(zmap,wall,0,198,impact=1)
staticZMap.append(wall)

column = Models.column2(R1=2., R2=1.7,N=10, h=15)
for i in range(10):
    staticZMap.append(World.placeRandomObject(zmap,column,impact=1,zshift=-0.5))
allObjects += staticZMap

# Create static neutral objects
star = Models.star(R1=0.5,R2=1.,N=10,h=0.5)
star = T.rotate(star, (0,0,0), (0,1,0), 90.)
nstars = 5
targetCoords = []
for i in range(nstars):
    st = World.placeRandomObject(zmap,star,impact=0)
    bb = G.bbox(st)
    targetCoords.append((bb[0]+(bb[3]-bb[0])*0.5,bb[1]+(bb[4]-bb[1])*0.5))
    staticNeutral.append(st)
allObjects += staticNeutral

# position du viewer (posCam, posEye, Height of viewer, deltaZ, speed)
pos = [(zmap[2]/2,zmap[3]/2,0), (zmap[2]/2+1,zmap[2]/2,0), 1., 1., 0.]
pos = World.placeViewer(zmap, pos)
CPlot.display(allObjects, displayInfo=0, bgColor=1, shadow=0, posCam=pos[0],
              posEye=pos[1], meshStyle=3) #, stereo=2, stereoDist=0.3)

#CPlot.setState(message='u: forward, j: backward, h: turn left, k: turn right')
CPlot.setState(message='%d stars to be found'%nstars)
CPlot.setState(activateShortCuts=0)
kstate = World.createKeyState()

t = 0
while 1 != 2:
    World.getKeys(kstate)
    if kstate['forward'] == 1: # forward
        pos = World.moveForward(zmap, pos, 0.5)
    if kstate['backward'] == 1: #backward
        pos = World.moveBackward(zmap, pos, 0.5)
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

    # petit jeu de chasse au tresor
    for i in range(5):
        tg = targetCoords[i]
        if abs(pos[0][0]-tg[0])<3 and abs(pos[0][1]-tg[1])<3:
            o = World.placeObject(zmap,star,tg[0],tg[1],zshift=10.,impact=0)
            targetCoords[i] = (-10000,-10000)
            no = len(allObjects)-len(staticNeutral)+i
            CPlot.replace(allObjects, no, o)
            nstars -= 1
            CPlot.setState(message='%d stars to be found'%nstars)
