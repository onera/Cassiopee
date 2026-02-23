import CPlot
import Modeler.WorldZ as World
import Modeler.Models as Models
import Transform as T
import Converter as C
import Generator as G
import Geom as D
import time

# Objects
staticNeutral=[]; staticZMap=[]; dynamicNeutral=[]; dynamicZMap=[]
allObjects=[]

# Create Ground
zmap = World.createZMap(50, 50)
zmap = World.initGround(zmap, type=0)
ground = C.copy(zmap)
allObjects += [ground]

# mur du bas
l = D.line((0,0,0), (49,0,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,0,0,impact=1)
staticZMap.append(box)
# mur du haut
l = D.line((0,0,0), (49,0,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,0,49,impact=1)
staticZMap.append(box)
# mur de droite
l = D.line((0,0,0), (0,49,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,49,0,impact=1)
staticZMap.append(box)
# demi mur de gauche
l = D.line((0,0,0), (0,20,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,0,0,impact=1)
staticZMap.append(box)
# demi mur de gauche
l = D.line((0,0,0), (0,20,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,0,40,impact=1)
staticZMap.append(box)

# couloir d'entree
l = D.line((0,0,0), (10,0,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,0,20,impact=1)
staticZMap.append(box)
l = D.line((0,0,0), (10,0,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,0,28,impact=1)
staticZMap.append(box)

# mur du fond (secret)
l = D.line((0,0,0), (24,0,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,0,4,impact=1)
staticZMap.append(box)
l = D.line((0,0,0), (19,0,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,28,4,impact=1)
staticZMap.append(box)
l = D.line((0,0,0), (27,0,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,5,8,impact=1)
staticZMap.append(box)

# mur vertical d'access
l = D.line((0,0,0), (0,8,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,5,8,impact=1)
staticZMap.append(box)
l = D.line((0,0,0), (5,0,0))
box = Models.wall(l, 2,1,1,5,0.08)
box = World.placeObject(zmap,box,5,16,impact=1)
staticZMap.append(box)
box = Models.panel(text='Ferme!',h=1.)
box = T.rotate(box, (0,0,0), (0,0,1), 90.)
box = World.placeObject(zmap,box,33.5,5,impact=1)
staticZMap.append(box)

# Paneau du milieu
box = Models.panel(text='Trouve 2\netoiles',h=0.2)
box = T.rotate(box, (0,0,0), (0,0,1), -90.)
box = World.placeObject(zmap,box,30.,20,impact=1)
staticZMap.append(box)

# Create stars
nstars = 2
targetCoords = []
star = Models.star(R1=0.5,R2=1.,N=10,h=0.5)
star = T.rotate(star, (0,0,0), (0,1,0), 90.)
st = World.placeObject(zmap,star,47,2,impact=0)
bb = G.bbox(st)
targetCoords.append((bb[0]+(bb[3]-bb[0])*0.5,bb[1]+(bb[4]-bb[1])*0.5))
staticNeutral.append(st)

st = World.placeObject(zmap,star,31,22,impact=0)
bb = G.bbox(st)
targetCoords.append((bb[0]+(bb[3]-bb[0])*0.5,bb[1]+(bb[4]-bb[1])*0.5))
staticNeutral.append(st)

allObjects += staticZMap+staticNeutral

# position du viewer (posCam, posEye, Height of viewer, deltaZ, speed)
#pos = [(zmap[2]/2,zmap[3]/2,0), (zmap[2]/2+1,zmap[2]/2,0), 1., 1., 0]
pos = [(0,25,0), (1,25,0),1,0.3,0]
pos = World.placeViewer(zmap, pos)
CPlot.display(allObjects, displayInfo=0, bgColor=1, shadow=0, posCam=pos[0],
              posEye=pos[1], meshStyle=3) #, stereo=2, stereoDist=0.3)

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

    for i in range(nstars):
        tg = targetCoords[i]
        if abs(pos[0][0]-tg[0])<3 and abs(pos[0][1]-tg[1])<3:
            o = World.placeObject(zmap,star,tg[0],tg[1],zshift=10.,impact=0)
            targetCoords[i] = (-10000,-10000)
            no = len(allObjects)-len(staticNeutral)+i
            CPlot.replace(allObjects, no, o)
            nstars -= 1
            CPlot.setState(message='%d stars to be found'%nstars)
