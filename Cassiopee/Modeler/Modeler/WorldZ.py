#
# Cassiopee's world
#
import math
import Converter
import Generator
import Transform
import CPlot
import KCore.Vector as Vector
import random
import time

# - create a zmap -
def createZMap(ni,nj):
    a = Generator.cart((0,0,0), (1,1,1), (ni,nj,1))
    return a

# - initGround -
def initGround(a, type=0):
    if type == 0:
        Converter._initVars(a, '{z}=0.')
    elif type == 1:
        Converter._initVars(a, '{z}=2*cos(0.1*{x})*sin(0.1*{y})')
    elif type == 2:
        Converter._initVars(a, '{z}=3.*(3*cos(0.01*{x})+cos(0.1*{x}))*sin(0.01*{y})')
    else:
        Converter._initVars(a, '{z}=floor(5*cos(0.1*{x})*sin(0.01*{y}))')
    return a

# -getZ-
# Retourne le z de (x,y) sur la zmap
# IN: xy
# OUT: z of zmap (interpolated)
def getZ(zmap, x, y):
    ni = zmap[2]; nj = zmap[3]
    zm = zmap[1][2]
    vi = int(x); vj = int(y)
    vi = max(vi,0); vi = min(vi,ni-2)
    vj = max(vj,0); vj = min(vj,nj-2)
    vi1 = vi+1; vj1 = vj+1
    alpha = x-vi; beta = y-vj
    alpha1 = 1-alpha; beta1 = 1-beta
    z1 = zm[vi+ni*vj]
    z2 = zm[vi1+ni*vj]
    z3 = zm[vi+ni*vj1]
    z4 = zm[vi1+ni*vj1]
    z = alpha1*beta1*z1+alpha*beta1*z2+alpha1*beta*z3+alpha*beta*z4
    return z

# - place viewer on zmap -
# Compute z par interpolation + posEye
def placeViewer(zmap, pos):
    (x,y,z) = pos[0]
    (xe,ye,ze) = pos[1]
    H = pos[2]
    deltaZ = pos[3]
    sn = pos[4]
    z = getZ(zmap,x,y)
    ze = getZ(zmap,xe,ye)
    posCam = (x,y,z+H)
    posEye = (xe,ye,ze+H)
    CPlot.setState(posCam=posCam, posEye=posEye, dirCam=(0,0,1))
    pos = [posCam, posEye, H, deltaZ, sn]
    return pos

def placeViewerAtCenter(zmap, h=1., deltaZ=1., sn=0.):
    pos = [(zmap[2]/2,zmap[3]/2,0), (zmap[2]/2+1,zmap[2]/2,0), h, deltaZ, sn]
    pos = placeViewer(zmap, pos)
    return pos

# - create keyState -
def createKeyState():
    return {'forward':0, 'backward':0, 'left':0, 'right':0, 'jump':0}

# - getkeys -
# Modifie le key state en fonction des touches pressees ou relachees
def getKeys(state):
    key = CPlot.getKeyboard()
    if len(key) == 0: return
    for k in range(len(key)):
        try: v = ord(key[k])
        except: v = key
        if v == 1: state['forward']=1
        elif v == 5: state['forward']=0
        elif v == 2: state['backward']=1
        elif v == 6: state['backward']=0
        elif v == 3: state['left']=1
        elif v == 7: state['left']=0
        elif v == 4: state['right']=1
        elif v == 8: state['right']=0
        elif v == 106: state['jump']=1
    CPlot.resetKeyboard()

# - move viewer forward -
# IN: pos: position
# IN: unit: de combien on se deplace
# IN: headstop: limite les mouvements de tete
def moveForward(zmap, pos, unit=1., headStop=0.1):
    (x,y,z) = pos[0] # position of viewer
    posCam = pos[0]
    (xe,ye,ze) = pos[1] # posEye
    H = pos[2] # hauteur viewer
    deltaZ = pos[3] # deltaZ permis
    sn = pos[4] # speed
    sn = max(sn, unit)
    sn = min(0.001+sn, 2*unit)
    #CPlot.setState(message='speed=%f'%sn)
    unit = sn

    delta = (xe-x,ye-y,0.)
    delta = Vector.normalize(delta)
    delta = Vector.mul(unit, delta)
    delta2 = Vector.mul(2., delta)
    (xp,yp,zp) = Vector.add(posCam, delta)
    (xpe,ype,zpe) = Vector.add(posCam, delta2)
    zp = getZ(zmap,xp,yp)
    zpe = getZ(zmap,xpe,ype)

    if zpe-zp>headStop*unit: zpe = zp+headStop*unit
    elif zpe-zp<-headStop*unit: zpe = zp-headStop*unit
    if abs(zp+H-z) < deltaZ*unit:
        posCam = (xp,yp,zp+H)
        posEye = (xpe,ype,zpe+H)
        CPlot.setState(posCam=posCam, posEye=posEye, dirCam=(0,0,1))
        pos = [posCam, posEye, H, deltaZ, sn]
    return pos

def moveBackward(zmap, pos, unit=1.):
    (x,y,z) = pos[0]
    posCam = pos[0]
    (xe,ye,ze) = pos[1]
    H = pos[2]
    deltaZ = pos[3]
    sn = 0.
    delta = (xe-x,ye-y,0.)
    delta = Vector.normalize(delta)
    delta = Vector.mul(unit, delta)
    (xp,yp,zp) = Vector.sub(posCam, delta)
    (xpe,ype,zpe) = posCam
    zp = getZ(zmap,xp,yp)
    if abs(zp+H-z) < deltaZ*unit:
        posCam = (xp,yp,zp+H)
        posEye = (xpe,ype,zpe)
        CPlot.setState(posCam=posCam, posEye=posEye, dirCam=(0,0,1))
        pos = [posCam, posEye, H, deltaZ, sn]
    return pos

def turnLeft(zmap, pos, angle=5., unit=1., headStop=0.1):
    (x,y,z) = pos[0]
    (xe,ye,ze) = pos[1]
    H = pos[2]
    deltaZ = pos[3]
    sn = 0.
    vcos = math.cos(angle*math.pi/180.)
    vsin = math.sin(angle*math.pi/180.)
    (deltax,deltay,deltaz) = (xe-x,ye-y,ze-z)
    px = vcos*deltax-vsin*deltay
    py = vsin*deltax+vcos*deltay
    posCam = (x,y,z)
    zpe = getZ(zmap,x+px,y+py)
    if zpe+H-z>headStop*unit: zpe = z-H+headStop*unit; limit = True
    elif zpe+H-z<-headStop*unit: zpe = z-H-headStop*unit; limit = True
    posEye = (x+px,y+py,zpe+H)
    CPlot.setState(posCam=posCam, posEye=posEye, dirCam=(0,0,1))
    pos = [posCam, posEye, H, deltaZ, sn]
    return pos

def turnRight(zmap, pos, angle=5., unit=1., headStop=0.1):
    (x,y,z) = pos[0]
    (xe,ye,ze) = pos[1]
    H = pos[2]
    deltaZ = pos[3]
    sn = 0.
    vcos = math.cos(angle*math.pi/180.)
    vsin = math.sin(angle*math.pi/180.)
    (deltax,deltay,deltaz) = (xe-x,ye-y,ze-z)
    px = vcos*deltax+vsin*deltay
    py = -vsin*deltax+vcos*deltay
    posCam = (x,y,z)
    zpe = getZ(zmap,x+px,y+py)
    if zpe+H-z>headStop*unit: zpe = z-H+headStop*unit; limit = True
    elif zpe+H-z<-headStop*unit: zpe = z-H-headStop*unit; limit = True
    posEye = (x+px,y+py,zpe+H)
    CPlot.setState(posCam=posCam, posEye=posEye, dirCam=(0,0,1))
    pos = [posCam, posEye, H, deltaZ, sn]
    return pos

def jump(zmap, pos, unit=10.):
    (x,y,z) = pos[0]
    (xe,ye,ze) = pos[1]
    z0 = z
    H = pos[2]
    sn = pos[4]
    delta = (xe-x,ye-y,0.)
    delta = Vector.normalize(delta)
    delta = Vector.mul(sn, delta)
    z += unit
    x += delta[0]; y += delta[1]
    posCam = (x,y,z)
    posEye = (xe+delta[0],ye+delta[1],ze+z-z0)
    pos = [posCam, posEye, H, pos[3], sn]
    CPlot.setState(posCam=posCam, posEye=posEye, dirCam=(0,0,1))
    return pos

# Finis le mouvement du viewer pour retrouver le sol
def completeMotion(zmap, pos, unit=0.3):
    (x,y,z) = pos[0]
    z0 = z
    H = pos[2]
    sn = pos[4]
    zp = getZ(zmap, x, y)
    if z == zp+H: return pos # deja au sol
    (xe,ye,ze) = pos[1]
    z -= unit
    delta = (xe-x,ye-y,0.)
    delta = Vector.normalize(delta)
    delta = Vector.mul(0.8*sn, delta)
    x += delta[0]; y += delta[1]
    if z < zp+H: z = zp+H
    posCam = (x,y,z)
    posEye = (xe+delta[0],ye+delta[1],ze+z-z0)
    CPlot.setState(posCam=posCam, posEye=posEye, dirCam=(0,0,1))
    pos = [posCam, posEye, H, pos[3], sn]
    return pos

# - placeObject -
# Place object on zmap at given position
def placeObject(zmap, o, x, y, zshift=0., impact=0):
    bb = Generator.bbox(o)
    ni1 = zmap[2]-1; nj1 = zmap[3]-1
    if x > ni1-bb[3]+bb[0]: x = ni1-bb[3]+bb[0]
    if y > nj1-bb[4]+bb[1]: y = nj1-bb[4]+bb[1]
    dx = x-bb[0]; dy = y-bb[1]
    z = getZ(zmap,x+0.5*(bb[3]-bb[0]),y+0.5*(bb[4]-bb[1]))
    dz = z-bb[2]+zshift
    o2 = Transform.translate(o, (dx,dy,dz))
    if impact == 1: # impact zmap with bbox
        bb2 = [bb[0]+dx,bb[1]+dy,bb[2]+dz,bb[3]+dx,bb[4]+dy,bb[5]+dz]
        vi = int(bb2[0]); vj = int(bb2[1])
        vi2 = int(bb2[3])+1; vj2 = int(bb2[4])+1
        ni = zmap[2]; nj = zmap[3]
        vi = max(vi,0); vi2 = min(vi2,ni)
        vj = max(vj,0); vj2 = min(vj2,nj)
        zm = zmap[1][2]
        zm = zm.reshape((nj,ni))
        zm[vj:vj2,vi:vi2] = bb2[5]
    return o2

# - place at random on zmap -
# IN: zmap: zmap object
# IN: o: object to place
# IN: zshift (if any)
# IN: impact: if 1 modify zmap
# IN: rotateZ: if true, rotate o at random in z
def placeRandomObject(zmap, o, zshift=0., impact=0, rotateZ=False):
    ni1 = zmap[2]-1; nj1 = zmap[3]-1
    posX = random.random()*ni1
    posY = random.random()*nj1
    if rotateZ == True:
        alpha = random.random()*90.
        o = Transform.rotate(o, (0,0,0), (0,0,1), alpha)
    return placeObject(zmap, o, posX, posY, zshift, impact)

# - simple loop -
# Only enables user to move on zmap
# No jump, just move
def simpleLoop(zmap, pos, unit=0.3):
    CPlot.setState(activateShortCuts=0)
    t = 0
    while 1 != 2:
        key = CPlot.getKeyboard()
        key = key[-1:]
        try: v = ord(key)
        except: v = key
        if v == 1: # forward
            pos = moveForward(zmap, pos, unit)
            CPlot.resetKeyboard()
        elif v == 2: #backward
            pos = moveBackward(zmap, pos, unit)
            CPlot.resetKeyboard()
        elif v == 3: # left
            pos = turnLeft(zmap, pos, 5.)
            CPlot.resetKeyboard()
        elif v == 4: # right
            pos = turnRight(zmap, pos, 5.)
            CPlot.resetKeyboard()
        t += 1

# Meilleur gestion du clavier
def simpleLoop2(zmap, pos, unit=0.3):
    CPlot.setState(activateShortCuts=0)
    kstate = createKeyState()
    t = 0
    while 1 != 2:
        getKeys(kstate)
        if kstate['forward'] == 1: # forward
            pos = moveForward(zmap, pos, unit)
        if kstate['backward'] == 1: #backward
            pos = moveBackward(zmap, pos, unit)
        if kstate['left'] == 1: # left
            pos = turnLeft(zmap, pos, 5.)
        if kstate['right'] == 1: # right
            pos = turnRight(zmap, pos, 5.)
        if kstate['jump'] == 1: # jump
            pos = jump(zmap, pos, 1.+4*unit)
            kstate['jump']=0
        pos = completeMotion(zmap, pos)
        time.sleep(0.05)
        t += 1
