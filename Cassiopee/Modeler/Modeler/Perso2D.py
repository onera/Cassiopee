# Perso class

import Geom.PyTree as D
import Transform.PyTree as T
import CPlot.PyTree as CPlot
import Converter.Internal as Internal
import Converter.PyTree as C
from . import Sound
from . import Animator2D

animator = Animator2D.Animator2D()

class perso:
    def __init__(self, name='serge', file='Images/serge2.png'):
        # name
        self.name = name
        # position in world
        self.x = 0
        self.y = 0
        # local speed
        self.vx = 0.
        self.vy = 0.
        # speed : vitesse actuelle (avec buff)
        self.speed = 0
        # base speed : vitesse de base (sans buff)
        self.baseSpeed = 100
        # time of last buff
        self.timeBuff = 0
        # model
        self.piece = None
        # sounds
        self.gob = None
        self.faster = None
        # images for billboards (sprites)
        self.images = self.registerImage(file)
        # count. Un compteur pour les perso se deplacant automatiquement
        self.count = 0

    # register sounds
    def registerSounds(self):
        self.gob = Sound.registerSound("gob2.wav")  
        self.faster = Sound.registerSound("gettin' faster.wav")

    # register images
    def registerImage(self, file):
        animator.registerImage(self.name, file, 4, 4, speed=1./16.)

    # play a sound
    def playSound(self, s):
        if s is not None:
            Sound.closeAllSounds()
            Sound.playSound(s)

    # Positionne perso init
    def place(self, world, i=None, j=None):
        m = world.map
        nx = world.nx
        ny = world.ny
        if i is None:
            if world.ix == -1: i = int(nx*0.5)
            else: i = world.ix
        if j is None:
            if world.iy == -1: j = int(ny*0.5)
            else: j = world.ny-world.iy-1
        if m[i,j] == 0:
            self.x = i+0.5
            self.y = world.ny-j-0.5
            self.baseSpeed = world.hardSpeed*world.baseSpeed/100.
            self.speed = self.baseSpeed
            self.vx = 0.
            self.vy = 0.
        else:
            print ('>>init pb')

    # Build the model (sphere)
    def build(self, world):
        a = D.sphere((0,0,0),0.5,N=8)
        CPlot._addRender2Zone(a, color='Blue')
        self.piece = a
        world.all[2][3][2].append(self.piece)

    def build1(self, world):
        a = D.point((0,0,0))
        CPlot._addRender2Zone(a, material='Sphere', shaderParameters=[30.*1.,animator.getShaderNo(self.name)])
        C._initVars(a, 'TBB__', 13./16.)
        animator.images[self.name][5] = 13./16.
        Internal._createChild(a, 'zpos', 'DataArray_t', value=0.001)
        a[0] = self.name
        self.piece = a
        world.all[2][3][2].append(self.piece)

    # Update x and y from vx,vy
    def move(self, world):
        self.x += self.vx
        self.y += self.vy
        if self.x > world.nx: self.x = 0.
        elif self.x < 0: self.x = world.nx
        if self.y > world.ny: self.y = 0.
        elif self.y < 0: self.y = world.ny

    # Checking collision with world
    def checkCollision(self, world):
        # Check buffs
        if self.timeBuff > 0 and world.time - self.timeBuff > 100:
            self.timeBuff = 0
            self.speed = self.baseSpeed

        # Check collision
        m = world.map
        x = self.x; y = self.y
        # decalage lie a la taille du personnage
        if self.vx > 0: x += 0.3
        if self.vx < 0: x -= 0.3
        if self.vy > 0: y += 0.3
        if self.vy < 0: y -= 0.3

        i = int(x)
        j = world.ny-int(y)-1
        i = max(i,0); i = min(i,world.nx-1)
        j = max(j,0); j = min(j,world.ny-1)
        #print i,j,m[i,j]
        if m[i,j] == 1: # wall 
            self.playSound(self.gob)
            self.x -= self.vx
            self.y -= self.vy
            return 1
        elif m[i,j] == 2: # orange pill
            self.playSound(self.gob)
            # supprime la pastille de world
            world.map[i,j] = 0
            pc = world.object[i,j]
            #i = CPlot.getCPlotNumber(world.all, 'Pastilles', world.all[2][2][2][pc][0])
            #CPlot.delete([i])
            a = D.point( (i+0.5,world.ny-j-0.5,0.) )
            CPlot.replace(world.all, 2, pc, a)
            world.object[i,j] = 0
            return 2
        elif m[i,j] == 3: # Slow blue pill
            self.playSound(self.gob)
            # supprime la pastille de world
            world.map[i,j] = 0
            pc = world.object[i,j]
            a = D.point( (i+0.5,world.ny-j-0.5,0.) )
            CPlot.replace(world.all, 2, pc, a)
            world.object[i,j] = 0
            self.timeBuff = world.time
            self.speed = 0.7*self.speed
            return 3
        elif m[i,j] == 4: # Fast red pill
            self.playSound(self.gob)
            # supprime la pastille de world
            world.map[i,j] = 0
            pc = world.object[i,j]
            a = D.point( (i+0.5,world.ny-j-0.5,0.) )
            CPlot.replace(world.all, 2, pc, a)
            world.object[i,j] = 0
            self.timeBuff = world.time
            self.speed = 1.3*self.speed
            Sound.playSound(self.faster)
            return 4
        return 0

    # Render a rotating sphere
    def render(self, world, init=False, centerCam=True):
        if centerCam:
            CPlot.setState(posCam=(self.x,self.y,20), posEye=(self.x,self.y,0))
        if self.vx > 0:
            T._rotate(self.piece, (0,0,0), (0,1,0), self.speed*100.)
        if self.vy > 0:
            T._rotate(self.piece, (0,0,0), (1,0,0), -self.speed*100.)
        if self.vx < 0:
            T._rotate(self.piece, (0,0,0), (0,1,0), -self.speed*100.)
        if self.vy < 0:
            T._rotate(self.piece, (0,0,0), (1,0,0), self.speed*100.)
        b = T.translate(self.piece, (self.x, self.y,0.))
        if init:
            CPlot.add(world.all, 3, 0, b)
        else:
            CPlot.replace(world.all, 3, 0, b)
        self.vx = 0.
        self.vy = 0.

    # Render a billboard
    def render1(self, world, init=False, centerCam=True):
        if centerCam:
            CPlot.setState(posCam=(self.x,self.y,10), posEye=(self.x,self.y,0))
        if self.vx > 0:
            self.animate(5./16.,8./16.)
        if self.vy > 0:
            self.animate(1./16.,4./16.)
        if self.vx < 0:
            self.animate(9./16.,12./16.)
        if self.vy < 0: 
            self.animate(13./16.,16./16.)
        b = T.translate(self.piece, (self.x, self.y, 0.001))
        zones = Internal.getZones(world.all[2][3])
        if self.name == 'hero': noz = 0
        elif self.name == 'serge': noz = 1
        if init:
            CPlot.add(world.all, 3, noz, b)
        else:
            CPlot.replace(world.all, 3, noz, b)
        self.vx = 0.; self.vy = 0.

    def animate(self, curMin, curMax):
        im = animator.images[self.name]
        speed = im[4]
        cur = im[5]
        cur += speed
        if cur < curMin: cur = curMin
        if cur > curMax: cur = curMin
        im[5] = cur
        C._initVars(self.piece, 'TBB__', cur)

    def getCur(self):
        im = animator.images[self.name]
        return im[5]

    def moveUp(self):
        self.vx = 0.
        self.vy = self.speed

    def moveDown(self):
        self.vx = 0.
        self.vy = -self.speed

    def moveLeft(self):
        self.vx = -self.speed
        self.vy = 0.

    def moveRight(self):
        self.vx = self.speed
        self.vy = 0.

    def autoMove(self):
        if self.count < 3: self.moveRight()
        elif self.count < 8: pass
        elif self.count < 10: self.moveDown()
        elif self.count < 13: self.moveLeft()
        elif self.count < 17: pass
        elif self.count < 19: self.moveUp()
        self.count += 1
        if self.count == 19: self.count = 0 
