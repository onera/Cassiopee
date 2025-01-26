# a 2D world
# 0: O: nothing
# 1: X: wall
# 2: P: pastille orange (finish level)
# 3: S: pastille bleu (slow)
# 4: F: pastille rouge (fast)

import Generator.PyTree as G
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import numpy

# Global base speed (depend des perfos)
SPEED = 0.35

class World:
    def __init__(self, fileName=None):
        # World title
        self.title = '2D World'
        # World author
        self.author = 'Ash'
        # Music file
        self.music = 'world.wav'
        # Hard ware speed (units)
        self.hardSpeed = SPEED
        # Base world speed (en %)
        self.baseSpeed = 100.
        # World time
        self.time = 0
        # World size (nbre de cases)
        self.nx = 0
        self.ny = 0
        # Word init pos
        self.ix = -1
        self.iy = -1
        # World map
        self.map = numpy.zeros( (0,0), dtype=numpy.int32 )
        # pyTree world
        self.world = None
        # pyTree all
        self.all = C.newPyTree(['World', 'Objects', 'Persos'])
        # object i,j-> no de la zone correspondante
        self.object = None
        # build map
        if fileName is None: self.buildMapTest()
        else: self.buildMapFromFile(fileName)
        # build models
        self.build()

    # build sample empty world for test
    def buildMapTest(self):
        self.title = 'Test'
        self.author = 'Ash'
        self.music = None
        self.baseSpeed = 100.
        self.ix = 5; self.iy = 5
        self.nx = 10; self.ny = 10
        self.map = numpy.zeros( (self.nx, self.ny), dtype=numpy.int32 )
        self.map[5,5] = 0

    # build objs dict from file
    def buildMapFromFile(self, file):
        f = open(file, 'r')
        data = f.read()
        f.close()
        lines = data.splitlines()
        self.title = lines[0]
        self.author = lines[1]
        self.music = lines[2]
        r = lines[3].split(';')
        self.baseSpeed = float(r[0])
        if r[1] == 'm': self.ix = -1
        else: self.ix = int(r[1])
        if r[2] == 'm': self.iy = -1
        else: self.iy = int(r[2])
        out = []
        for l in lines[4:]:
            if l != '' and l != ' ': out.append(l)
        lines = out
        self.ny = len(lines)
        self.nx = 0
        for l in lines:
            self.nx = max(len(l), self.nx)
        self.map = numpy.zeros( (self.nx, self.ny), dtype=numpy.int32 )
        self.object = numpy.zeros( (self.nx, self.ny), dtype=numpy.int32 )
        j = 0; pc = 0
        for l in lines:
            nl = len(l)
            for i in range(nl):
                if l[i] == 'X' or l[i] == 'x': # wall
                    self.map[i,j] = 1
                elif l[i] == 'P' or l[i] == 'p': # Pastille orange
                    self.map[i,j] = 2
                    self.object[i,j] = pc; pc += 1
                elif l[i] == 'S' or l[i] == 's': # pastille slow
                    self.map[i,j] = 3
                    self.object[i,j] = pc; pc += 1
                elif l[i] == 'F' or l[i] == 'f': # pastille fast
                    self.map[i,j] = 4
                    self.object[i,j] = pc; pc += 1
                elif l[i] == 'H' or l[i] == 'h': # hollow wall
                    self.map[i,j] = 5
                elif l[i] == 'Z' or l[i] == 'z': # rocket spring
                    self.map[i,j] = 5
            j += 1


    # Create objects from map
    def build(self):
        # Create square
        a = G.cart((0,0,0), (1,1,1), (self.nx+1,self.ny+1,1))
        CPlot._addRender2Zone(a, color='White', meshOverlay=True)
        self.world = [a]
        for j in range(self.ny):
            for i in range(self.nx):
                if self.map[i,j] == 1: # wall
                    a = G.cart((i,self.ny-j-1,0), (1,1,2), (2,2,2))
                    CPlot._addRender2Zone(a, color='White', meshOverlay=True)
                    self.world += [a]
                elif self.map[i,j] == 2: # orange pill
                    a = G.cylinder((i+0.5,self.ny-j-0.5,0.), 0., 0.5, 0., 360., 0.3, (20,2,2))
                    CPlot._addRender2Zone(a, color='Orange')
                    self.all[2][2][2].append(a)
                elif self.map[i,j] == 3: # blue pill
                    a = G.cylinder((i+0.5,self.ny-j-0.5,0.), 0., 0.5, 0., 360., 0.3, (20,2,2))
                    CPlot._addRender2Zone(a, color='Blue')
                    self.all[2][2][2].append(a)
                elif self.map[i,j] == 4: # red pill
                    a = G.cylinder((i+0.5,self.ny-j-0.5,0.), 0., 0.5, 0., 360., 0.3, (20,2,2))
                    CPlot._addRender2Zone(a, color='Red')
                    self.all[2][2][2].append(a)
                elif self.map[i,j] == 5: # hollow wall
                    a = G.cart((i,self.ny-j-1,0), (1,1,2), (2,2,2))
                    CPlot._addRender2Zone(a, color='Grey', meshOverlay=True)
                    self.world += [a]
                    self.map[i,j] = 0
                elif self.map[i,j] == 6: # rocket spring
                    a = G.cart((i,self.ny-j-1,0), (1,1,0.5), (2,2,2))
                    CPlot._addRender2Zone(a, color='Blue')
                    self.world += [a]
        self.all[2][1][2] = self.world
        C._initVars(self.all, 'TBB__', 0.)

    # Static first display of static world
    def render(self):
        CPlot.display(self.all, bgColor=1, shadow=0, mode='render',
                      posCam=(self.nx*0.5, self.ny*0.5, 20.),
                      posEye=(self.nx*0.5,self.ny*0.5,0),
                      dirCam=(0,1,0), meshStyle=3)
        CPlot.setState(activateShortCuts=0, displayInfo=0)

    def advanceTime(self):
        self.time += 0.1
