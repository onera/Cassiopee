# Animator2D 
# Billboard animation
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import CPlot.PyTree as CPlot
from operator import itemgetter

class Animator2D:
    def __init__(self):
        self.t = C.newPyTree(['Images'])
        self.images = {}
        self.nimages = 0
        self.line = []

    # tri t suivant zpos
    def sort(self):
        zones = Internal.getZones(self.t)
        out = []
        for z in zones:
            zpos = Internal.getNodeFromName1(z, 'zpos')
            zpos = Internal.getValue(zpos)
            out.append( (zpos, z) )
        so = sorted(out, key=itemgetter(0))
        zs = []
        for s in so:
            zs.append(so[1])
        self.t[2][1][2] = zs

    # Register an image
    def registerImage(self, key, fileName, Ni=1, Nj=1, speed=0.1):
        # nom du fichier image, Ni, Nj, no de l'image courante, 
        # vitesse de l'animation, ptr courant de l'animation
        self.images[key] = [fileName,Ni,Nj,self.nimages,speed,0.]
        self.line += [fileName,Ni,Nj]
        self.nimages += 1

    # Return [0.,2] from key
    def getShaderNo(self, key):
        n = self.images[key][3]
        sp = (n+1)/(self.nimages*0.5)
        sp = max(sp, 0.1)
        sp = min(sp, 2.)
        return sp

    # set CPlot from registered images
    def prepare(self):
        # load all billboards images
        CPlot.setState(billBoards=self.line, billBoardSize=0.8)

    # Update all shaderparam (no need to call it anymore)
    def updateShaderParam(self):
        zones = Internal.getZones(self.t)
        for z in zones:
            param = Internal.getNodeFromName(z, 'ShaderParameters')
            print (z[0], self.getShaderNo(z[0]))
            param[1][1] = self.getShaderNo(z[0])
            
    # display all t
    # you must register all images before display
    def openDisplay(self):
        self.prepare()
        CPlot.display(self.t, mode='Render', bgColor=1, posCam=(0.,0,1), 
                      posEye=(0.,0,0), dirCam=(0,1,0))
        CPlot.setState(activateShortCuts=0, displayInfo=0)

    # move key right
    def moveRight(self, key, xr=0.01):
        self.move(key, (xr,0,0))

    def moveLeft(self, key, xr=0.01):
        self.move(key, (-xr,0,0))

    def moveUp(self, key, xr=0.01):
        self.move(key, (0,xr,0))
    
    def moveDown(self, key, xr=0.01):
        self.move(key, (0,-xr,0))

    # move a key
    def move(self, key, v):
        if len(v) == 2: v = (v[0],v[1],0.)
        zones = Internal.getZones(self.t)
        noz = 0
        for z in zones:
            if z[0] == key: break
            noz += 1
        z = Internal.getNodeFromName2(self.t, key)
        zp = T.translate(z, v)
        CPlot.replace(self.t,1,noz,zp)
        CPlot.render()

    # Return zone of key
    def getZone(self, key):
        z = Internal.getNodeFromName2(self.t, key)
        return z

    # Return the coord of the first point of key
    def getPos(self, key):
        z = self.getZone(key)
        px = Internal.getNodeFromName2(z, 'CoordinateX')
        py = Internal.getNodeFromName2(z, 'CoordinateY')
        pz = Internal.getNodeFromName2(z, 'CoordinateZ')
        return (px[1].ravel()[0],py[1].ravel()[0],pz[1].ravel()[0])

    # draw an image
    def drawImage(self, key, imageKey=None, pos=(0,0,0), scale=1.):
        if imageKey is None: imageKey = key
        a = D.point(pos)
        CPlot._addRender2Zone(a, material='Sphere', shaderParameters=[30.*scale,self.getShaderNo(imageKey)])
        a[0] = key
        C._initVars(a, 'TBB__', 0.)
        Internal._createChild(a, 'zpos', 'DataArray_t', value=pos[2])
        CPlot.add(self.t, 1, -1, a)
        CPlot.render()

    # draw a text
    def drawText(self, key, pos=(0,0,0), scale=1., text="youpi", color='Black',
                 h=0.):
        if h == 0.: a = D.text2D(text, font='vera')
        else: 
            a = D.text3D(text, font='vera')
            T._contract(a, (0,0,0), (1,0,0), (0,1,0), h)
        T._homothety(a, (0,0,0), scale)
        b = G.bbox(a)
        v = (pos[0]-b[0],pos[1]-b[1],pos[2]-b[2])
        T._translate(a, v)
        a[0] = key
        C._initVars(a, 'TBB__', 0.)
        CPlot._addRender2Zone(a, color=color)
        Internal._createChild(a, 'zpos', 'DataArray_t', value=pos[2])
        CPlot.add(self.t, 1, -1, a)
        CPlot.render()

    # erase a key
    def erase(self, key, render=True):
        # Replace par un point hors cadre
        a = D.point((0,0,0))
        noz = 0
        zones = Internal.getZones(self.t)
        for z in zones:
            if z[0] == key: break
            noz += 1
        CPlot.replace(self.t, 1, noz, a)
        if render: CPlot.render()

    # Animate one step of an image
    # si curMin et curMax fournis, loop entre curMin et curMax
    def animate(self, key, curMin=None, curMax=None):
        a = self.getZone(key)
        im = self.images[key]
        speed = im[4]
        cur = im[5]
        cur += speed
        if curMin is None:
            if cur > 1.: cur = 0. # loop
        else:
            if cur >= curMax: cur = curMin # loop
        im[5] = cur
        C._initVars(a, 'TBB__', cur)
        noz = 0
        zones = Internal.getZones(self.t)
        for z in zones:
            if z[0] == key: break
            noz += 1
        CPlot.replace(self.t, 1, noz, a)
        CPlot.render()
