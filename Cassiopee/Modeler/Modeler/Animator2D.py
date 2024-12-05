"""Animation module in 2D."""
# Billboard animation
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import CPlot.PyTree as CPlot
from operator import itemgetter
import math
import time

class Animator2D:
    def __init__(self):
        self.t = C.newPyTree(['Images'])
        self.images = {}
        self.nimages = 0
        self.line = []
        self.firstFit = False
        self.createKeyState()

    # tri t suivant zpos
    # normalement, cette fonction ne devrait pas etre appelee
    # car le tri doit etre fait au drawImage
    def sort(self):
        zones = Internal.getZones(self.t)
        out = []
        for z in zones:
            zpos = Internal.getNodeFromName1(z, 'zpos')
            if zpos is not None: zpos = Internal.getValue(zpos)
            else: zpos = 0.
            out.append( (zpos, z) )
        #so = sorted(out, key=itemgetter(0))
        out.sort(key=itemgetter(0))
        zs = []
        for s in out: zs.append(s[1])
        self.t[2][1][2] = zs
        # Must redisplay to keep ploter coherent
        CPlot.display(self.t)

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
        sp = n*1.9/max(self.nimages-1,1)+0.1
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
            #print(z[0], self.getShaderNo(z[0]))
            param[1][1] = self.getShaderNo(z[0])

    # display all t
    # you must register all images before display
    def openDisplay(self):
        self.prepare()
        CPlot.display(self.t, mode='Render', bgColor=1, posCam=(0.,0,1), 
                      posEye=(0.,0,0), dirCam=(0,1,0))
        CPlot.setState(activateShortCuts=0, displayInfo=0)

    # get noz of key
    def getNozOfKey(self, key, zones):
        for noz, z in enumerate(zones):
            if z[0] == key: return noz
        return -1

    # get noz of zpos
    def getNozOfZPos(self, zpos, zones):
        for noz, z in enumerate(zones):
            zp = Internal.getNodeFromName1(z, 'zpos')
            if zp is not None: zp = Internal.getValue(zp)
            else: zp = 0.
            if zp >= zpos: return noz
        return -1

    # move key right
    def moveRight(self, key, xr=0.01, render=True):
        self.move(key, (xr,0,0), render)

    def moveLeft(self, key, xr=0.01, render=True):
        self.move(key, (-xr,0,0), render)

    def moveUp(self, key, xr=0.01, render=True):
        self.move(key, (0,xr,0), render)

    def moveDown(self, key, xr=0.01, render=True):
        self.move(key, (0,-xr,0), render)

    # move a key from an offset v
    def move(self, key, v, render=True):
        if len(v) == 2: v = (v[0],v[1],0.)
        zones = Internal.getZones(self.t)
        noz = self.getNozOfKey(key, zones)
        z = Internal.getNodeFromName2(self.t, key)
        zp = T.translate(z, v)
        CPlot.replace(self.t,1,noz,zp)
        if render: CPlot.render()

    # place a key in x,y
    def place(self, key, pos, render=True):
        if len(pos) == 2: v = (pos[0],pos[1],0.)
        else: v = pos   
        zones = Internal.getZones(self.t)
        noz = self.getNozOfKey(key, zones)
        z = Internal.getNodeFromName2(self.t, key)
        sp = Internal.getNodeFromName2(z, 'ShaderParameters')[1]
        zp = D.point(v)
        CPlot._addRender2Zone(zp, material='Sphere', shaderParameters=[sp[0],sp[1]])
        zp[0] = key
        C._initVars(zp, 'TBB__', 0.)
        Internal._createChild(zp, 'zpos', 'DataArray_t', value=v[2])
        CPlot.replace(self.t,1,noz,zp)
        if render: CPlot.render()

    # Return zone of key in tree
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
    def drawImage(self, key, imageKey=None, pos=(0,0,0), scale=1., render=True):
        if imageKey is None: imageKey = key
        a = D.point(pos)
        CPlot._addRender2Zone(a, material='Sphere', shaderParameters=[30.*scale,self.getShaderNo(imageKey)])
        a[0] = key
        C._initVars(a, 'TBB__', 0.)
        Internal._createChild(a, 'zpos', 'DataArray_t', value=pos[2])
        # tri suivant zpos
        zones = Internal.getZones(self.t)
        noz = self.getNozOfZPos(pos[2], zones)
        CPlot.add(self.t, 1, noz, a)
        if not self.firstFit: self.fitView()
        if render: CPlot.render()

    # draw multiple time the same image in one key
    def drawMultipleImage(self, key, imageKey=None, pos=[(0,0,0)], scale=1., render=True):
        if imageKey is None: imageKey = key
        #a = D.polyline(pos)
        #a = C.convertArray2Hexa(a)
        #a = C.convertArray2Node(a)
        a = []
        for i in pos: a.append(D.point(i))
        a =  T.join(a)

        CPlot._addRender2Zone(a, material='Sphere', shaderParameters=[30.*scale,self.getShaderNo(imageKey)])
        a[0] = key
        C._initVars(a, 'TBB__', 0.)
        Internal._createChild(a, 'zpos', 'DataArray_t', value=pos[0][2])
        # tri suivant zpos
        zones = Internal.getZones(self.t)
        noz = self.getNozOfZPos(pos[0][2], zones)
        CPlot.add(self.t, 1, noz, a)
        if not self.firstFit: self.fitView()
        if render: CPlot.render()        

    # draw a particle system from a point
    def drawParticles(self, key, imageKey=None, pos=(0,0,0), Np=10, scale=1., render=True):
        import random
        if imageKey is None: imageKey = key
        out = []
        for i in range(Np):
            xpos = (pos[0]-i*0.01/Np,pos[1],pos[2])
            b = D.point(xpos)
            alpha = i*math.pi/(Np-1)
            C._initVars(b, 'vx=%g'%math.cos(alpha))
            C._initVars(b, 'vy=%g'%math.sin(alpha))
            C._initVars(b, 'vz=0.')
            r = random.randint(-10, 10)
            if r > 0: r = 0
            C._initVars(b, 'life=%d'%r)
            out.append(b)
        a = T.join(out, tol=1.e-12)
        #CPlot._addRender2Zone(a, material='Sphere', shaderParameters=[30.*scale,self.getShaderNo(imageKey)])
        a[0] = key
        C._initVars(a, 'TBB__', 0.)
        Internal._createChild(a, 'zpos', 'DataArray_t', value=pos[2])
        # Tri suivant zpos
        zones = Internal.getZones(self.t)
        noz = self.getNozOfZPos(pos[2], zones)
        CPlot.add(self.t, 1, noz, a)
        if not self.firstFit: self.fitView()
        if render: CPlot.render()        

    # draw a text
    def drawText(self, key, pos=(0,0,0), scale=1., text="youpi", color='Black',
                 h=0., render=True):
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
        if not self.firstFit: self.fitView()
        if render: CPlot.render()

    # draw a selector
    def drawSelector(self, key, pos=(0,0,0), w=1., h=1., e=0.01, color='Black', render=True):
        x0 = pos[0]; y0 = pos[1]; z0 = pos[2]
        P0 = (x0,y0,z0); P1 = (x0+w,y0,z0); P2 = (x0+w,y0+h,z0); P3 = (x0,y0+h,z0)
        P4 = (x0+e,y0+e,z0); P5 = (x0+w-e,y0+e,z0); P6 = (x0+w-e,y0+h-e,z0); P7 = (x0+e,y0+h-e,z0)
        Q1 =  D.quadrangle(P0,P1,P5,P4)
        Q2 =  D.quadrangle(P5,P1,P2,P6)
        Q3 =  D.quadrangle(P7,P6,P2,P3)
        Q4 =  D.quadrangle(P0,P4,P7,P3)
        a = T.join([Q1,Q2,Q3,Q4])
        a[0] = key
        C._initVars(a, 'TBB__', 0.)
        CPlot._addRender2Zone(a, color=color)
        Internal._createChild(a, 'zpos', 'DataArray_t', value=pos[2])
        CPlot.add(self.t, 1, -1, a)
        if not self.firstFit: self.fitView()
        if render: CPlot.render()

    # erase a key
    def erase(self, key, render=True):
        # Replace par un point hors cadre
        a = D.point((-100,-100,0))
        zones = Internal.getZones(self.t)
        noz = self.getNozOfKey(key, zones)
        CPlot.replace(self.t, 1, noz, a)
        if render: CPlot.render()

    # Animate one step of an image (pour les images contenant plusieurs images)
    # si curMin et curMax fournis, loop entre curMin et curMax
    def animate(self, key, curMin=None, curMax=None, imageKey=None, render=True):
        if imageKey is None: imageKey = key
        a = self.getZone(key)
        im = self.images[imageKey]
        speed = im[4]
        cur = im[5]
        cur += speed
        if curMin is None:
            if cur > 1.: cur = 0. # loop
        else:
            if cur >= curMax: cur = curMin # loop
        im[5] = cur
        C._initVars(a, 'TBB__', cur)
        zones = Internal.getZones(self.t)
        noz = self.getNozOfKey(key, zones)
        CPlot.replace(self.t, 1, noz, a)
        if render: CPlot.render()

    # Positionne l'animation de l'image key sur une case (posi,posj)
    # avec 0 < posi < Ni et 0 < posj < Nj 
    def setAnimation(self, key, posi, posj, imageKey=None, render=True):
        if imageKey is None: imageKey = key
        a = self.getZone(key)
        im = self.images[imageKey]
        Ni = im[1]; Nj = im[2]
        cur = posi + Ni*posj
        cur = cur*1. / (Ni*Nj-0.5)
        im[5] = cur
        C._initVars(a, 'TBB__', cur)
        zones = Internal.getZones(self.t)
        noz = self.getNozOfKey(key, zones)
        CPlot.replace(self.t, 1, noz, a)
        if render: CPlot.render()

    # Animate particles
    def animateParticles(self, key, step, gravity, render=True):
        a = self.getZone(key)
        C._initVars(a, '{life}={life}+1')
        C._initVars(a, '{vy}={vy}-%g'%gravity)
        C._initVars(a, '{ux}=%g * {vx} * ({life}>0)'%step)
        C._initVars(a, '{uy}=%g * {vy} * ({life}>0)'%step)
        a = T.deform(a, ['ux','uy','vz'])
        zones = Internal.getZones(self.t)
        noz = self.getNozOfKey(key, zones)
        CPlot.replace(self.t, 1, noz, a)       
        if render: CPlot.render()

    # Clone une key (jamais teste)
    def clone(self, key, newKey):
        a = self.getZone(key)
        b = Internal.copyTree(a); b[0] = newKey
        im = self.images[key]
        self.images[newKey] = [im[0],im[1],im[2],self.nimages,im[4],im[5]]
        self.line += [im[0],im[1],im[2]]
        self.nimages += 1
        CPlot.add(self.t, 1, -1, b)  

    # forcefitView
    # A appeler eventuellement si openDisplay est appele sans aucune image
    # et si certains bitmaps ne s'affichent pas
    def fitView(self, pos=(0,0), render=True):
        CPlot.fitView()
        CPlot.setState(posCam=(pos[0],pos[1],1),posEye=(pos[0],pos[1],0),dirCam=(0,1,0))
        self.firstFit = True
        if render: CPlot.render()

    # Equivalent de time sleep
    def sleep(self, seconds):
        time.sleep(seconds)

    # - create keyState -
    def createKeyState(self):
        self.keyState = {'up':0, 'down':0, 'left':0, 'right':0, 'space':0, 'escape':0}

    # - getkeys -
    # Modifie le key state en fonction des touches pressees ou relachees
    def getKeys(self):
        """Modifies keyState depending on pressed keys."""
        state = self.keyState
        key = CPlot.getKeyboard()
        if len(key) == 0: return
        CPlot.resetKeyboard()
        for k in range(len(key)):
            try: v = ord(key[k])
            except: v = key
            if v == 1: state['up']=1
            elif v == 5: state['up']=0
            elif v == 2: state['down']=1
            elif v == 6: state['down']=0
            elif v == 3: state['left']=1
            elif v == 7: state['left']=0
            elif v == 4: state['right']=1
            elif v == 8: state['right']=0
            elif v == 32: state['space']=1
            elif v == 27: state['escape']=1

