# ice mesher 2D
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Post.PyTree as P
import Geom.PyTree as D
import KCore.Vector as Vector
import Converter.Internal as Internal
import numpy

#===============================================================================================
# verifie l'order de profile (1D STRUCT)
# retourne 1 si normales exterieures et 0 sinon
#===============================================================================================
def checkOrder(profile):
    distrib = D.line((0,0,0), (0.01,0,0), N=2)
    lay = G.addNormalLayers(profile, distrib, niter=0, check=False)
    lay1 = T.subzone(lay, (1,-1,1), (-1,-1,1))
    l0 = D.getLength(profile)
    l1 = D.getLength(lay1)
    if l1 > l0: return 1
    else: return 0

#===============================================================================================
# IN: BAR1: profile with ice, correctly meshed, in XY plane, simple loop, no branching
# IN: BAR3: outer domain boundary, correctly meshed, in XY plane, simple curve with no branching
# IN: ht: total height of boundary layer mesh
# IN: hf: height of first wall cell
# OUT: mesh: ready for coda
#===============================================================================================
def mesher(BAR1, BAR3, ht, hf, extruder=0):

    #===============
    # BAR1 (profile)
    #===============
    BAR1 = C.convertArray2Tetra(BAR1)
    profile = C.initVars(BAR1, 'CoordinateZ', 0.) # force xy
    profile = C.convertBAR2Struct(profile) # profile is STRUCT
    # check order
    order = checkOrder(profile)
    print('order of profile (BAR1)', order)
    # reverse order
    if order == 0: T._reorder(profile, (-1,2,3))
    C.convertPyTree2File(profile, 'profileOrdered.cgns')

    #================
    # BAR3 (exterior)
    #================
    BAR3 = C.convertBAR2Struct(BAR3)
    BAR3 = C.initVars(BAR3, 'CoordinateZ', 0.) # force xy
    order = checkOrder(BAR3)
    print('order if exterior (BAR3)', order)
    if order == 0: T._reorder(BAR3, (-1,2,3))
    BAR3 = C.convertArray2Tetra(BAR3)
    C.convertPyTree2File(BAR3, 'exterior.cgns')

    #====================
    # Boundary layer mesh - avoid addNL fail a automatiser
    #====================
    G._getVolumeMap(profile)
    hmean = C.getMeanValue(profile, 'centers:vol')
    nLayers = int(ht/hmean)
    nLayers = max(nLayers, 10)
    nLayers = min(nLayers, 100)

    # choice extruder
    if extruder == 0:
        # version addNormalLayer
        print("extruder addNormalLayers")
        # number of iteration of normal smoothing
        niter = 100
        distrib = D.line((0,0,0), (ht,0,0), N=nLayers)
        goOn = True
        while goOn:
            bl = G.addNormalLayers(profile, distrib, niter=niter, check=True)
            zoneDim = Internal.getZoneDim(bl)
            nk = zoneDim[2]
            if nk < nLayers:
                print('addNL was stopped %d out of %d'%(nk, nLayers))
                niter += 200; #nLayers += 10
                if niter > 1000: goOn = False
            else: goOn = False
        T._reorder(bl, (-1,2,3))
    else:
        # version hyper2D
        print("extruder hyper2D")
        zoneDim = Internal.getZoneDim(profile)
        distrib = G.cart((0,0,0), (1,ht/(nLayers-1),1), (zoneDim[1],nLayers,1))
        profile2 = T.reorder(profile, (-1,2,3))
        bl = G.hyper2D(profile2, distrib, "O", etaStart=5, etaEnd=100, beta=0., forced=True)
    C.convertPyTree2File(bl, 'bl.cgns')

    # k-remeshing - weakness: no control of outer cell height
    #d = D.line((0,0,0), (1,0,0), N=40)
    #d = G.enforcePlusX(d, hf/ht, 40, 50)
    d = G.cartr2((0,0,0), (hf,1,1), (1.1,1.1,1.1), (ht,1,1))
    d = T.subzone(d, (1,1,1), (-1,1,1))
    C._initVars(d, '{CoordinateX}={CoordinateX}/%g'%ht)
    bl = G.map(bl, d, dir=2)
    C.convertPyTree2File(bl, 'bl2.cgns')

    # conversion en quad
    bl = C.convertArray2Hexa(bl)
    bl = G.close(bl)

    #=========================
    # exterior of bl mesh - OK
    #=========================
    ext = P.exteriorFaces(bl)
    ext = T.splitConnexity(ext)
    ext = ext[1] # always true
    T._reorder(ext, (-1,))

    #=============
    # T3 mesh - OK
    #=============
    borders = T.join(ext, BAR3)
    C.convertPyTree2File(borders, 'borders.cgns')
    m2 = G.T3mesher2D(borders, triangulateOnly=0, grading=1.1, metricInterpType=0)
    T._reorder(m2, (1,))

    # All mesh
    m = [bl, m2]

    # addkplane
    dz = 1.e-6
    T._addkplane(m)
    T._contract(m, (0,0,0), (1,0,0), (0,1,0), dz)

    # merge in one zone
    m2 = C.mergeConnectivity(m[0], m[1])

    #==========
    # BCs - OK
    #==========
    e0 = P.exteriorFaces(m)
    e0 = T.breakElements(e0)

    # recover BAR1
    e1 = T.addkplane(BAR1)
    T._contract(e1, (0,0,0), (1,0,0), (0,1,0), dz)
    e1[0] = 'wall'
    C._addBC2Zone(m2, 'wall', 'BCWall', subzone=e1)

    # recover BAR3
    e3 = T.addkplane(BAR3)
    T._contract(e3, (0,0,0), (1,0,0), (0,1,0), dz)
    e3[0] = 'far'
    C._addBC2Zone(m2, 'far', 'BCFarfield', subzone=e3)

    # recover symmetry planes
    es1 = P.selectCells(e0, '{CoordinateZ}<%g'%(0.5*dz), strict=True)
    es2 = P.selectCells(e0, '{CoordinateZ}>%g'%(dz-0.5*dz), strict=True)

    c = 0
    for e in es1+es2:
        if C.getNCells(e) > 0:
            e[0] = 'sym%d'%c
            C._addBC2Zone(m2, 'sym%d'%c, 'BCSymmetryPlane', subzone=e); c += 1

    return m2

#============================================================================
# IN: a: 1D struct or BAR (input geometry)
# IN: hx: mean h step to apply on profile
# IN: hmin/hmax: minimum and maximum h
# IN: fp: factor for peaks (pics)
# IN: ft: factor for throughs (creux)
# OUT: 1D struct remeshed, ready for extrusion and consSmooth
#============================================================================
def mapper(BAR1, hmin, hmax, hx, fp=0.5, ft=2.):

    BAR1 = C.convertArray2Tetra(BAR1)
    profile = C.initVars(BAR1, 'CoordinateZ', 0.) # force xy
    profile = C.convertBAR2Struct(profile) # profile is STRUCT
    # check order
    order = checkOrder(profile)
    if order == 0: T._reorder(profile, (-1,2,3))

    # closed?
    x0, y0, z0 = C.getValue(profile, 'GridCoordinates', 0)
    x1, y1, z1 = C.getValue(profile, 'GridCoordinates', -1)
    nn = Vector.norm( (x1-x0,y1-y0,z1-z0))
    closed = False
    if nn < 1.e-10: closed = True
    print("closed", closed)

    # split
    bs = T.splitSharpEdges(profile, alphaRef=30.)
    for c, b in enumerate(bs):
        bs[c] = C.convertBAR2Struct(b)

    # compute left and right vectors
    dxLeft = []; dxRight = []
    for b in bs:
        # left
        x0, y0, z0 = C.getValue(b, 'GridCoordinates', 0)
        x1, y1, z1 = C.getValue(b, 'GridCoordinates', 1)
        dxLeft.append((x1-x0,y1-y0,z1-z0))

        # right
        x0, y0, z0 = C.getValue(b, 'GridCoordinates', -1)
        x1, y1, z1 = C.getValue(b, 'GridCoordinates', -2)
        dxRight.append((x1-x0,y1-y0,z1-z0))

    # compute angle left
    angLeft = []
    l = len(bs)
    for c, b in enumerate(bs):
        if c == 0:
            dx1 = dxLeft[c]
            if closed: dx2 = dxRight[l-1]
            else: dx2 = (-dx1[0],-dx1[1],-dx1[2])
        else:
            dx1 = dxLeft[c]
            dx2 = dxRight[c-1]

        h1 = Vector.norm(dx1)
        h2 = Vector.norm(dx2)
        co = Vector.dot(dx1,dx2)
        co = co/(h1*h2)
        si = Vector.cross(dx1,dx2)
        si = si[2]/(h1*h2)
        asi = numpy.asin(si)
        aco = numpy.acos(co)
        if asi > 0: angle = aco # pic
        else: angle = -aco # creux
        angLeft.append(angle)

    # angle right
    angRight = []
    for c, b in enumerate(bs):
        if c < l-1: angle = angLeft[c+1]
        else:
            if closed: angle = angLeft[0]
            else: angle = numpy.pi
        angRight.append(angle)

    # creux: angle > 0, pointe: angle < 0
    for c, b in enumerate(bs):
        print(c, 'angleLeft=', angLeft[c]*180./numpy.pi, 'angleRight=',  angRight[c]*180./numpy.pi)

    # set h
    out = []
    critc = 100.*numpy.pi/180.
    critp = -100.*numpy.pi/180.
    for c, b in enumerate(bs):
        h1 = hx
        h2 = hx

        if angLeft[c] > 0 and angLeft[c] < critc: # creux critic
            h1 = h1*ft
        elif angLeft[c] < 0 and angLeft[c] > critp: # pic critic
            h1 = h1*fp
        if angRight[c] > 0 and angRight[c] < critc: # creux critic
            h2 = h2*ft
        elif angRight[c] < 0 and angRight[c] > critp: # pic critic
            h2 = h2*fp

        h1 = min(h1, hmax)
        h1 = max(h1, hmin)
        h2 = min(h2, hmax)
        h2 = max(h2, hmin)

        print(c, 'h1=', h1, 'h2=', h2)
        d = D.distrib2(b, h1, h2, normalized=True, algo=1)
        out.append(G.map(b, d))
    all = T.join(out)
    C.convertPyTree2File(all, 'profileMapped.cgns')
    return all

#==============================
# check quality of output mesh
#==============================
def checkQuality(m):
    # check positive volumes
    G._getVolumeMap(m)
    volmin = C.getMinValue(m, 'centers:vol')
    volmax = C.getMaxValue(m, 'centers:vol')
    volmean = C.getMeanValue(m, 'centers:vol')
    if volmin > 0: print("INFO: all positive volumes.")
    else: print("Warning: negative volume detected in final mesh.")

    return 0

#===============================================
# check deviation of profile from input profile
#===============================================
def checkDeviation(m):
    return 0.

def checkArea(m):
    return 0