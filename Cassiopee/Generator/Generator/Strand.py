# Extrusion for STRAND (TRI, pyTree, parallel, new version)
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Connector.Mpi as Xmpi
import Converter.Internal as Internal
import Generator.generator as generator
import Converter
import numpy

#===============================================
# IN: t: multibloc distributed TRI mesh
# IN: tc: connectivity of t
# IN: distrib: 1D mesh of heights (1D struct)
# IN: niter: number of smoothing iterations
# IN: eps: eps used in smoothing
#===============================================
def addNormalLayers(t, tc, distrib, niter=0, eps=0.4):

    kmax = C.getNPts(distrib)
    if kmax < 2: raise ValueError("addNormalLayers: distribution must contain at least 2 points.")

    # output tree
    strand = Internal.copyRef(t)
    coords = {} # all coordinates of planes
    for z in Internal.getZones(strand):
        zx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
        zy = Internal.getNodeFromName2(z, 'CoordinateY')[1]
        zz = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
        coords[z[0]] = []
        coords[z[0]].append((zx,zy,zz))

    # advancing surfaces
    surf = Internal.copyRef(t)

    for k1 in range(kmax-1):

        print("Generating layer %d"%k1)

        hloc = C.getValue(distrib, 'CoordinateX', k1+1) - C.getValue(distrib, 'CoordinateX', k1)

        # Keep pure normal in sx0
        surf = getSmoothNormalMap(surf, tc, niter=0, eps=eps)
        surf = modifyNormalWithMetric(surf, tc)
        C._initVars(surf, '{sx0} = {sx}')
        C._initVars(surf, '{sy0} = {sy}')
        C._initVars(surf, '{sz0} = {sz}')

        # Get smooth normal map
        surf = getSmoothNormalMap(surf, tc, niter=niter, eps=eps)

        # Modifiy sx with ht
        surf = modifyNormalWithMetric(surf, tc)

        # Mix sx0 and sx depending on height
        if kmax == 2: beta0 = 0.1
        else: beta0 = float((kmax-2-k1))/float(kmax-2); beta0 = beta0*beta0
        C._initVars(surf, '{sx} = %20.6g * {sx} + %20.6g * {sx0}'%(1-beta0, beta0))
        C._initVars(surf, '{sy} = %20.6g * {sy} + %20.6g * {sy0}'%(1-beta0, beta0))
        C._initVars(surf, '{sz} = %20.6g * {sz} + %20.6g * {sz0}'%(1-beta0, beta0))

        # multiply by the height
        C._initVars(surf, '{sx} = %20.6g * {sx}'%hloc)
        C._initVars(surf, '{sy} = %20.6g * {sy}'%hloc)
        C._initVars(surf, '{sz} = %20.6g * {sz}'%hloc)

        # get layers (prisms) - for volume check
        layers = G.grow(surf, vector=['sx','sy','sz'])
        G._getVolumeMap(layers)
        #C.convertPyTree2File(layers, "out.cgns")
        volmin = C.getMinValue(layers, "centers:vol")
        print("INFO: volmin of layer=", volmin)
        if volmin < 0: C.convertPyTree2File(layers, 'layers%03d.cgns'%k1)

        # deform surf to get new TRI surf
        T._deform(surf, vector=['sx','sy','sz'])

        # keep track of coordinates
        for z in Internal.getZones(surf):
            zx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
            zy = Internal.getNodeFromName2(z, 'CoordinateY')[1]
            zz = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
            coords[z[0]].append((zx,zy,zz))

    # build strand mesh
    for z in Internal.getZones(strand):
        coord = coords[z[0]]
        nk = len(coord)
        np = C.getNPts(z)
        fx = numpy.empty(np*nk, dtype=numpy.float64)
        fy = numpy.empty(np*nk, dtype=numpy.float64)
        fz = numpy.empty(np*nk, dtype=numpy.float64)
        for k, i in enumerate(coord):
            px = i[0]; py = i[1]; pz = i[2]
            fx[k*np:(k+1)*np] = px[:]
            fy[k*np:(k+1)*np] = py[:]
            fz[k*np:(k+1)*np] = pz[:]

        # Repush coordinates in strand grid
        x1 = Internal.getNodeFromName2(z, 'CoordinateX')
        y1 = Internal.getNodeFromName2(z, 'CoordinateY')
        z1 = Internal.getNodeFromName2(z, 'CoordinateZ')

        x1[1] = fx.ravel('k')
        y1[1] = fy.ravel('k')
        z1[1] = fz.ravel('k')

        # modify strand grid dimensions
        dim = z[1]; dim[0] = fx.size

    return strand

#===============================================
# IN: t: multibloc distributed TRI mesh
# IN: tc: connectivity of t
# IN: distrib: 1D mesh of heights (1D struct)
# IN: niter: number of smoothing iterations
# IN: eps: eps used in smoothing
# this version uses blanking
#===============================================
def addNormalLayers2(t, tc, distrib, niter=0, eps=0.4):

    kmax = C.getNPts(distrib)
    if kmax < 2: raise ValueError("addNormalLayers: distribution must contain at least 2 points.")

    # output tree
    strand = Internal.copyRef(t)
    coords = {} # all coordinates of planes
    for z in Internal.getZones(strand):
        zx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
        zy = Internal.getNodeFromName2(z, 'CoordinateY')[1]
        zz = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
        coords[z[0]] = []
        coords[z[0]].append((zx,zy,zz))

    # advancing surfaces
    surf = Internal.copyRef(t)
    #C._initVars(surf, 'cellN', 1.)

    for k1 in range(kmax-1):

        print("Generating layer %d"%k1)

        hloc = C.getValue(distrib, 'CoordinateX', k1+1) - C.getValue(distrib, 'CoordinateX', k1)

        # Keep pure normal in sx0
        surf = getSmoothNormalMap(surf, tc, niter=0, eps=eps)
        surf = modifyNormalWithMetric(surf, tc)
        C._initVars(surf, '{sx0} = {sx}')
        C._initVars(surf, '{sy0} = {sy}')
        C._initVars(surf, '{sz0} = {sz}')

        # Get smooth normal map
        surf = getSmoothNormalMap(surf, tc, niter=niter, eps=eps)

        # Modifiy sx with ht
        surf = modifyNormalWithMetric(surf, tc)

        # Mix sx0 and sx depending on height
        #if kmax == 2: beta0 = 0.1
        #else: beta0 = float((kmax-2-k1))/float(kmax-2); beta0 = beta0*beta0
        #C._initVars(surf, '{sx} = %20.6g * {sx} + %20.6g * {sx0}'%(1-beta0, beta0))
        #C._initVars(surf, '{sy} = %20.6g * {sy} + %20.6g * {sy0}'%(1-beta0, beta0))
        #C._initVars(surf, '{sz} = %20.6g * {sz} + %20.6g * {sz0}'%(1-beta0, beta0))

        # multiply by the height
        C._initVars(surf, '{sx} = %20.6g * {sx}'%hloc)
        C._initVars(surf, '{sy} = %20.6g * {sy}'%hloc)
        C._initVars(surf, '{sz} = %20.6g * {sz}'%hloc)

        # get layers (prisms) - for volume check
        layers = G.grow(surf, vector=['sx','sy','sz'])
        G._getVolumeMap(layers)
        #C.convertPyTree2File(layers, "out.cgns")
        volmin = C.getMinValue(layers, "centers:vol")
        if volmin > 0: print("INFO: volmin of layer=", volmin)
        else: print("FAIL: volmin of layer=", volmin)
        if volmin <= 0:
            C._initVars(layers, '{centers:cellN}={centers:vol}>0')
            C.convertPyTree2File(layers, 'layers%03d.cgns'%k1)

        # deform surf to get new TRI surf
        T._deform(surf, vector=['sx','sy','sz'])

        # keep track of coordinates
        for z in Internal.getZones(surf):
            zx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
            zy = Internal.getNodeFromName2(z, 'CoordinateY')[1]
            zz = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
            coords[z[0]].append((zx,zy,zz))

    # build strand mesh
    for z in Internal.getZones(strand):
        coord = coords[z[0]]
        nk = len(coord)
        np = C.getNPts(z)
        fx = numpy.empty(np*nk, dtype=numpy.float64)
        fy = numpy.empty(np*nk, dtype=numpy.float64)
        fz = numpy.empty(np*nk, dtype=numpy.float64)
        for k, i in enumerate(coord):
            px = i[0]; py = i[1]; pz = i[2]
            fx[k*np:(k+1)*np] = px[:]
            fy[k*np:(k+1)*np] = py[:]
            fz[k*np:(k+1)*np] = pz[:]

        # Repush coordinates in strand grid
        x1 = Internal.getNodeFromName2(z, 'CoordinateX')
        y1 = Internal.getNodeFromName2(z, 'CoordinateY')
        z1 = Internal.getNodeFromName2(z, 'CoordinateZ')

        x1[1] = fx.ravel('k')
        y1[1] = fy.ravel('k')
        z1[1] = fz.ravel('k')

        # modify strand grid dimensions
        dim = z[1]; dim[0,0] = fx.size

    return strand

#===================================
# Compute a smooth normal map on t
#===================================
def getSmoothNormalMap(surf, tc, niter=2, eps=0.4):
    it = 1
    G._getNormalMap(surf)
    C._normalize(surf, ['centers:sx','centers:sy','centers:sz'])
    surf = C.center2Node(surf, 'centers:sx')
    surf = C.center2Node(surf, 'centers:sy')
    surf = C.center2Node(surf, 'centers:sz')
    C._normalize(surf, ['sx','sy','sz'])
    #Xmpi._setInterpTransfers(surf, tc, variables=['sx','sy','sz'], storage=1, compact=0))

    while it < niter:
        surf = C.node2Center(surf, 'sx')
        surf = C.node2Center(surf, 'sy')
        surf = C.node2Center(surf, 'sz')
        C._normalize(surf, ['centers:sx','centers:sy','centers:sz'])
        surf = C.center2Node(surf, 'centers:sx')
        surf = C.center2Node(surf, 'centers:sy')
        surf = C.center2Node(surf, 'centers:sz')
        C._normalize(surf, ['sx','sy','sz'])
        #Xmpi._setInterpTransfers(surf, tc, variables=['sx','sy','sz'], compact=0)
        it += 1
    return surf

#==========================================
# modify the normal with local step factor
#==========================================
def modifyNormalWithMetric(surf, tc):
    surf = getLocalStepFactor(surf, tc)
    C._initVars(surf, '{sx} = {ht}*{sx}')
    C._initVars(surf, '{sy} = {ht}*{sy}')
    C._initVars(surf, '{sz} = {ht}*{sz}')
    return surf

#===========================
# get the local step factor
#===========================
def getLocalStepFactor(surf, tc):

    for z in Internal.getZones(surf):
        s = C.getFields(Internal.__GridCoordinates__, z, api=1)[0]
        sn = C.getFields(Internal.__FlowSolutionNodes__, z, api=1)[0]
        sc = Converter.addVars([s,sn])
        sc = Converter.node2Center(sc)
        ht = generator.getLocalStepFactor(s, sc)
        C.setFields([ht], z, 'nodes', writeDim=False)

    # Lissage hauteur du pas
    niter = 10; it = 0
    while it < niter:
        surf = C.node2Center(surf, 'ht')
        surf = C.center2Node(surf, 'centers:ht')
        Xmpi._setInterpTransfers(surf, tc, variables=['ht'], compact=0)
        it += 1

    return surf

def _remap(t, d):
    for z in Internal.getZones(t):
        _remap__(z, d)
    return None

# remap strand zone in k
def _remap__(a, d):
    # Get npxnk from strand grid
    dim = Internal.getZoneDim(a)
    cn = Internal.getNodeFromName2(a, 'ElementConnectivity')
    np = numpy.max(cn[1])
    nk = dim[1]//np
    x1 = Internal.getNodeFromName2(a, 'CoordinateX')[1]
    y1 = Internal.getNodeFromName2(a, 'CoordinateY')[1]
    z1 = Internal.getNodeFromName2(a, 'CoordinateZ')[1]

    # create a struct 2D containing strands
    b = G.cart((0,0,0), (1,1,1), (np,nk,1))
    x = Internal.getNodeFromName2(b, 'CoordinateX')[1]
    y = Internal.getNodeFromName2(b, 'CoordinateY')[1]
    z = Internal.getNodeFromName2(b, 'CoordinateZ')[1]
    x.ravel('k')[:] = x1[:]
    y.ravel('k')[:] = y1[:]
    z.ravel('k')[:] = z1[:]
    #C.convertPyTree2File(b, 'out.cgns')

    # Remap the struct grid
    #h = 0.06
    #dy = 0.0002
    #d = D.line((0,0,0), (1,0,0), N=40)
    #d = G.enforcePlusX(d, dy/h, 40, 50)
    b = G.map(b, d, dir=2)
    #C.convertPyTree2File(b, 'out.cgns')

    x = Internal.getNodeFromName2(b, 'CoordinateX')[1]
    y = Internal.getNodeFromName2(b, 'CoordinateY')[1]
    z = Internal.getNodeFromName2(b, 'CoordinateZ')[1]

    # Repush coordinates in strand grid
    x1 = Internal.getNodeFromName2(a, 'CoordinateX')
    y1 = Internal.getNodeFromName2(a, 'CoordinateY')
    z1 = Internal.getNodeFromName2(a, 'CoordinateZ')

    x1[1] = x.ravel('k')
    y1[1] = y.ravel('k')
    z1[1] = z.ravel('k')

    # modify strand grid dimensions
    dim = a[1]; dim[0] = x.size

    return None

# subzone strand grid (keep only a given k plane)
def subzonePlane(t, k):
    tp = Internal.copyRef(t)
    for z in Internal.getZones(tp):
        _subzonePlane__(z, k)
    return tp

def _subzonePlane__(a, k):
    dim = Internal.getZoneDim(a)
    cn = Internal.getNodeFromName2(a, 'ElementConnectivity')
    np = numpy.max(cn[1])
    nk = dim[1]//np

    x1 = Internal.getNodeFromName2(a, 'CoordinateX')[1]
    y1 = Internal.getNodeFromName2(a, 'CoordinateY')[1]
    z1 = Internal.getNodeFromName2(a, 'CoordinateZ')[1]
    x2 = x1[np*k: np*(k+1)]
    y2 = y1[np*k: np*(k+1)]
    z2 = z1[np*k: np*(k+1)]

    x1 = Internal.getNodeFromName2(a, 'CoordinateX')
    y1 = Internal.getNodeFromName2(a, 'CoordinateY')
    z1 = Internal.getNodeFromName2(a, 'CoordinateZ')

    x1[1] = x2.ravel('k')
    y1[1] = y2.ravel('k')
    z1[1] = z2.ravel('k')

    # modify strand grid dimensions
    dim = a[1]; dim[0] = x2.size
    return None

# rm ghost cells from strand grid
def rmGhostCells(a):
    # a faire
    return None