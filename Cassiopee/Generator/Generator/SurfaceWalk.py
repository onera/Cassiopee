"""Surface walk module. Extension of Generator.
"""
from . import Generator as G
from . import generator
__version__ = G.__version__

try: range = xrange
except: pass

#=============================================================================
# Python Interface to create surface grids by marching on surfaces
#=============================================================================
try: 
    import Transform as T
    import Converter as C
    import Post as P
    import Geom as D
    import math
except: raise ImportError("SurfaceWalk: requires Converter, Geom, Transform, Post modules.")

#=============================================================================
# build an extension with respect to normals to surfaces (smoothed niter times)
# starting from a curve c on the surface
#=============================================================================
def buildExtension__(c, surfaces, dh, niter=0):
    vars = c[0]
    c = C.convertBAR2Struct(c)
    imax = c[1].shape[1]
    for nos in range(len(surfaces)):
        if len(surfaces[nos]) == 5: surfaces[nos] = C.convertArray2Hexa(surfaces[nos])
    surfaces = T.join(surfaces)
    surfaces = G.close(surfaces)
    normals = G.getSmoothNormalMap(surfaces, niter=niter)
    surfaces2 = C.addVars([surfaces,normals])
    res = P.extractMesh([surfaces2], c)
    res = C.normalize(res, ['sx','sy','sz'])
    n1 = C.extractVars(res, ['sx','sy','sz'])
    # contour
    jmax = dh[2]
    cp = c
    coords = C.array(vars,imax,jmax,1)
    coords[1][0,0:imax] = cp[1][0,:]
    coords[1][1,0:imax] = cp[1][1,:]
    coords[1][2,0:imax] = cp[1][2,:]
    sn = ['sx','sy','sz']
    for j1 in range(1,jmax):
        hloc = dh[1][0,j1]-dh[1][0,j1-1]
        istart = j1*imax; iend = istart+imax
        ht = G.getLocalStepFactor__(cp, sn, smoothType=0, nitLocal=0, kappaType=0, kappaS=0, algo=0)
        ht[1][0,:] = hloc*ht[1][0,:]
        n2 = C.copy(n1)
        n2[1][0,:] = ht[1][0,:]*n1[1][0,:]
        n2[1][1,:] = ht[1][0,:]*n1[1][1,:]
        n2[1][2,:] = ht[1][0,:]*n1[1][2,:]
        cp = C.addVars([cp,n2])
        cn = T.deform(cp, ['sx','sy','sz'])
        coords[1][0,istart:iend] = cn[1][0,:]
        coords[1][1,istart:iend] = cn[1][1,:]
        coords[1][2,istart:iend] = cn[1][2,:]
        cp = cn
    return coords

def surfaceWalk__(surfaces, c, distrib, constraints, niter,alphaRef, check, toldist):
    tolsurf = (0.1*toldist)**2
    vn = ['sx','sy','sz'] # normales a la surface
    veta = ['etax','etay','etaz'] # tangente a la surface

    if len(distrib) != 5: distrib = C.convertBAR2Struct(distrib)
    c0 = G.close(c, tol=toldist)
    if len(c) == 4: c0 = C.convertBAR2Struct(c0) 

    surfaces = C.convertArray2Tetra(surfaces, split='withBarycenters')
    surfaces = T.join(surfaces)
    surfaces = G.close(surfaces)
    # 1. Projection ortho du contour sur surfaces
    c2 = T.projectOrtho(c0, [surfaces]); c2 = G.close(c2, tol=toldist)

    #2. Le contour est il une boucle ?
    jmax = distrib[2]; imax = c2[1].shape[1]; loop = 0
    dx = c2[1][0,0]-c2[1][0,imax-1]
    dy = c2[1][1,0]-c2[1][1,imax-1]
    dz = c2[1][2,0]-c2[1][2,imax-1]
    if (abs(dx) < toldist and abs(dy) < toldist and abs(dz) <toldist): loop = 1

    # 3. Si contraintes : determination du pas le plus petit pour redistribuer les pts sur les contraintes
    #                     suivie de la determination des points de contraintes
    constraints2 = []; constrainedPts = []
    if constraints != []:
        hp0 = 1.e12
        for i in range(1,distrib[2]): hp0 = min(hp0, distrib[1][0,i]-distrib[1][0,i-1])
        hook = C.createHook(c2, function='nodes')
        # identification des pts de la contrainte avec le contour
        for noc in range(len(constraints)):
            cons = constraints[noc]
            if len(cons) == 4: cons = C.convertBAR2Struct(cons)
            nodesc = C.identifyNodes(hook, cons, tol=toldist)

            indc = -1
            if nodesc[0] != -1: indc = nodesc[0]
            elif nodesc[-1] != -1: indc = nodesc[-1]; cons = T.reorder(cons,(-1,2,3))
            if indc != -1:
                # redistribution des contraintes
                L1 = D.getLength(cons); hp = hp0/L1 # ramene a [0,1]
                nds = int(1./hp)+1
                ds = G.cart((0,0,0),(hp,1,1),(nds,1,1))
                if ds[2] > 1: constraints2.append(G.map(cons, ds))
                else: constraints2.append(cons)
                constrainedPts.append(int(indc-1))# demarre a 1

        # Free hook
        C.freeHook(hook)

    # 4.  Surface walk
    # 4.1 tableau de coordonnees de la surface : la premiere  rangee correspond a c0 (et non c0 projete)
    coords = C.array(c0[0], imax, jmax,1)
    coords[1][:,0:imax] = c0[1][:,0:imax]
    if loop == 1: coords[1][:,imax-1] = coords[1][:,0]

    # 4.2 calcul des normales
    normales = G.getSmoothNormalMap(surfaces,niter=niter)
    normales = C.addVars([surfaces,normales])
    hook = C.createHook(normales, function='nodes')
    # Check de l angle de courbure max autorise
    cosalphaRef = math.cos(alphaRef*math.pi/180.)
    alphaMin = abs(alphaRef); alphaCheck = 1
    if abs(alphaRef) >= 180.: alphaCheck = 0
    stop = 0; j1 = 1; jmaxout = -1
    while j1 < jmax and stop == 0:
        # Interpolation des normales sur le contour
        indices, dists = C.nearestNodes(hook, c2) 
        c2 = C.addVars(c2, vn)
        alpn2 = C.extractVars(c2,vn)
        n1 = C.extractVars(normales,vn)
        alpn21 = alpn2[1]; n11 = n1[1]
        alpn21[:,:] = n11[:,indices[:]-1]

        alpn2 = C.normalize(alpn2, vn)
        alpn = alpn2
        # Calcul de eta sur le contour : ksi  x n2
        eta = generator.computeEta(c2, alpn, niter, loop)
        if j1 == 1: alpp = alpn

        if alphaCheck == 1 and j1 != 1:
            alpxn = alpn[1][0,:]; alpyn = alpn[1][1,:]; alpzn = alpn[1][2,:]
            alpxp = alpp[1][0,:]; alpyp = alpp[1][1,:]; alpzp = alpp[1][2,:]

            for ieta in range(eta[1].shape[1]):
                ps = alpxn[ieta]*alpxp[ieta]+alpyn[ieta]*alpyp[ieta]+alpzn[ieta]*alpzp[ieta]                
                if abs(ps) < cosalphaRef: stop = 1; jmaxout = j1; break
        #----------------------------------
        # add layer a partir de c et de eta
        #----------------------------------
        if stop == 0:
            alpp = alpn
            etap = C.normalize(eta, veta)
            if constraints2 != []: 
                eta = generator.straightenVector(c2, etap, constrainedPts, constraints2, loop, niter, toldist)
            else: eta = etap
            if loop == 1: # meme eta en imin et imax
                eta1 = (eta[1][0,0]+eta[1][0,imax-1])/2.
                eta2 = (eta[1][1,0]+eta[1][1,imax-1])/2.
                eta3 = (eta[1][2,0]+eta[1][2,imax-1])/2.
                eta[1][0,0] = eta1; eta[1][0,imax-1] = eta1
                eta[1][1,0] = eta2; eta[1][1,imax-1] = eta2
                eta[1][2,0] = eta3; eta[1][2,imax-1] = eta3
            # calcul du h eventuellement corrige
            hloc = distrib[1][0,j1]-distrib[1][0,j1-1]
            eta1 = eta[1]
            eta1[0,:] = hloc*eta1[0,:]
            eta1[1,:] = hloc*eta1[1,:]
            eta1[2,:] = hloc*eta1[2,:]
            c2 = C.addVars([c2,eta])
            c2 = T.deform(c2, veta)

            # projection orthogonale de c2 sur surfaces
            c2 = T.projectOrtho(c2, [surfaces])

            nim = j1*imax; nip2 = nim + imax
            coords[1][0,nim:nip2] = c2[1][0,0:imax]
            coords[1][1,nim:nip2] = c2[1][1,0:imax]
            coords[1][2,nim:nip2] = c2[1][2,0:imax]
            j1 += 1
            jmaxout = j1

    C.freeHook(hook)

    # Retourne les coordonnees de la surface creee
    j1prev = 1
    volmin = -1.e6
    if check == 1:
        jminout = jmaxout
        for j1 in range(1,jmaxout-1): 
            subc = T.subzone(coords,(1,j1,1),(imax,j1+1,1))            
            vol = G.getVolumeMap(subc)
            if C.getMinValue(vol,'vol') < 0.1*volmin: 
                jminout = min(jminout,j1prev)
                break
            if j1 == 1: volmin = C.getMinValue(vol,'vol')
            j1prev += 1
        return T.subzone(coords,(1,1,1),(imax,jminout,1))
    else: #check == 0: 
        if jmaxout == jmax: return coords
        else: return T.subzone(coords,(1,1,1),(imax,jmaxout,1))
    return None
