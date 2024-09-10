"""Mapedge functions."""
try:
    import Generator as G
    import Transform as T
    import Converter as C
    import Geom as D
    import KCore.kcore as KCore
    import numpy
    importOK = True
except ImportError: 
    importOK = False

try: range = xrange
except: pass

def checkImport():
    if not importOK: 
        raise ImportError("mapEdge requires Converter, Generator, Transform.")

# uniformize a 1D mesh (regular step) - OK
def uniformize(a, N=100, h=-1., factor=-1, density=-1., sharpAngle=30.):
    """Uniformize the distribution of points on a 1D-mesh."""
    if checkImport is None: return None
    if isinstance(a[0], list):
        L = D.getLength(a)
        if density > 0: N = int(L*density)+1
        elif factor > 0: N = int(factor*(C.getNPts(a)-1))+1
        elif h > 0.: N = int(L/h)+1
        b = []; Nt = 0; rest = 0
        for ct, i in enumerate(a):
            Li = D.getLength(i)
            Ni = int(T.kround(N*1.*(Li/L)))+1
            rest += Ni-1-N*1.*(Li/L)
            if rest < -0.9: Ni += 1; rest = 0.
            elif rest > 0.9: Ni -= 1; rest = 0.
            Ni = max(Ni,2)
            Nt += Ni-1
            if ct == len(a)-1: Ni = Ni-Nt+N-1; Ni = max(Ni,2)
            b.append(uniformize__(i, Ni, -1, -1, -1, sharpAngle))
        return b
    else:
        return uniformize__(a, N, h, factor, density, sharpAngle)

def uniformize__(a, N, h, factor, density, sharpAngle):
    if len(a) == 4: ntype = 1
    else: ntype = 0
    L = D.getLength(a)
    if density > 0: N = int(L*density)+1
    elif factor > 0: N = int(factor*(C.getNPts(a)-1))+1
    elif h > 0: N = int(L/h)+1
    N = max(N, 2)
    
    # I can add splitConnexity if needed
    b = T.splitSharpEdges(a, alphaRef=sharpAngle) # b is BAR closed
    if ntype == 1: b = T.splitTBranches(b)
    
    out = []; Nt = 0; rest = 0
    for ct, i in enumerate(b):
        if C.getNPts(i) >= 2: # avoid degenerated cases 
            Li = D.getLength(i)
            Ni = int(T.kround(N*1.*(Li/L)))+1
            rest += Ni-1-N*1.*(Li/L)
            if rest < -0.9: Ni += 1; rest = 0.
            elif rest > 0.9: Ni -= 1; rest = 0.
            Ni = max(Ni,2)
            Nt += Ni-1
            if ct == len(b)-1: Ni = Ni-Nt+N-1; Ni = max(Ni,2)
            i = C.convertBAR2Struct(i)
            distrib = G.cart((0,0,0), (1./(Ni-1.),1.,1.), (Ni,1,1))
            c = G.map(i, distrib) # c is STRUCT
            out.append(c)
    if ntype == 1: out = C.convertArray2Hexa(out)
    if len(out) > 1: b = T.join(out)
    else: b = out[0]
    if ntype == 1: b = G.close(b) # because of closed BARs
    return b

# Refine a 1D mesh (refine the initial distribution) - OK
def refine(a, N=10, factor=-1, sharpAngle=30.):
    """Refine the point distribution on a 1D-mesh."""
    if checkImport is None: return None
    if isinstance(a[0], list): 
        b = []
        for i in a:
            b.append(refine__(i, N, factor, sharpAngle))
        return b
    else:
        return refine__(a, N, factor, sharpAngle)

def refine__(a, N, factor, sharpAngle):
    if len(a) == 4: ntype = 1
    else: ntype = 0
    NPts = C.getNPts(a)
    L = D.getLength(a)
    N = max(N, 2)
    b = T.splitSharpEdges(a, alphaRef=sharpAngle)
    if ntype == 1: b = T.splitTBranches(b)

    out = []; Nt = 0; rest = 0
    for ct, i in enumerate(b):
        i = C.convertBAR2Struct(i)
        Li = D.getLength(i)
        if factor < 0.:
            Ni = int(T.kround(N*(C.getNPts(i)-1)/(NPts-1))+1)
            rest += Ni-1-N*1.*(Li/L)
            if rest < -0.9: Ni += 1; rest = 0.
            elif rest > 0.9: Ni -= 1; rest = 0.
            Ni = max(Ni,2)
            Nt += Ni-1
            if ct == len(b)-1: Ni = Ni-Nt+N-1
            factori = (Ni-1.)/(C.getNPts(i)-1.)
        else: factori = factor
        out.append(G.refine(i, factori))
    if ntype == 0:
        b = T.join(out)
    else:
        out = C.convertArray2Hexa(out)
        b = T.join(out)
        b = G.close(b)
    return b

# -- converti une chaine en range (a Structure seulement)
# IN: 'imin', 'jmin', ...
# OUT: [imin,imax,jmin,jmax,kmin,kmax]
def convertStringRange2Range__(wrange, a):
  ni = a[2]; nj = a[3]; nk = a[4]
  if wrange == 'imin': wrange = [1,1,1,nj,1,nk]
  elif wrange == 'imax': wrange = [ni,ni,1,nj,1,nk]
  elif wrange == 'jmin': wrange = [1,ni,1,1,1,nk]
  elif wrange == 'jmax': wrange = [1,ni,nj,nj,1,nk]
  elif wrange == 'kmin': wrange = [1,ni,1,nj,1,1]
  else: wrange = [1,ni,1,nj,nk,nk]
  return wrange

# set h step in place
# ind : indice ou [imin,imax,jmin,jmax,kmin,kmax] ou 'imin'
def setH(a, ind, h):
    """Set the mesh step for a given index."""
    pos = KCore.isNamePresent(a, 'h')
    if pos == -1:
        C._initVars(a, 'h', 0.)
        pos = KCore.isNamePresent(a, 'h')
    if isinstance(ind, int):
        if ind == -1: ind = a[1].shape[1]-1
        a[1][pos,ind] = h
    elif isinstance(ind, str):
        [imin,imax,jmin,jmax,kmin,kmax] = convertStringRange2Range__(ind, a)
        b = a[1][pos]
        b = b.reshape((a[4],a[3],a[2]))
        imin = imin-1; jmin = jmin-1; kmin = kmin-1
        b[kmin:kmax,jmin:jmax,imin:imax] = h
    else: # suppose range
        if len(ind) == 2:
            [imin,imax] = ind
            jmin = 1; jmax = 1; kmin = 1; kmax = 1
        elif len(ind) == 4:
            [imin,imax,jmin,jmax] = ind
            kmin = 1; kmax = 1
        else: [imin,imax,jmin,jmax,kmin,kmax] = ind
        b = a[1][pos]
        b = b.reshape((a[4],a[3],a[2]))
        imin = imin-1; jmin = jmin-1; kmin = kmin-1
        b[kmin:kmax,jmin:jmax,imin:imax] = h

# set factor in place
def setF(a, ind, f):
    """Set the mesh size factor for a given index."""
    pos = KCore.isNamePresent(a, 'f')
    if pos == -1:
        C._initVars(a, 'f', 0.)
        pos = KCore.isNamePresent(a, 'f')
    if isinstance(ind, int):
        if ind == -1: ind = a[1].shape[1]-1
        a[1][pos,ind] = f
    elif isinstance(ind, str):
        [imin,imax,jmin,jmax,kmin,kmax] = convertStringRange2Range__(ind, a)
        b = a[1][pos]
        b = b.reshape((a[4],a[3],a[2]))
        imin = imin-1; jmin = jmin-1; kmin = kmin-1
        b[kmin:kmax,jmin:jmax,imin:imax] = f
    else: # suppose range
        if len(ind) == 2:
           [imin,imax] = ind
           jmin = 1; jmax = 1; kmin = 1; kmax = 1
        elif len(ind) == 4:
           [imin,imax,jmin,jmax] = ind
           kmin = 1; kmax = 1
        else: [imin,imax,jmin,jmax,kmin,kmax] = ind
        b = a[1][pos]
        b = b.reshape((a[4],a[3],a[2]))
        imin = imin-1; jmin = jmin-1; kmin = kmin-1
        b[kmin:kmax,jmin:jmax,imin:imax] = f

# Build a distrib between 0. and 1. with h1 left, h2 right with N-1 points
def buildDistrib(h1, h2, N):
    Ni = int(T.kround(1./h2))+1
    a = G.cart((0,0,0), (h2,1,1), (Ni,1,1))
    a[1][0,Ni-2] = 1.-h2
    a[1][0,Ni-1] = 1.
    b = G.enforcePlusX(a, h1, (Ni-2,N-Ni-1), verbose=False)
    return b

# moyenne ponderee de hl
def moyenne(a, hl):
    L = D.getLength(a)
    npts = hl.size
    h1 = -1; h2 = -1; i1 = -1; i2 = -1
    href = 0.
    for i in range(npts):
        hi = hl[i]
        if hi == 0. and i == npts-1: hi = h1
        if hi > 1.e-12:
            if i == 0: i1 = 0; h1 = hi
            if i > 0:
                if i1 == -1:
                    i1 = 0; i2 = i; h1 = hi; h2 = hi
                if h1 > 0:
                    i2 = i; h2 = hi
                sub = T.subzone(a, (i1+1,1,1), (i2+1,1,1))                
                Li = D.getLength(sub)
                href += (h1+h2)*0.5*Li
                # loop
                i1 = i2; h1 = h2
    href = href / L
    return href

def enforceh3D(a, N=100, h=-1, dir=1):
    """Enforce mesh size in a 3D-mesh."""
    if checkImport is None: return None
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(enforceh3D_(i, N, h, dir))
        return b
    else:
        return enforceh3D_(a, N, h, dir)

def enforceh3D_(array, N, h, dir):
    if dir == 2: m = T.reorder(array, (2,1,3))
    elif dir == 3: m = T.reorder(array, (3,2,1))
    elif dir == 1: m = array
    ni = m[2]; nj = m[3]; nk = m[4]
    # first line
    l = T.subzone(m, (1,1,1), (ni,1,1))
    am = enforceh(l, N, h)
    ndi = am[2]
    a = C.array('x,y,z', ndi, nj, nk)
    for k in range(nk):
        for j in range(nj):
            l = T.subzone(m, (1,j+1,k+1), (ni,j+1,k+1))
            am = enforceh(l, N)
            ind = j*ndi+k*ndi*nj
            a[1][:,ind:ndi+ind] = am[1][:,0:ndi]
    if dir == 2: a = T.reorder(a, (2,1,3))
    elif dir == 3: a = T.reorder(a, (3,2,1))
    return a

# h must be set at sharp angles
# Set: N + h field (step size)
# Or: h + f field (factor)
def enforceh(a, N=100, h=-1):
    """Enforce mesh size in a 1D-mesh."""
    if checkImport is None: return None
    if isinstance(a[0], list):
        L = D.getLength(a)
        b = []
        for i in a:
            Li = D.getLength(i)
            Ni = int(T.kround(N*1.*(Li/L)))+1
            b.append(enforceh_(i, Ni, h))
        return b
    else:
        return enforceh_(a, N, h)

def enforceh_(a, N, h):
    if len(a) == 4: ntype = 1
    else: ntype = 0
    L = D.getLength(a)
    if ntype == 1:
        a = T.splitTBranches(a)
        a = C.convertBAR2Struct(a)
        out = []
        for i in a:
            Li = D.getLength(i)
            Ni = Li/L*N
            out.append(enforceh__(i, Ni, h))
        out = C.convertArray2Hexa(out)
        out = T.join(out)
    else:
        out = enforceh__(a, N, h)
    return out

def enforceh__(a, N, h):
    pos = KCore.isNamePresent(a, 'h')
    if pos == -1:
        pos = KCore.isNamePresent(a, 'f')
        if pos == -1:
            #raise ValueError("enforceh: h or f is not present.")
            return a
        else: 
            hl = numpy.copy(a[1][pos])
            factorMoy = moyenne(a, hl)
            if h > 0:
                href = h
                N = D.getLength(a)/(factorMoy*href)+1
                N = int(T.kround(N))
                hl[:] = hl[:]*href
            else:
                href = D.getLength(a)/((N-1.)*factorMoy)
                hl[:] = hl[:]*href
    else: hl = a[1][pos]
    npts = hl.size
    # Calcul h1s, h2s, i1s, i2s
    h1s=[]; h2s=[]; Ps=[]; i1s=[]; i2s=[]
    h1 = -1; h2 = -1; i1 = -1; i2 = -1
    
    # compute h
    vol = G.getVolumeMap(a)
    vol = C.center2Node(vol)
    vol = vol[1].ravel('k')

    # introduce delta h step (input:delta)
    #firstH = -1
    #for i in range(0, npts):
    #    hi = hl[i]
    #    if hi > 1.e-12: firstH = i; break
    #lastH = -1
    #for i in range(npts-1, -1, -1):
    #    hi = hl[i]
    #    if hi > 1.e-12: lastH = i; break
    #a1 = None; a2 = None
    #if firstH > delta and lastH < npts-delta:
    #    a1 = T.subzone(a, (1,1,1), (firstH-delta,-1,-1))
    #    a2 = T.subzone(a, (lastH+delta,1,1), (-1,-1,-1))
    #    a = T.subzone(a, (firstH-delta,1,1), (lastH+delta,-1,-1))
    #elif firstH > delta:
    #    a1 = T.subzone(a, (1,1,1), (firstH-delta,-1,-1))
    #    a = T.subzone(a, (firstH-delta,1,1), (-1,-1,-1))
    #elif lastH < npts-delta:
    #    a = T.subzone(a, (firstH-delta,1,1), (-1,-1,-1))
    #    a2 = T.subzone(a, (lastH+delta,1,1), (-1,-1,-1))

    for i in range(npts):
        hi = hl[i]
        if i == 0: # set first h1
            i1 = 0; h1 = hi
            if hi == 0.: h1 = vol[0]
            continue
        if hi == 0. and i == npts-1:
            hi = h1 # extrapolation
            hi = vol[npts-1]
        
        if hi > 1.e-12:
            #if i1 == -1: i1 = 0; i2 = i; h1 = hi; h2 = hi # propagate hi from inside if not set on borders
            # keep border h if not set
            i2 = i; h2 = hi
            sub = T.subzone(a, (i1+1,1,1), (i2+1,1,1))             
            Li = D.getLength(sub)
            Pi = Li*1./(0.5*(h1+h2)+1.e-12)
            i1s.append(i1); i2s.append(i2)
            h1s.append(h1/Li); h2s.append(h2/Li)
            Ps.append(Pi)
            i1 = i2; h1 = h2

    Pt = 0.
    for x in range(len(h1s)):
        Pi = Ps[x]
        Pt += Pi

    for x in range(len(h1s)):
        Pi = Ps[x]
        Pi = Pi/Pt*N
        Ps[x] = int(T.kround(Pi))+1

    out = []
    for x in range(len(h1s)):
        i1 = i1s[x]; i2 = i2s[x]
        h1 = h1s[x]; h2 = h2s[x]
        sub = T.subzone(a, (i1+1,1,1), (i2+1,1,1))
        d = buildDistrib(h1, h2, Ps[x])
        c = G.map(sub, d)
        #if h < -0.5: # pourquoi?
        #    setH(c, 0, sub[1][pos,0]); setH(c, -1, sub[1][pos,-1])
        out.append(c)
    out = T.join(out)
    out = C.rmVars(out, ['h'])
    return out

# Enforce h at ind (STRUCT)
def enforce(a, h, ind, supp, add):
    """Enforce h in a line. Return line."""
    c = D.getDistribution(a)
    L = D.getLength(a)
    if ind == 0: b = G.enforceMoinsX(c, h/L, supp, add)
    elif ind == a[2]-1: b = G.enforcePlusX(c, h/L, supp, add)
    else: b = G.enforceX(c, c[1][0,ind], h/L, supp, add)
    a = G.map(a, c)
    return a

# Enforce h everywhere in a line
def distrib1(a, h, normalized=True):
    """Enforce h in line. Return distribution."""
    L = D.getLength(a)
    N = int(T.kround(L / h))+1 # correspondant a hf
    if N <= 2: print('Error: distrib1: not enough point to remesh.')
    a = G.cart((0,0,0), (h,1,1), (N,1,1)) # match h
    if normalized: a = D.getDistribution(a)
    return a

# Enforce h1 and h2 at extremities in a dimensionned line
# algo: 0 (tanh), 1 (geom)
# Return distribution ready to map
def distrib2(a, h1, h2, add=20, forceAdd=False, normalized=True, algo=0, verbose=0):
    """Enforce h1,h2 in line. Return distribution."""
    L = D.getLength(a)

    if algo == 0: # tanh
        if h1 > h2:
            N = int(T.kround(L / h1))+1 # correspondant a hf
            if N <= 2: print('Error: distrib2: not enough point to remesh.')
            a = G.cart((0,0,0), (h1,1,1), (N,1,1)) # match h1
            if forceAdd: a = G.enforceMoinsX(a, h2, (N-1, add)) # match h2
            else: a = G.enforceMoinsX(a, h2, N-1, add) # match h2
            if normalized: a = D.getDistribution(a)
            return a
        else:
            N = int(T.kround(L / h2))+1 # correspondant a hf
            if N <= 2: print('Error: distrib2: not enough point to remesh.')
            a = G.cart((0,0,0), (h2,1,1), (N,1,1)) # match h2
            if forceAdd: a = G.enforcePlusX(a, h1, (N-1, add)) # match h1
            else: a = G.enforcePlusX(a, h1, N-1, add)
            if normalized: a = D.getDistribution(a)
            return a
    else: # geometrique
        if abs(h2-h1) < 1.e-6: # constant step
            N = int(T.kround(L / h1))
            N = max(N, 2)
            h = L/N
            a = G.cart((0,0,0), (h,1,1), (N,1,1))
        else:
            q = (L-h1)/(L-h2)
            if verbose > 0: print("Info: distrib2: geometric progression: %g"%q)
            N = numpy.log(h2/h1) / numpy.log(q)
            N = int(T.kround(N))+2
            #if N <= 2: print('Error: distrib2: not enough point to remesh.')
            N = max(N, 2)
            a = G.cartr1((0,0,0), (h1,1,1), (q,1,1), (N,1,1))
        if normalized: a = D.getDistribution(a)
        return a

# Super smooth - OK
def smooth(a, eps=0.1, niter=5, sharpAngle=30.):
    """Smooth point distribution in a 1D-mesh."""
    if checkImport is None: return None
    if isinstance(a[0], list): 
        b = []
        for i in a:
            b.append(smooth__(i, eps, niter, sharpAngle))
        return b
    else:
        return smooth__(a, eps, niter, sharpAngle)

def smooth__(a, eps, niter, sharpAngle):
    if len(a) == 4: ntype = 1
    else: ntype = 0
    b = T.splitSharpEdges(a, alphaRef=sharpAngle)
    if ntype == 1: b = T.splitTBranches(b)
    out = []
    P0 = D.point((0,0,0)); P1 = D.point((1,0,0))
    for i in b:
        i = C.convertBAR2Struct(i)
        c = D.getDistribution(i)
        d = T.smooth(c, eps, niter, fixedConstraints=[P0,P1])
        i = G.map(i, d)
        out.append(i)
    if ntype == 1: out = C.convertArray2Hexa(out)
    out = T.join(out)
    return out

# map curvature
# h prop a alpha si angle vif
# h impose sur radius a 1/3 ou 2/3
def mapCurvature(a, N, factor=1., sharpAngle=30.):
    L = D.getLength(a)
    h = L/(N-1)
    ang = D.getSharpestAngle(a)
    radius = D.getCurvatureRadius(a)
    # split at local max of radius + angles
    alpha = ang[1][0]
    #rad = radius[1][0]
    npts = alpha.size
    out = []
    #radiusp = 0.
    split = [0]
    for i in range(npts):
        alpha0 = abs(alpha[i]-180.)
        #radius0 = rad[i]
        if alpha0 > 30.: split.append(i)
    split.append(npts-1)
    
    nsplit = len(split)
    for x in range(nsplit):
        isp = split[x]
        alp = abs(alpha[isp]-180.)/180.
        alp = max(1.-alp,0.1)
        # set in split
        setH(a, split[x], h*alp)
        # set middle
        if x < nsplit-1: setH(a, 0.5*(split[x]+split[x+1]), h) 
    out = enforceh(a, N)
    return out
