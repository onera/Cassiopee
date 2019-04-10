# mapEdge
try:
    import Generator as G
    import Transform as T
    import Converter as C
    import Geom as D
    import KCore
    import numpy
    importOK = True
except: importOK = False

try: range = xrange
except: pass

def checkImport():
    if not importOK: 
        raise ImportError("mapEdge requires Converter, Generator, Transform.")
        return None

# uniformize a 1D mesh (regular step) - OK
def uniformize(a, N=100, h=-1., factor=-1, density=-1., sharpAngle=30.):
    """Unformize the distribution of points on a 1D-mesh."""
    if checkImport is None: return None
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(uniformize__(i, N, h, factor, density, sharpAngle))
        return b
    else:
        return uniformize__(a, N, h, factor, density, sharpAngle)

def uniformize__(a, N, h, factor, density, sharpAngle):
    if len(a) == 4: ntype = 1
    else: ntype = 0
    L = D.getLength(a)
    if density > 0: N = int(L*density)+1
    elif factor > 0: N = int(factor*(C.getNPts(a)-1))+1
    elif h > 0: N = int(L/(N-1))
    N = max(N, 2)
    # I can add splitConnexity if needed
    b = T.splitSharpEdges(a, alphaRef=sharpAngle) # b is BAR closed
    if ntype == 1: b = T.splitTBranches(b)
    out = []; Nt = 0
    nt = len(b); ct = 0; rest = 0
    for i in b:
        ct += 1
        Li = D.getLength(i)
        Ni = int(round(N*1.*(Li/L)))+1
        if Ni-1-N*1.*(Li/L) < 0:
           if rest == -1: rest = 0; Ni += 1
           elif rest == 0: rest = -1
           elif rest == +1: rest = 0
        elif Ni-1-N*1.*(Li/L) > 0: 
           if rest == +1: rest = 0; Ni -= 1
           elif rest == 0: rest = +1
           elif rest == -1: rest = 0
        Ni = max(Ni,2)
        Nt += Ni-1
        if ct == nt: Ni = Ni-Nt+N-1; Ni = max(Ni,2)
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

    out = []; Nt = 0
    nt = len(b); ct = 0; rest = 0
    for i in b:
        ct += 1
        i = C.convertBAR2Struct(i)
        Li = D.getLength(i)
        if factor < 0.:
            Ni = int(round(N*(C.getNPts(i)-1)/(NPts-1))+1)
            if Ni-1-N*1.*(Li/L) < 0:
             if rest == -1: rest = 0; Ni += 1
             elif rest == 0: rest = -1
             elif rest == +1: rest = 0
            elif Ni-1-N*1.*(Li/L) > 0:
             if rest == +1: rest = 0; Ni -= 1
             elif rest == 0: rest = +1
             elif rest == -1: rest = 0
            Ni = max(Ni,2)
            Nt += Ni-1
            if ct == nt: Ni = Ni-Nt+N-1
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

# set h step in place
def setH(a, ind, h):
    """Set the mesh step for a given index."""
    pos = KCore.isNamePresent(a, 'h')
    if pos == -1:
        C._initVars(a, 'h', 0.)
        pos = KCore.isNamePresent(a, 'h')
    if ind == -1: ind = a[1].shape[1]-1
    a[1][pos,ind] = h

# set factor in place
def setF(a, ind, f):
    """Set the mesh size factor for a given index."""
    pos = KCore.isNamePresent(a, 'f')
    if pos == -1:
        C._initVars(a, 'f', 0.)
        pos = KCore.isNamePresent(a, 'f')
    if ind == -1: ind = a[1].shape[1]-1
    a[1][pos,ind] = f

# Build a dstrib between 0. and 1. with h1 left, h2 right
def buildDistrib(h1, h2, N):
    Ni = int(round(1./h2))+1
    a = G.cart((0,0,0), (h2,1,1), (Ni,1,1))
    a[1][0,Ni-2] = 1.-h2
    b = G.enforcePlusX(a, h1, (Ni-2,N-Ni))
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

# h must be set at sharp angles
# Set: N + h field (step size)
# Or: h + f field (factor)
def enforceh(a, N=100, h=-1):
    """Enforce mesh size in a 1D-mesh."""
    if checkImport is None: return None
    if isinstance(a[0], list): 
        b = []
        for i in a:
            b.append(enforceh_(i, N, h))
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
            raise ValueError("h or f is not present.")
        else: 
            hl = numpy.copy(a[1][pos])
            factorMoy = moyenne(a, hl)
            if h > 0:
                href = h
                N = D.getLength(a)/(factorMoy*href)+1
                N = int(round(N))
            else:
                href = D.getLength(a)/((N-1.)*factorMoy)
            hl[:] = hl[:]*href
    else: hl = a[1][pos]
    npts = hl.size
    # Calcul h1s, h2s, i1s, i2s
    h1s=[]; h2s=[]; Ps=[]; i1s=[]; i2s=[]
    h1 = -1; h2 = -1; i1 = -1; i2 = -1
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
                #print 'h1',Li,h1,h2
                Pi = Li*1./(0.5*(h1+h2)+1.e-12)
                i1s.append(i1); i2s.append(i2)
                h1s.append(h1/Li); h2s.append(h2/Li)
                Ps.append(Pi)
                # loop
                i1 = i2; h1 = h2
    Pt = 0.
    for x in range(len(h1s)):
        Pi = Ps[x]
        Pt += Pi

    for x in range(len(h1s)):
        Pi = Ps[x]
        Pi = Pi/Pt*N
        #print Ps[x], Pi
        Ps[x] = int(round(Pi))+1

    out = []
    for x in range(len(h1s)):
        i1 = i1s[x]; i2 = i2s[x]
        h1 = h1s[x]; h2 = h2s[x]
        # subzone
        #print i1, i2, h1, h2
        sub = T.subzone(a, (i1+1,1,1), (i2+1,1,1))
        d = buildDistrib(h1, h2, Ps[x])            
        # remap
        c = G.map(sub, d)
        out.append(c)
    out = T.join(out)
    return out

# Enforce h at ind (STRUCT)
def enforce(a, h, ind, supp, add):
    c = D.getDistribution(a)
    L = D.getLength(a)
    if ind == 0:
        b = G.enforceMoinsX(b, h/L, supp, add)
    elif ind == a[2]-1:
        b = G.enforcePlusX(b, h/L, supp, add)
    else:
        b = G.enforceX(b, c[1][0,ind], h/L, supp, add)
    a = G.map(a, b)
    return a

# Super smooth - not working
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
    for i in b:
        i = C.convertBAR2Struct(i)
        c = D.getDistribution(i)
        c = T.smooth(c, eps, niter)
        i = G.map(i, c)
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
    rad = radius[1][0]
    npts = alpha.size
    out = []
    radiusp = 0.
    split = [0]
    for i in range(npts):
        alpha0 = abs(alpha[i]-180.)
        radius0 = rad[i]
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
