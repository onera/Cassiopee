# - tkMapEdge -
"""Remap mesh edges."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import CPlot.iconics as iconics
import Converter.Internal as Internal
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import Post.PyTree as P
import KCore.Vector as Vector

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# CopyDistrib pour des zones edges
#==============================================================================
def copyDistrib1D(source):
    fail = False
    nzs = CPlot.getSelectedZones()
    errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            zp = G.map(z, source, 1)
            cad = Internal.getNodeFromName1(z, 'CAD')
            if cad is not None: zp[2].append(cad)
            CTK.replace(CTK.t, nob, noz, zp)
        except Exception as e:
            fail = True; errors += [0,str(e)]
    if len(errors)>0: Panels.displayErrors(errors, header='Error: copyDistrib1D')
    return fail

#==============================================================================
# Stretch pour une zone edge
#==============================================================================
def stretch1D(h):
    fail = False
    nzs = CPlot.getSelectedZones()
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    dims = Internal.getZoneDim(z)
    try:
        if dims[0] == 'Unstructured': a = C.convertBAR2Struct(z)
        else: a = Internal.copyTree(z)
    except Exception as e:
        #print('Error: stretch1D: %s.'%str(e))
        Panels.displayErrors([0,str(e)], header='Error: stretch1D')
        return True # Fail

    ind = CPlot.getActivePointIndex()
    if ind == []: return True # Fail
    ind = ind[0]

    zp = D.getCurvilinearAbscissa(z)
    l = D.getLength(a)
    a = D.getCurvilinearAbscissa(a)
    distrib = C.cpVars(a, 's', a, 'CoordinateX')
    C._initVars(distrib, 'CoordinateY', 0.)
    C._initVars(distrib, 'CoordinateZ', 0.)
    C._rmVars(distrib, 's')

    N = dims[1]
    val = C.getValue(zp, 's', ind)

    Xc = CPlot.getActivePoint()
    valf = val
    Pind = C.getValue(z, 'GridCoordinates', ind)
    if ind < N-1: # cherche avec indp1
        Pindp1 = C.getValue(z, 'GridCoordinates', ind+1)
        v1 = Vector.sub(Pindp1, Pind)
        v2 = Vector.sub(Xc, Pind)
        if Vector.dot(v1,v2) >= 0:
            val2 = C.getValue(zp, 's', ind+1)
            alpha = Vector.norm(v2)/Vector.norm(v1)
            valf = val+alpha*(val2-val)
    if ind > 0 and val == valf: # cherche avec indm1
        Pindm1 = C.getValue(z, 'GridCoordinates', ind-1)
        v1 = Vector.sub(Pindm1, Pind)
        v2 = Vector.sub(Xc, Pind)
        if Vector.dot(v1,v2) >= 0:
            val2 = C.getValue(zp, 's', ind-1)
            alpha = Vector.norm(v2)/Vector.norm(v1)
            valf = val+alpha*(val2-val)

    if h < 0: # enforce point
        distrib = G.enforcePoint(distrib, valf)
    else: # enforce h
        if val == 0: distrib = G.enforcePlusX(distrib, h/l, N//10, 1)
        elif val == 1: distrib = G.enforceMoinsX(distrib, h/l, N//10, 1)
        else: distrib = G.enforceX(distrib, valf, h/l, N//10, 1)
    try:
        a1 = G.map(z, distrib)
        cad = Internal.getNodeFromName1(z, 'CAD')
        if cad is not None: a1[2].append(cad)
        CTK.replace(CTK.t, nob, noz, a1)
    except Exception as e:
        fail = True
        Panels.displayErrors([0,str(e)], header='Error: stretch1D')
    return fail

#==============================================================================
# Set enforce h as a sizemap
def setEnforce(event=None):
    nzs = CPlot.getSelectedZones()
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    h = CTK.varsFromWidget(VARS[1].get(), 1)
    setEnforceZ(z, h, VARS[8].get())

# Perform a setH on z with clicked point
def setEnforceZ(z, h, mode):
    if len(h) != 1:
        CTK.TXT.insert('START', 'Invalid spacing.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    h = h[0]

    npts = C.getNPts(z)
    ind = CPlot.getActivePointIndex()
    if ind == []: return True
    ind = ind[0]
    if mode == 'HFactor':
        P0 = C.getValue(z, 'GridCoordinates', ind)
        if ind == npts-1: P1 = C.getValue(z, 'GridCoordinates', ind-1)
        else: P1 = C.getValue(z, 'GridCoordinates', ind+1)
        hloc = Vector.norm(Vector.sub(P1,P0))
        h = h*hloc
        #print("setting %f"%hloc)
    D.setH(z, ind, h)
    CTK.TXT.insert('START', 'Spacing set to %f.\n'%h)

# Pass in clik mode for setting H
def setEnforceMode(event=None):
    import time
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    W = WIDGETS['enforceMode']
    if not CTK.__BUSY__:
        CPlot.unselectAllZones()
        CTK.__BUSY__ = True
        TTK.sunkButton(W)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            l = []
            while l == []:
                nz = CPlot.getSelectedZone()
                l = CPlot.getActivePointIndex()
                time.sleep(CPlot.__timeStep__)
                W.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                nob = CTK.Nb[nz]+1
                noz = CTK.Nz[nz]
                z = CTK.t[2][nob][2][noz]
                setEnforceZ(z, CTK.varsFromWidget(VARS[1].get(), 1), VARS[8].get())

        CTK.__BUSY__ = False
        TTK.raiseButton(W)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(W)
        CPlot.setState(cursor=0)

def enforceH(event=None):
    v = CTK.varsFromWidget(VARS[10].get(), 1)
    if len(v) != 1:
        CTK.TXT.insert('START', 'Invalid number of points or factor.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    v = v[0]
    CTK.saveTree()
    CPlot.setState(cursor=1)
    nzs = CPlot.getSelectedZones()
    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        zones.append(z)
    npts = C.getNPts(zones)
    if VARS[9].get() == 'NFactor': N = v*npts
    else: N = v
    D._enforceh(zones, N=N)
    for c, nz in enumerate(nzs):
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        CTK.replace(CTK.t, nob, noz, zones[c])
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.TXT.insert('START', 'Spacings enforced.\n')

    # CAD remesh if possible
    edges = getSelection(nzs)
    remeshCAD(edges)
    CPlot.setState(cursor=0)

#==============================================================================
# Smooth pour les zones edges
#==============================================================================
def smooth1D(niter, eps):
    fail = False
    nzs = CPlot.getSelectedZones()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            D._smooth(z, eps, niter)
            CTK.replace(CTK.t, nob, noz, z)
        except Exception as e:
            fail = True
            Panels.displayErrors([0,str(e)], header='Error: smooth1D')
    return fail

#==============================================================================
# Refine pour les zones edges
#==============================================================================
def refine1D(density, npts, factor):
    fail = False
    nzs = CPlot.getSelectedZones()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dims = Internal.getZoneDim(z)
        try:
            if dims[0] == 'Unstructured':
                a = C.convertBAR2Struct(z); np = dims[1]
            else: a = z; np = dims[1]*dims[2]*dims[3]
            if factor < 0: factor = (npts-1.)/(np-1)
            G._refine(a, factor, 1)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True
            Panels.displayErrors([0,str(e)], header='Error: refine1D')
    return fail

#==============================================================================
# Uniformize pour les zones edges
#==============================================================================
def uniformize1D(density, h, npts, gnpts, factor):
    fail = False
    nzs = CPlot.getSelectedZones()
    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        zones.append(z)
    try:
        if npts > 0:
            for z in zones: D._uniformize(z, npts, h, factor, density)
        else:
            D._uniformize(zones, gnpts, h, factor, density)
        for c, nz in enumerate(nzs):
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            CTK.replace(CTK.t, nob, noz, zones[c])
    except Exception as e:
        fail = True
        Panels.displayErrors([0,str(e)], header='Error: uniformize1D')
    return fail

#==============================================================================
# Applique une fonction pour une zone structuree 2D
# ntype=0: uniformize
# ntype=1: refine
# ntype=2: stretch
# ntype=3: copyDistrib
# ntype=4: smooth
# Retourne True (failed), False (success)
#==============================================================================
def apply2D(density, h, npts, factor, ntype=0):
    nzs = CPlot.getSelectedZones()
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    ret = getEdges2D(zone, 0.)
    if ret is None: return True
    (m, r, u, ind) = ret
    out = []
    # Applique la fonction sur m[0] (edge a modifier)
    i = m[0]
    dims = Internal.getZoneDim(i)
    np = dims[1]*dims[2]*dims[3]
    if ntype == 0: # uniformize
        if density > 0: npts = D.getLength(i)*density
        if factor > 0: npts = np*factor
        if h > 0: npts = D.getLength(i)/h+1
        npts = int(max(npts, 2))
        distrib = G.cart((0,0,0), (1./(npts-1.),1,1), (npts,1,1))
        b = G.map(i, distrib)
    elif ntype == 1: # refine
        if factor < 0: factor = (npts-1.)/(np-1)
        else: npts = factor*(np-1)+1
        b = G.refine(i, factor, 1)
    elif ntype == 2: # stretch (factor=h)
        h = factor
        l = D.getLength(i)
        a = D.getCurvilinearAbscissa(i)
        distrib = C.cpVars(a, 's', a, 'CoordinateX')
        C._initVars(distrib, 'CoordinateY', 0.)
        C._initVars(distrib, 'CoordinateZ', 0.)
        distrib = C.rmVars(distrib, 's')
        N = dims[1]
        val = C.getValue(a, 's', ind)
        Xc = CPlot.getActivePoint()
        valf = val
        Pind = C.getValue(i, 'GridCoordinates', ind)
        if ind < N-1: # cherche avec indp1
            Pindp1 = C.getValue(i, 'GridCoordinates', ind+1)
            v1 = Vector.sub(Pindp1, Pind)
            v2 = Vector.sub(Xc, Pind)
            if Vector.dot(v1,v2) >= 0:
                val2 = C.getValue(a, 's', ind+1)
                alpha = Vector.norm(v2)/Vector.norm(v1)
                valf = val+alpha*(val2-val)
        if ind > 0 and val == valf: # cherche avec indm1
            Pindm1 = C.getValue(i, 'GridCoordinates', ind-1)
            v1 = Vector.sub(Pindm1, Pind)
            v2 = Vector.sub(Xc, Pind)
            if Vector.dot(v1,v2) >= 0:
                val2 = C.getValue(a, 's', ind-1)
                alpha = Vector.norm(v2)/Vector.norm(v1)
                valf = val+alpha*(val2-val)
        if h < 0: distrib = G.enforcePoint(distrib, valf)
        else:
            if val == 0: distrib = G.enforcePlusX(distrib, h/l, N/10, 1)
            elif val == 1: distrib = G.enforceMoinsX(distrib, h/l, N/10, 1)
            else: distrib = G.enforceX(distrib, valf, h/l, N/10, 1)
        b = G.map(i, distrib)
    elif ntype == 3: # copyDistrib (factor=source=edge pour l'instant)
        source = factor
        b = G.map(i, source, 1)
    elif ntype == 4: # smooth (factor=eps, npts=niter)
        niter = npts
        eps = factor
        a = D.getCurvilinearAbscissa(i)
        distrib = C.cpVars(a, 's', a, 'CoordinateX')
        C._initVars(distrib, 'CoordinateY', 0.)
        C._initVars(distrib, 'CoordinateZ', 0.)
        distrib = C.rmVars(distrib, 's')
        bornes = P.exteriorFaces(distrib)
        distrib = T.smooth(distrib, eps=eps, niter=niter,
                           fixedConstraints=[bornes])
        b = G.map(i, distrib, 1)
    dimb = Internal.getZoneDim(b)
    npts = dimb[1]
    out.append(b)
    # Raffine les edges si necessaires
    if npts != np:
        ret = getEdges2D(zone, 2.)
        if ret is None: return True
        (m, r, u, ind) = ret
    for i in r:
        dims = Internal.getZoneDim(i)
        np = dims[1]*dims[2]*dims[3]
        factor = (npts-1.)/(np-1)  # npts de m
        b = G.refine(i, factor, 1)
        out.append(b)
    # Garde les autres
    out += u
    #tp = C.newPyTree(['Base', out])
    #C.convertPyTree2File(tp, 'edges.cgns')

    # Rebuild
    try:
        b = G.TFI(out)
        # Projection du patch interieur
        #dimsb = Internal.getZoneDim(b)
        #bs = T.subzone(b, (2,2,1), (dimsb[1]-1,dimsb[2]-1,1))
        #bs = T.projectOrtho(bs, [zone])
        #b = T.patch(b, bs, position=(2,2,1))

        #tp = C.newPyTree(['Base'])
        #tp[2][1][2] += [b, zone]
        #C.convertPyTree2File(tp, 'face.cgns')
        b = T.projectOrtho(b, [zone])
        CTK.replace(CTK.t, nob, noz, b)
        return False
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: apply2D')
        return True

#==============================================================================
# Applique une fonction pour une zone structuree 3D
# ntype=0: uniformize
# ntype=1: refine
# ntype=2: stretch
# ntype=3: copyDistrib
# ntype=4: smooth
# Retourne True (failed), False (success)
#==============================================================================
def apply3D(density, h, npts, factor, ntype):
    nzs = CPlot.getSelectedZones()
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    ret = getEdges3D(zone, 0.)
    if ret is None: return True
    (m, r, f, ue, uf, ind) = ret
    out = []
    # Applique la fonction sur m
    i = m[0]
    dims = Internal.getZoneDim(i)
    np = dims[1]*dims[2]*dims[3]
    if ntype == 0: # uniformize
        if density > 0: npts = D.getLength(i)*density
        if factor > 0: npts = np*factor
        if h > 0: npts = D.getLength(i)/h+1
        npts = int(max(npts, 2))
        distrib = G.cart((0,0,0), (1./(npts-1.),1,1), (npts,1,1))
        b = G.map(i, distrib)
    elif ntype == 1: # refine
        if factor < 0: factor = (npts-1.)/(np-1)
        else: npts = factor*(np-1)+1
        b = G.refine(i, factor, 1)
    elif ntype == 2: # stretch (factor=h)
        h = factor
        l = D.getLength(i)
        a = D.getCurvilinearAbscissa(i)
        distrib = C.cpVars(a, 's', a, 'CoordinateX')
        C._initVars(distrib, 'CoordinateY', 0.)
        C._initVars(distrib, 'CoordinateZ', 0.)
        distrib = C.rmVars(distrib, 's')
        N = dims[1]
        val = C.getValue(a, 's', ind)
        Xc = CPlot.getActivePoint()
        valf = val
        Pind = C.getValue(i, 'GridCoordinates', ind)
        if ind < N-1: # cherche avec indp1
            Pindp1 = C.getValue(i, 'GridCoordinates', ind+1)
            v1 = Vector.sub(Pindp1, Pind)
            v2 = Vector.sub(Xc, Pind)
            if Vector.dot(v1,v2) >= 0:
                val2 = C.getValue(a, 's', ind+1)
                alpha = Vector.norm(v2)/Vector.norm(v1)
                valf = val+alpha*(val2-val)
        if ind > 0 and val == valf: # cherche avec indm1
            Pindm1 = C.getValue(i, 'GridCoordinates', ind-1)
            v1 = Vector.sub(Pindm1, Pind)
            v2 = Vector.sub(Xc, Pind)
            if Vector.dot(v1,v2) >= 0:
                val2 = C.getValue(a, 's', ind-1)
                alpha = Vector.norm(v2)/Vector.norm(v1)
                valf = val+alpha*(val2-val)
        if h < 0: distrib = G.enforcePoint(distrib, valf)
        else:
            if val == 0: distrib = G.enforcePlusX(distrib, h/l, N/10, 1)
            elif val == 1: distrib = G.enforceMoinsX(distrib, h/l, N/10, 1)
            else: distrib = G.enforceX(distrib, valf, h/l, N/10, 1)
        b = G.map(i, distrib)
    elif ntype == 3:
        source = factor
        b = G.map(i, source, 1)
    elif ntype == 4: # smooth (factor=eps, npts=niter)
        niter = npts
        eps = factor
        a = D.getCurvilinearAbscissa(i)
        distrib = C.cpVars(a, 's', a, 'CoordinateX')
        C._initVars(distrib, 'CoordinateY', 0.)
        C._initVars(distrib, 'CoordinateZ', 0.)
        distrib = C.rmVars(distrib, 's')
        bornes = P.exteriorFaces(distrib)
        distrib = T.smooth(distrib, eps=eps, niter=niter,
                           fixedConstraints=[bornes])
        b = G.map(i, distrib, 1)
    dimb = Internal.getZoneDim(b)
    npts = dimb[1]
    out.append(b)
    # Raffine les edges si necessaires
    if npts != np:
        ret = getEdges3D(zone, 2.)
        if ret is None: return True
        (m, r, f, ue, uf, ind) = ret
    for i in r:
        dims = Internal.getZoneDim(i)
        np = dims[1]*dims[2]*dims[3]
        factor = (npts-1.)/(np-1)  # npts de m
        b = G.refine(i, factor, 1)
        out.append(b)
    # Garde les autres
    out += ue
    outf = []
    # Rebuild les faces
    for i in f:
        # trouve les edges de la face
        edges = P.exteriorFacesStructured(i)
        match = []
        for e in edges:
            dime = Internal.getZoneDim(e)
            np = dime[1]-1
            P0 = C.getValue(e, Internal.__GridCoordinates__, 0)
            P1 = C.getValue(e, Internal.__GridCoordinates__, np)
            for ei in out: # retrouve les edges par leurs extremites
                dimei = Internal.getZoneDim(ei)
                npi = dimei[1]-1
                Q0 = C.getValue(ei, Internal.__GridCoordinates__, 0)
                Q1 = C.getValue(ei, Internal.__GridCoordinates__, npi)
                t1 = Vector.norm2(Vector.sub(P0,Q0))
                t2 = Vector.norm2(Vector.sub(P1,Q1))
                if t1 < 1.e-12 and t2 < 1.e-12: match.append(ei)
        if len(match) == 4: # OK
            fn = G.TFI(match)
            # Projection du patch interieur
            #dimsf = Internal.getZoneDim(fn)
            #fns = T.subzone(fn, (2,2,1), (dimsf[1]-1,dimsf[2]-1,1))
            #fns = T.projectOrtho(fns, [i])
            #fn = T.patch(fn, fns, position=(2,2,1))
            #fn = T.projectOrtho(fn, [i])
            outf.append(fn)
        else: return True
    outf += uf
    try:
        b = G.TFI(outf)
        CTK.replace(CTK.t, nob, noz, b)
        return False
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: apply3D')
        return True

#==============================================================================
def forceIJK(i,j,k,ni,nj,nk):
    bc = 0
    if i == 1 or i == ni: bc += 1
    if j == 1 or j == nj: bc += 1
    if k == 1 or k == nk: bc += 1
    if bc < 2: # force bc
        di1 = i-1; di2 = ni-i
        dj1 = j-1; dj2 = nj-j
        dk1 = k-1; dk2 = nk-k
        se = [(di1,'i1'),(di2,'i2'),(dj1,'j1'),(dj2,'j2'),(dk1,'dk1'),(dk2,'k2')]
        se.sort()
        for s in se:
            if s[0] != 0: bcf = s[1]; break

        if bcf == 'i1': i = 1
        elif bcf == 'i2': i = ni
        elif bcf == 'j1': j = 1
        elif bcf == 'j2': j = nj
        elif bcf == 'k1': k = 1
        elif bcf == 'k2': k = nk
    return (i,j,k)

#==============================================================================
# Recupere les edges concernes par le remaillage
# Retourne None si failed
# Sinon retourne les edges a modifier, a raffiner ou inchanges
# Et l'indice du point clicke sur l'edge a modifier
#==============================================================================
def getEdges2D(zone, factor):
    dim = Internal.getZoneDim(zone)
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    l = CPlot.getActivePointIndex()
    if l == []: return None
    i = l[2]; j = l[3]; k = l[4]
    e1 = T.subzone(zone, (1,1,1), (ni,1,1))
    e2 = T.subzone(zone, (1,nj,1), (ni,nj,1))
    e3 = T.subzone(zone, (1,1,1), (1,nj,1))
    e4 = T.subzone(zone, (ni,1,1), (ni,nj,1))
    u = [e1,e2,e3,e4]
    (i,j,k) = forceIJK(i,j,k,ni,nj,nk)
    r = []
    if i == 1:
        ind = j-1
        m = [e3]; u.remove(e3)
        if factor != 1.: r = [e4]; u.remove(e4)
    elif i == ni:
        ind = j-1
        m = [e4]; u.remove(e4)
        if factor != 1.: r = [e3]; u.remove(e3)
    elif j == 1:
        ind = i-1
        m = [e1]; u.remove(e1)
        if factor != 1.: r = [e2]; u.remove(e2)
    elif j == nj:
        ind = i-1
        m = [e2]; u.remove(e2)
        if factor != 1.: r = [e1]; u.remove(e1)
    return m, r, u, ind

#==============================================================================
# Recupere les edges/faces concernes par le remaillage
# Retourne None si failed
# Sinon retourne les edges a modifier, a raffiner ou inchanges
# et les faces a modifier ou inchangees
# Et l'indice du point clicke sur l'edge a modifier
#==============================================================================
def getEdges3D(zone, factor):
    dim = Internal.getZoneDim(zone)
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    l = CPlot.getActivePointIndex()
    if l == []: return None
    i = l[2]; j = l[3]; k = l[4]
    e1 = T.subzone(zone, (1,1,1), (ni,1,1))
    e2 = T.subzone(zone, (1,nj,1), (ni,nj,1))
    e3 = T.subzone(zone, (1,1,1), (1,nj,1))
    e4 = T.subzone(zone, (ni,1,1), (ni,nj,1))

    e5 = T.subzone(zone, (1,1,nk), (ni,1,nk))
    e6 = T.subzone(zone, (1,nj,nk), (ni,nj,nk))
    e7 = T.subzone(zone, (1,1,nk), (1,nj,nk))
    e8 = T.subzone(zone, (ni,1,nk), (ni,nj,nk))

    e9 = T.subzone(zone, (1,1,1), (1,1,nk))
    e10 = T.subzone(zone, (ni,1,1), (ni,1,nk))
    e11 = T.subzone(zone, (1,nj,1), (1,nj,nk))
    e12 = T.subzone(zone, (ni,nj,1), (ni,nj,nk))

    f1 = T.subzone(zone, (1,1,1), (ni,nj,1))
    f2 = T.subzone(zone, (1,1,nk), (ni,nj,nk))
    f3 = T.subzone(zone, (1,1,1), (ni,1,nk))
    f4 = T.subzone(zone, (1,nj,1), (ni,nj,nk))
    f5 = T.subzone(zone, (1,1,1), (1,nj,nk))
    f6 = T.subzone(zone, (ni,1,1), (ni,nj,nk))

    ue = [e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12] # unmodified edges
    uf = [f1,f2,f3,f4,f5,f6] # unmodified faces

    #temp = C.newPyTree(['Base'])
    #temp[2][1][2] += ue
    #temp[2][1][2] += uf
    #C.convertPyTree2File(temp, 'edges.cgns')

    (i,j,k) = forceIJK(i,j,k,ni,nj,nk)

    r = [] # refined edges
    if i == 1 and j == 1:
        ind = k-1
        m = [e9]
        f = [f3,f5]
        if factor != 1.: r = [e11,e12,e10]; f += [f4,f6]

    elif i == ni and j == 1:
        ind = k-1
        m = [e10]
        f = [f6,f3]
        if factor != 1.:
            r = [e9,e11,e12]
            f += [f5,f4]

    elif i == 1 and j == nj:
        ind = k-1
        m = [e11]
        f = [f5,f4]
        if factor != 1.: r = [e9,e10,e12]; f += [f3,f6]

    elif i == ni and j == nj:
        ind = k-1
        m = [e12]
        f = [f6,f4]
        if factor != 1.: r = [e10,e9,e11]; f += [f3,f5]

    elif i == 1 and k == 1:
        ind = j-1
        m = [e3]
        f = [f1,f5]
        if factor != 1.: r = [e4,e8,e7]; f += [f2,f6]

    elif i == ni and k == 1:
        ind = j-1
        m = [e4]
        f = [f6,f1]
        if factor != 1.: r = [e3,e8,e7]; f += [f2,f5]

    elif i == 1 and k == nk:
        ind = j-1
        m = [e7]
        f = [f5,f2]
        if factor != 1.: r = [e8,e4,e3]; f += [f6,f1]

    elif i == ni and k == nk:
        ind = j-1
        m = [e8]
        f = [f6,f2]
        if factor != 1.: r = [e7,e3,e4]; f += [f5,f1]

    elif j == 1 and k == 1:
        ind = i-1
        m = [e1]
        f = [f1,f3]
        if factor != 1.: r = [e2,e6,e5]; f += [f2,f4]

    elif j == nj and k == 1:
        ind = i-1
        m = [e2]
        f = [f4,f1]
        if factor != 1.: r = [e1,e6,e5]; f += [f2,f3]

    elif j == 1 and k == nk:
        ind = i-1
        m = [e5]
        f = [f2,f3]
        if factor != 1.: r = [e6,e2,e1]; f += [f1,f4]

    elif j == nj and k == nk:
        ind = i-1
        m = [e6]
        f = [f2,f4]
        if factor != 1.: r = [e5,e1,e2]; f += [f3,f1]

    for i in m+r: ue.remove(i)
    for i in f: uf.remove(i)
    return m, r, f, ue, uf, ind

#==============================================================================
def uniformize(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    rtype = VARS[2].get()
    density = -1; npts = -1; gnpts = -1; factor = -1; h = -1
    if rtype == 'Density':
        density = CTK.varsFromWidget(VARS[0].get(), 1)
        if len(density) != 1:
            CTK.TXT.insert('START', 'Invalid points density.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        density = density[0]
    elif rtype == 'Npts':
        npts = CTK.varsFromWidget(VARS[0].get(), 2)
        if len(npts) != 1:
            CTK.TXT.insert('START', 'Invalid number of points.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        npts = npts[0]
    elif rtype == 'GNpts': # global number of points
        gnpts = CTK.varsFromWidget(VARS[0].get(), 2)
        if len(gnpts) != 1:
            CTK.TXT.insert('START', 'Invalid number of points.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        gnpts = gnpts[0]
    elif rtype == 'NFactor':
        factor = CTK.varsFromWidget(VARS[0].get(), 1)
        if len(factor) != 1:
            CTK.TXT.insert('START', 'Invalid number factor.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        factor = factor[0]
    elif rtype == 'H':
        h = CTK.varsFromWidget(VARS[0].get(), 1)
        if len(h) != 1:
            CTK.TXT.insert('START', 'Invalid h step.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        h = h[0]

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['frame'])

    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1:
            fail = apply3D(density, h, npts, factor, ntype=0)
        elif dim[2] != 1 and dim[3] == 1:
            fail = apply2D(density, h, npts, factor, ntype=0)
        else: fail = uniformize1D(density, h, npts, gnpts, factor)
    else: fail = uniformize1D(density, h, npts, gnpts, factor) # all zones

    if not fail:
        CTK.TXT.insert('START', 'Uniformize successfull.\n')
    else:
        CTK.TXT.insert('START', 'Uniformize edge fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

    # CAD remesh if possible
    edges = getSelection(nzs)
    remeshCAD(edges)
    CTK.setCursor(0, WIDGETS['frame'])

#==============================================================================
def enforce(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    h = CTK.varsFromWidget(VARS[1].get(), 1)
    if len(h) != 1:
        CTK.TXT.insert('START', 'Invalid spacing.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    h = h[0]

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['frame'])

    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1:
            fail = apply3D(1., -1, 1, h, ntype=2)
        elif dim[2] != 1 and dim[3] == 1:
            fail = apply2D(1., -1, 1, h, ntype=2)
        else: fail = stretch1D(h)
    else: fail = stretch1D(h)

    if not fail:
        CTK.TXT.insert('START', 'Stretch successfull.\n')
    else:
        CTK.TXT.insert('START', 'stretch failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    #C.fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['frame'])

#==============================================================================
def refine(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    rtype = VARS[4].get()
    factor = -1; npts = 2
    if rtype == 'NFactor':
        factor = CTK.varsFromWidget(VARS[5].get(), 1)
        if len(factor) != 1:
            CTK.TXT.insert('START', 'Invalid refinement factor.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
        factor = factor[0]
    else:
        npts = CTK.varsFromWidget(VARS[5].get(), 2)
        if len(npts) != 1:
            CTK.TXT.insert('START', 'Invalid number of points.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
        npts = npts[0]

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['frame'])

    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1:
            fail = apply3D(1., -1, npts, factor, ntype=1)
        elif dim[2] != 1 and dim[3] == 1:
            fail = apply2D(1., -1, npts, factor, ntype=1)
        else: fail = refine1D(1., npts, factor) # all zones
    else: fail = refine1D(1., npts, factor) # all zones

    if not fail:
        CTK.TXT.insert('START', 'Refine successfull.\n')
    else:
        CTK.TXT.insert('START', 'Refine fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

    # CAD remesh if possible
    edges = getSelection(nzs)
    remeshCAD(edges)
    CTK.setCursor(0, WIDGETS['frame'])

#==============================================================================
def smooth(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    eps = CTK.varsFromWidget(VARS[7].get(), 1)
    if len(eps) != 1:
        CTK.TXT.insert('START', 'Invalid smoother power.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    eps = eps[0]
    niter = CTK.varsFromWidget(VARS[6].get(), 2)
    if len(niter) != 1:
        CTK.TXT.insert('START', 'Invalid number of smoother iterations.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    niter = niter[0]

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['frame'])

    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1:
            fail = apply3D(1., -1, niter, eps, ntype=4)
        elif dim[2] != 1 and dim[3] == 1:
            fail = apply2D(1., -1, niter, eps, ntype=4)
        else: fail = smooth1D(niter, eps)
    else: fail = smooth1D(niter, eps) # all zones

    if not fail:
        CTK.TXT.insert('START', 'Smooth successfull.\n')
    else:
        CTK.TXT.insert('START', 'Smooth fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

    # CAD remesh if possible
    edges = getSelection(nzs)
    remeshCAD(edges)
    CTK.setCursor(0, WIDGETS['frame'])

#==============================================================================
def setSourceEdge():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    selected = ''
    if len(nzs) > 1:
        CTK.TXT.insert('START', 'Only first zone is taken.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    selected = CTK.t[2][nob][0]+'/'+z[0]
    # ajoute l'indice du point actif
    pt = CPlot.getActivePointIndex()
    selected += ';'+str(pt[2:])
    VARS[3].set(selected)

#==============================================================================
def copyDistrib():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['frame'])

    # source edge
    v = VARS[3].get()
    v = v.split(';')
    try: pt = eval(v[1])
    except: pt = None
    v = v[0]
    v = v.lstrip(); v = v.rstrip()
    sname = v.split('/', 1)
    edge = []
    bases = Internal.getNodesFromName1(CTK.t, sname[0])
    if bases != []:
        nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
        for z in nodes:
            if z[0] == sname[1]: edge.append(z)
    if edge == []:
        CTK.TXT.insert('START', 'Invalid source edge.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    source = edge[0]
    dimSource = Internal.getZoneDim(source)
    if dimSource[0] == 'Unstructured':
        CTK.TXT.insert('START', 'Invalid source edge.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)

    # subzone de source eventuellement
    if dimSource[4] == 2: # source est 2D
        ind = CPlot.getActivePointIndex()
        ni = dimSource[1]; nj = dimSource[2]
        deltai = min(ind[2]-ni,ind[2]-1)
        deltaj = min(ind[3]-nj,ind[3]-1)
        if deltai < deltaj:
            source = T.subzone(source, (1,pt[1],pt[2]),(ni,pt[1],pt[2]))
        else:
            source = T.subzone(source, (pt[0],1,pt[2]),(pt[0],nj,pt[2]))
    elif dimSource[4] == 3: # source est 3D
        ind = CPlot.getActivePointIndex()
        ni = dimSource[1]; nj = dimSource[2]; nk = dimSource[3]
        deltai = min(ind[2]-ni,ind[2]-1)
        deltaj = min(ind[3]-nj,ind[3]-1)
        deltak = min(ind[4]-nk,ind[4]-1)
        if deltai < deltaj and deltai < deltak:
            source = T.subzone(source, (1,pt[1],pt[2]),(ni,pt[1],pt[2]))
        elif deltaj < deltai and deltaj < deltak:
            source = T.subzone(source, (pt[0],1,pt[2]),(pt[0],nj,pt[2]))
        else:
            source = T.subzone(source, (pt[0],pt[1],1),(pt[0],pt[1],nk))
    else: # source is 1D
        source = Internal.copyTree(source)

    # Extrait la distribution en i
    source = D.getCurvilinearAbscissa(source)
    C._initVars(source, '{CoordinateX}={s}')
    C._initVars(source, 'CoordinateY', 0)
    C._initVars(source, 'CoordinateZ', 0)
    source = C.rmVars(source, 's')

    # Traitement
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1:
            fail = apply3D(1., -1, 1, source, ntype=3)
        elif dim[2] != 1 and dim[3] == 1:
            fail = apply2D(1., -1, 1, source, ntype=3)
        else: fail = copyDistrib1D(source)
    else: fail = copyDistrib1D(source) # all zones

    if not fail:
        CTK.TXT.insert('START', 'Distribution copy done.\n')
    else:
        CTK.TXT.insert('START', 'Distribution copy fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

    # CAD remesh if possible
    edges = getSelection(nzs)
    remeshCAD(edges)
    CTK.setCursor(0, WIDGETS['frame'])


#==============================================================================
# get selection from CTK.t
def getSelection(nzs):
    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        zones.append(CTK.t[2][nob][2][noz])
    return zones

#==============================================================================
# remesh CAD when an edge is modified
#==============================================================================
def remeshCAD(edges):
    try: import OCC.PyTree as OCC
    except: return
    if CTK.CADHOOK is None: return
    valids = []
    for e in edges:
        zdim = Internal.getZoneDim(e)
        if zdim[0] != 'Structured': continue
        if zdim[4] != 1: continue
        CAD = Internal.getNodeFromName1(e, 'CAD')
        if CAD is None: continue
        D._getCurvilinearAbscissa(e)
        valids.append(e)
    OCC._remeshTreeFromEdges(CTK.CADHOOK, CTK.t, valids)
    CTK.display(CTK.t)

#==============================================================================
# enforce h in edge locally
#==============================================================================
def enforceLocal(event=None):
    v = CTK.varsFromWidget(VARS[12].get(), 1)
    if len(v) != 1:
        CTK.TXT.insert('START', 'Invalid h or hfactor.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    v = v[0]
    width = WIDGETS['widthScale'].get() / 100.
    width = max(width, 0.1)
    mode = VARS[11].get()

    # get edge
    nzs = CPlot.getSelectedZones()
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(z)
    if dim[4] != 1:
        CTK.TXT.insert('START', 'Zone must be an edge.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    npts = C.getNPts(z)

    CTK.saveTree()
    CPlot.setState(cursor=1)

    # split zone of width
    ind = CPlot.getActivePointIndex()
    if ind == []: return
    ind = ind[0]

    P0 = C.getValue(z, 'GridCoordinates', ind)
    if ind == npts-1: P1 = C.getValue(z, 'GridCoordinates', ind-1)
    else: P1 = C.getValue(z, 'GridCoordinates', ind+1)
    hloc = Vector.norm(Vector.sub(P1,P0))
    if mode == 'HFactor':
        h = v*hloc; factor = v
    else: h = v; factor = h / hloc

    delta = int(width*npts)+1
    imin = ind+1-delta
    imax = ind+1+delta

    CAD = Internal.getNodeFromName1(z, 'CAD')
    render = Internal.getNodeFromName1(z, '.RenderInfo')

    #print("imin=",imin,"imax",imax,"ind",ind,"npts",npts)
    if imin > 1: z0 = T.subzone(z, (1,1,1), (imin,-1,-1))
    else: z0 = None
    if ind > 1: zp1 = T.subzone(z, (max(imin,1),1,1), (ind+1,-1,-1))
    else: zp1 = None
    if ind < npts-1: zp2 = T.subzone(z, (ind+1,1,1), (min(imax, npts),-1,-1))
    else: zp2 = None
    if imax < npts: z1 = T.subzone(z, (imax,1,1), (npts,-1,-1))
    else: z1 = None

    if zp1 is None and zp2 is None:
        if z0 is not None: zp1 = z0; z0 = None
        elif z1 is not None: zp2 = z1; z1 = None
        else: print("Error: can not remesh edge.")

    if zp1 is not None:
        P0 = C.getValue(zp1, 'GridCoordinates', 0)
        P1 = C.getValue(zp1, 'GridCoordinates', 1)
        h1 = Vector.norm(Vector.sub(P1,P0))
        L = D.getLength(zp1)
        if h+h1 > L: h1 = h

    if zp2 is not None:
        P0 = C.getValue(zp2, 'GridCoordinates', -1)
        P1 = C.getValue(zp2, 'GridCoordinates', -2)
        h2 = Vector.norm(Vector.sub(P1,P0))
        L = D.getLength(zp2)
        if h+h2 > L: h2 = h

    if zp1 is None: h1 = h2
    if zp2 is None: h2 = h1

    #if z0 is None: h1 = (h1+h)*0.5
    #if z1 is None: h2 = (h2+h)*0.5

    # D._setH(zp, ind-imin+1, h)
    # # guess a cool number of points
    # if factor == 1.:
    #     N = npts
    # elif factor < 1.:
    #     N = npts + (1./factor)/100.*npts
    #     N = int(N)+1
    # else:
    #     N = npts - (1./factor)/100.*npts
    #     N = int(N)+1
    # print("nbre de points=", N)
    # D._enforceh(zp, N=N)

    if zp1 is not None:
        d2 = D.distrib2(zp1, h1, h, algo=1)
        zp1 = G.map(zp1, d2)
    if zp2 is not None:
        d2 = D.distrib2(zp2, h, h2, algo=1)
        zp2 = G.map(zp2, d2)

    zo = None
    if zp1 is not None: zo = zp1
    if zp2 is not None:
        if zo is not None: zo = T.join(zo, zp2)
        else: zo = zp2
    if z0 is not None: zo = T.join(z0, zo)
    if z1 is not None: zo = T.join(zo, z1)
    zo[0] = z[0] # keep orig name and CAD
    zo[2].append(CAD)
    if render: zo[2].append(render)

    CTK.replace(CTK.t, nob, noz, zo)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.TXT.insert('START', 'Local spacing enforced.\n')

    # CAD remesh if possible
    edges = getSelection(nzs)
    remeshCAD(edges)
    CPlot.setState(cursor=0)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkMapEdge  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Map distributions on edges.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=2)
    Frame.columnconfigure(3, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkMapEdge')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Point density or Npts -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -1- Enforce step
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkMapEdgeEnforceHeight' in CTK.PREFS:
        V.set(CTK.PREFS['tkMapEdgeEnforceHeight'])
    # -2- Option for uniformize
    V = TK.StringVar(win); V.set('NFactor'); VARS.append(V)
    # -3- Source mesh for copy -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -4- Option for refine
    V = TK.StringVar(win); V.set('NFactor'); VARS.append(V)
    # -5-  Number of points/factor for refine
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -6- Smoothing iterations
    V = TK.StringVar(win); V.set('5'); VARS.append(V)
    if 'tkMapEdgeSmoothIt' in CTK.PREFS:
        V.set(CTK.PREFS['tkMapEdgeSmoothIt'])
    # -7- Smoothing eps
    V = TK.StringVar(win); V.set('0.5'); VARS.append(V)
    if 'tkMapEdgeSmoothEps' in CTK.PREFS:
        V.set(CTK.PREFS['tkMapEdgeSmoothEps'])
    # -8- Type of set enforce
    V = TK.StringVar(win); V.set('HFactor'); VARS.append(V)
    # -9- Type of h enforce
    V = TK.StringVar(win); V.set('NFactor'); VARS.append(V)
    # -10- Nbre de pts pour h enforce
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -11- Type of local enforce
    V = TK.StringVar(win); V.set('HFactor'); VARS.append(V)
    # -12- Factor for local enforce
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)

    # - Uniformize -
    B = TTK.Button(Frame, text="Uniformize", command=uniformize)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Uniformize an edge with regular spacing.')
    B = TTK.OptionMenu(Frame, VARS[2], 'NFactor', 'Density', 'Npts', 'GNpts', 'H')
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=7)
    B.grid(row=0, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Point multiplication factor/point density/number of points.')
    B.bind('<Return>', uniformize)

    # - Refine edge -
    B = TTK.Button(Frame, text="Refine", command=refine)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refine an edge keeping distribution.')
    B = TTK.OptionMenu(Frame, VARS[4], 'NFactor', 'Npts')
    B.grid(row=1, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White', width=7)
    B.grid(row=1, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refinement factor or number of points.')
    B.bind('<Return>', refine)

    # - Copy distribution -
    B = TTK.Button(Frame, command=setSourceEdge, text='Copy')
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get source edge defining distribution for copy.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=10)
    B.grid(row=2, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Source edge for distribution copy.')
    B = TTK.Button(Frame, text="Paste", command=copyDistrib)
    B.grid(row=2, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Paste distribution from source edge.')

    # - Smooth edge -
    B = TTK.Button(Frame, text="Smooth", command=smooth)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Smooth an edge distribution.')
    B = TTK.Entry(Frame, textvariable=VARS[7], background='White', width=7)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Smoother power.')
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White', width=7)
    B.grid(row=3, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of smoother iterations.')
    B.bind('<Return>', smooth)

    # - Enforce -
    #B = TTK.Button(Frame, text="Set", command=setEnforce)
    #B.grid(row=4, column=0, columnspan=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Set step on curve.')
    #B = TTK.OptionMenu(Frame, VARS[8], 'HFactor', 'H')
    #B.grid(row=4, column=1, sticky=TK.EW)
    #B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=7)
    #B.grid(row=4, column=2, columnspan=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Enforced spacing.')
    #B.bind('<Return>', setEnforce)
    #B = TTK.Button(Frame, command=setEnforceMode,
    #               image=iconics.PHOTO[8], padx=0, pady=0)
    #B.grid(row=4, column=3, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Set size mode.')
    #WIDGETS['enforceMode'] = B

    #B = TTK.Button(Frame, text="Enforce", command=enforceH)
    #B.grid(row=5, column=0, columnspan=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Enforce all given spacing.')
    #B = TTK.OptionMenu(Frame, VARS[9], 'NFactor', 'Npts')
    #B.grid(row=5, column=1, sticky=TK.EW)
    #B = TTK.Entry(Frame, textvariable=VARS[10], background='White', width=7)
    #B.grid(row=5, column=2, columnspan=2, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Enforced number of points.')
    #B.bind('<Return>', enforceH)

    # - Enforce local -
    B = TTK.Button(Frame, text="Enforce", command=enforceLocal)
    B.grid(row=6, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Enforce local spacing.')
    B = TTK.OptionMenu(Frame, VARS[11], 'HFactor', 'H')
    B.grid(row=6, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[12], background='White', width=7)
    B.grid(row=6, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Enforced variable.')
    B.bind('<Return>', enforceLocal)
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL,
                  showvalue=0, borderwidth=1, value=50)
    WIDGETS['widthScale'] = B
    B.grid(row=6, column=3, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['EdgeNoteBook'].add(WIDGETS['frame'], text='tkMapEdge')
    except: pass
    CTK.WIDGETS['EdgeNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['EdgeNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkMapEdgeEnforceHeight'] = VARS[1].get()
    CTK.PREFS['tkMapEdgeSmoothIt'] = VARS[6].get()
    CTK.PREFS['tkMapEdgeSmoothEps'] = VARS[7].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[1].set('1.e-6')
    VARS[6].set('5')
    VARS[7].set('0.5')
    CTK.PREFS['tkMapEdgeEnforceHeight'] = VARS[1].get()
    CTK.PREFS['tkMapEdgeSmoothIt'] = VARS[6].get()
    CTK.PREFS['tkMapEdgeSmoothEps'] = VARS[7].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkMapEdge '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
