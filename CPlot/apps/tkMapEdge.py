# - edge mapping -
try: import Tkinter as TK
except: import tkinter as TK
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
        else: a = z
    except Exception as e:
        #print 'Error: stretch1D: %s.'%str(e)
        Panels.displayErrors([0,str(e)], header='Error: stretch1D')
        return True # Fail

    ind = CPlot.getActivePointIndex()
    if ind == []: return True # Fail
    ind = ind[0]
    
    l = D.getLength(a)
    a = D.getCurvilinearAbscissa(a)
    zp = D.getCurvilinearAbscissa(z)
    distrib = C.cpVars(a, 's', a, 'CoordinateX')
    C._initVars(distrib, 'CoordinateY', 0.)
    C._initVars(distrib, 'CoordinateZ', 0.)
    distrib = C.rmVars(distrib, 's')
    
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
        if val == 0: distrib = G.enforcePlusX(distrib, h/l, N/10, 1)
        elif val == 1: distrib = G.enforceMoinsX(distrib, h/l, N/10, 1)
        else: distrib = G.enforceX(distrib, valf, h/l, N/10, 1)
    try:
        a1 = G.map(a, distrib)
        CTK.replace(CTK.t, nob, noz, a1)
    except Exception as e:
        fail = True
        Panels.displayErrors([0,str(e)], header='Error: stretch1D')
    return fail

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
        dims = Internal.getZoneDim(z)
        try:
            if dims[0] == 'Unstructured': a = C.convertBAR2Struct(z)
            else: a = z
            a = D.getCurvilinearAbscissa(a)
            distrib = C.cpVars(a, 's', a, 'CoordinateX')
            C._initVars(distrib, 'CoordinateY', 0.)
            C._initVars(distrib, 'CoordinateZ', 0.)
            distrib = C.rmVars(distrib, 's')
            bornes = P.exteriorFaces(distrib)
            distrib = T.smooth(distrib, eps=eps, niter=niter, 
                               fixedConstraints=[bornes])
            b = G.map(a, distrib)
            CTK.replace(CTK.t, nob, noz, b)
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
            b = G.refine(a, factor, 1)
            CTK.replace(CTK.t, nob, noz, b)
        except Exception as e:
            fail = True
            Panels.displayErrors([0,str(e)], header='Error: refine1D')
    return fail

#==============================================================================
# Uniformize pour les zones edges
#==============================================================================
def uniformize1D(density, npts, factor):
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
            if density > 0: npts = D.getLength(a)*density
            if factor > 0: npts = np*factor[0]
            npts = int(max(npts, 2))
            distrib = G.cart((0,0,0), (1./(npts-1.),1,1), (npts,1,1))
            b = G.map(a, distrib)
            CTK.replace(CTK.t, nob, noz, b)
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
def apply2D(density, npts, factor, ntype=0):
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
        if factor > 0: npts = np*factor[0]
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
    #tp = C.newPyTree(['Base'])
    #tp[2][1][2] += out
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
def apply3D(density, npts, factor, ntype):
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
        if factor > 0: npts = np*factor[0]
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
                if (t1 < 1.e-12 and t2 < 1.e-12): match.append(ei)
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
    if (i == 1 and j == 1):
        ind = k-1
        m = [e9]
        f = [f3,f5]
        if (factor != 1.): r = [e11,e12,e10]; f += [f4,f6]

    elif (i == ni and j == 1):
        ind = k-1
        m = [e10]
        f = [f6,f3]
        if (factor != 1.): 
            r = [e9,e11,e12]
            f += [f5,f4]

    elif (i == 1 and j == nj):
        ind = k-1
        m = [e11]
        f = [f5,f4]
        if (factor != 1.): r = [e9,e10,e12]; f += [f3,f6]

    elif (i == ni and j == nj):
        ind = k-1
        m = [e12]
        f = [f6,f4]
        if (factor != 1.): r = [e10,e9,e11]; f += [f3,f5]

    elif (i == 1 and k == 1):
        ind = j-1
        m = [e3]
        f = [f1,f5]
        if (factor != 1.): r = [e4,e8,e7]; f += [f2,f6]

    elif (i == ni and k == 1):
        ind = j-1
        m = [e4]
        f = [f6,f1]
        if (factor != 1.): r = [e3,e8,e7]; f += [f2,f5]

    elif (i == 1 and k == nk):
        ind = j-1
        m = [e7]
        f = [f5,f2]
        if (factor != 1.): r = [e8,e4,e3]; f += [f6,f1]

    elif (i == ni and k == nk):
        ind = j-1
        m = [e8]
        f = [f6,f2]
        if (factor != 1.): r = [e7,e3,e4]; f += [f5,f1]

    elif (j == 1 and k == 1):
        ind = i-1
        m = [e1]
        f = [f1,f3]
        if (factor != 1.): r = [e2,e6,e5]; f += [f2,f4]

    elif (j == nj and k == 1):
        ind = i-1
        m = [e2]
        f = [f4,f1]
        if (factor != 1.): r = [e1,e6,e5]; f += [f2,f3]

    elif (j == 1 and k == nk):
        ind = i-1
        m = [e5]
        f = [f2,f3]
        if (factor != 1.): r = [e6,e2,e1]; f += [f1,f4]

    elif (j == nj and k == nk):
        ind = i-1
        m = [e6]
        f = [f2,f4]
        if (factor != 1.): r = [e5,e1,e2]; f += [f3,f1]

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

    type = VARS[2].get()
    density = -1; npts = 2; factor = -1
    if type == 'Density':
        density = CTK.varsFromWidget(VARS[0].get(), 1)
        if len(density) != 1:
            CTK.TXT.insert('START', 'Invalid points density.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
        density = density[0]
    elif type == 'Npts':
        npts = CTK.varsFromWidget(VARS[0].get(), 2)
        if len(npts) != 1:
            CTK.TXT.insert('START', 'Invalid number of points.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
        npts = npts[0]
    elif type == 'Factor':
        factor = CTK.varsFromWidget(VARS[0].get(), 1)
        if len(factor) != 1:
            CTK.TXT.insert('START', 'Invalid number factor.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')

    CTK.saveTree()

    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1: 
            fail = apply3D(density, npts, factor, ntype=0)
        elif dim[2] != 1 and dim[3] == 1: 
            fail = apply2D(density, npts, factor, ntype=0)
        else: fail = uniformize1D(density, npts, factor)
    else: fail = uniformize1D(density, npts, factor) # all zones

    if not fail:
        CTK.TXT.insert('START', 'Uniformize successfull.\n')
    else:
        CTK.TXT.insert('START', 'Uniformize edge fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
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
    
    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1: 
            fail = apply3D(1., 1, h, ntype=2)
        elif dim[2] != 1 and dim[3] == 1: 
            fail = apply2D(1., 1, h, ntype=2)
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

    type = VARS[4].get()
    factor = -1; npts = 2
    if type == 'Factor':
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

    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1: 
            fail = apply3D(1., npts, factor, ntype=1)
        elif dim[2] != 1 and dim[3] == 1: 
            fail = apply2D(1., npts, factor, ntype=1)
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

    # Get first selected zone
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(zone)
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1: 
            fail = apply3D(1., niter, eps, ntype=4)
        elif dim[2] != 1 and dim[3] == 1: 
            fail = apply2D(1., niter, eps, ntype=4)
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
    # Extrait la distribution en i
    source = D.getCurvilinearAbscissa(source)
    C._initVars(source, '{CoordinateX}={s}')
    C._initVars(source, 'CoordinateY', 0)
    C._initVars(source, 'CoordinateZ', 0)
    source = C.rmVars(source, 's')
    
    # Traitement
    if dim[0] == 'Structured':
        if dim[2] != 1 and dim[3] != 1: 
            fail = apply3D(1., 1, source, ntype=3)
        elif dim[2] != 1 and dim[3] == 1: 
            fail = apply2D(1., 1, source, ntype=3)
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
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkMapEdge', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Map distributions on edges.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=2)
    Frame.columnconfigure(3, weight=0)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkMapEdge')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Point density or Npts -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -1- Enforce height
    V = TK.StringVar(win); V.set('1.e-6'); VARS.append(V)
    if 'tkMapEdgeEnforceHeight' in CTK.PREFS:
        V.set(CTK.PREFS['tkMapEdgeEnforceHeight'])
    # -2- Option for uniformize
    V = TK.StringVar(win); V.set('Factor'); VARS.append(V)
    # -3- Source mesh for copy -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -4- Option for refine
    V = TK.StringVar(win); V.set('Factor'); VARS.append(V)
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

    # - Uniformize -
    B = TTK.Button(Frame, text="Uniformize", command=uniformize)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Uniformize an edge with regular spacing.')
    B = TTK.OptionMenu(Frame, VARS[2], 'Factor', 'Density', 'Npts')
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=7)
    B.grid(row=0, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Point multiplication factor/point density/number of points.')
    B.bind('<Return>', uniformize)

    # - Enforce -
    B = TTK.Button(Frame, text="Enforce", command=enforce)
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Enforce given spacing in edge.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=7)
    B.grid(row=1, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Enforced spacing.')
    B.bind('<Return>', enforce)

    # - Refine edge -
    B = TTK.Button(Frame, text="Refine", command=refine)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refine an edge keeping distribution.')
    B = TTK.OptionMenu(Frame, VARS[4], 'Factor', 'Npts')
    B.grid(row=2, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White', width=7)
    B.grid(row=2, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refinement factor or number of points.')
    B.bind('<Return>', refine)

    # - Copy distribution -
    B = TTK.Button(Frame, command=setSourceEdge,
                   image=iconics.PHOTO[8], padx=0, pady=0)
    B.grid(row=3, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set source edge defining distribution to copy.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=10)
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Source edge for distribution copy.')
    B = TTK.Button(Frame, text="Copy", command=copyDistrib)
    B.grid(row=3, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Copy distribution from source edge.')

    # - Smooth edge -
    B = TTK.Button(Frame, text="Smooth", command=smooth)
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Smooth an edge distribution.')
    B = TTK.Entry(Frame, textvariable=VARS[7], background='White', width=7)
    B.grid(row=4, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Smoother power.')
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White', width=7)
    B.grid(row=4, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of smoother iterations.')
    B.bind('<Return>', smooth)
    
#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW)

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

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
if (__name__ == "__main__"):
    import sys
    if (len(sys.argv) == 2):
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
