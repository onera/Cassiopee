# - fonctions de bloc split -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Transform.PyTree as T
import Post.PyTree as P
import Generator.PyTree as G

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# zip: connects zones if points at their borders are distant from tol
# IN: t, cplot.selectedZones, eps
# OUT: t with modified zones and displayed
#==============================================================================
def zip():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    eps = CTK.varsFromWidget(VARS[4].get(), type=1)
    if len(eps) != 1:
        CTK.TXT.insert('START', 'Zip tolerance is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    eps = eps[0]

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    fail = False
    zones = []; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        zones.append(z)
    try: zones = G.zip(zones, eps)
    except Exception as e:
        fail = True; errors += [0,str(e)]

    if not fail:
        c = 0
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            a = zones[c]; c += 1
            CTK.replace(CTK.t, nob, noz, a)
            CTK.TXT.insert('START', 'Zones zipped.\n')
    else:
        Panels.displayErrors(errors, header='Error: zip')
        CTK.TXT.insert('START', 'Zip fails at least for one zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()


#==============================================================================
# Join la selection
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche
#==============================================================================
def join():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    eps = CTK.varsFromWidget(VARS[4].get(), type=1)
    if len(eps) != 1:
        CTK.TXT.insert('START', 'Join tolerance is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    eps = eps[0]

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    Z = []; dels = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        Z.append(CTK.t[2][nob][2][noz])
        dels.append(CTK.t[2][nob][0]+Internal.SEP1+CTK.t[2][nob][2][noz][0])
    nob0 = CTK.Nb[nzs[0]]+1
    try:
        j = T.join(Z, tol=eps)
        CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
        CPlot.delete(dels)
        CTK.add(CTK.t, nob0, -1, j)
        CTK.TXT.insert('START', 'Selection joined.\n')
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: join')
        CTK.TXT.insert('START', 'Join failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#=============================================================================
def merge():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    eps = CTK.varsFromWidget(VARS[5].get(), type=1)
    if len(eps) != 1:
        CTK.TXT.insert('START', 'Merge tolerance is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    eps = eps[0]
    try: dircons = int(VARS[6].get())
    except: dircons = VARS[6].get()

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    Z = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        Z.append(CTK.t[2][nob][2][noz])
    nob0 = CTK.Nb[nzs[0]]+1

    if dircons == 'Names': # fait des paquets par noms homogenes
        zoneNames = {}
        zoneNames['__Pool__'] = []
        for z in Z:
            n = z[0]
            s = n.rsplit('.',2)
            if len(s) == 2:
                try: no = int(s[1])
                except: no = -1
                if no > -1:
                    if s[0] in zoneNames: zoneNames[s[0]] += [n]
                    else: zoneNames[s[0]] = [n]
                else: zoneNames['__Pool__'] += [n]
            else: zoneNames['__Pool__'] += [n]
        # Reconstitue les listes de zone
        out = []
        for k in zoneNames:
            l = []
            for i in zoneNames[k]:
                l.append(Internal.getNodeFromName(Z, i))
            out.append(l)
        try:
            js = []
            for zones in out:
                js += T.merge(zones, tol=eps)
            CTK.TXT.insert('START', 'Selection merged.\n')
        except Exception as e:
            Panels.displayErrors([0,str(e)], header='Error: merge')
            CTK.TXT.insert('START', 'Merge failed.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')

        CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
        CTK.t[2][nob0][2] += js

    else:
        try:
            j = T.merge(Z, dir=dircons, tol=eps)
            CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
            CTK.t[2][nob0][2] += j
            CTK.TXT.insert('START', 'Selection merged.\n')
        except Exception as e:
            Panels.displayErrors([0,str(e)], header='Error: merge')
            CTK.TXT.insert('START', 'Merge failed.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def split():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    stype = VARS[1].get()
    CTK.saveTree()

    # View split
    if stype == 'X-view':
        import KCore.Vector as Vector
        point = CPlot.getActivePoint()
        if point == []: return
        posCam = CPlot.getState('posCam')
        posEye = CPlot.getState('posEye')
        dirCam = CPlot.getState('dirCam')
        v1 = dirCam
        v2 = Vector.sub(posEye, posCam)
        v3 = Vector.cross(v1, v2)
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            zp = T.rotate(z, point, (v1,v2,v3), ((1,0,0), (0,1,0), (0,0,1)))
            zp = C.initVars(zp, '{centers:__tag__}={centers:CoordinateZ}<=%f'%point[2])
            #zp = T.rotate(zp, point, ((1,0,0), (0,1,0), (0,0,1)), (v1,v2,v3))
            # replace initial coordinates
            gc1 = Internal.getNodeFromName(zp, Internal.__GridCoordinates__)
            gc2 = Internal.getNodeFromName(z, Internal.__GridCoordinates__)
            gc1[2][:] = gc2[2][:]

            z1 = P.selectCells2(zp, 'centers:__tag__')
            z1[0] = C.getZoneName(z[0]+'1')
            C._initVars(zp, '{centers:__tag__} = 1.-{centers:__tag__}')
            z2 = P.selectCells2(zp, 'centers:__tag__')
            z2[0] = C.getZoneName(z[0]+'2')
            CTK.replace(CTK.t, nob, noz, z1)
            CTK.add(CTK.t, nob, -1, z2)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
        return

    # Coordinate Split
    if stype == 'X-coords' or stype == 'Y-coords' or stype == 'Z-coords':
        point = CPlot.getActivePoint()
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            if stype == 'X-coords':
                zp = C.initVars(z, '{centers:__tag__}={centers:CoordinateX}<=%f'%point[0])
            elif stype == 'Y-coords':
                zp = C.initVars(z, '{centers:__tag__}={centers:CoordinateY}<=%f'%point[1])
            else:
                zp = C.initVars(z, '{centers:__tag__}={centers:CoordinateZ}<=%f'%point[2])
            z1 = P.selectCells2(zp, 'centers:__tag__')
            z1[0] = C.getZoneName(z[0]+'1')
            C._initVars(zp, '{centers:__tag__} = 1.-{centers:__tag__}')
            z2 = P.selectCells2(zp, 'centers:__tag__')
            z2[0] = C.getZoneName(z[0]+'2')
            CTK.replace(CTK.t, nob, noz, z1)
            CTK.add(CTK.t, nob, -1, z2)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
        return

    # Indice splits
    ind = CPlot.getActivePointIndex()
    if ind == []:
        CTK.TXT.insert('START', 'No active point.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
        return
    i1 = ind[2]; j1 = ind[3]; k1 = ind[4]

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dims = Internal.getZoneDim(z)
        try:
            if dims[0] == 'Structured':
                ni = dims[1]; nj = dims[2]; nk = dims[3]
                if stype == 'i-indices':
                    z1 = T.subzone(z, (1,1,1), (i1,nj,nk))
                    z2 = T.subzone(z, (i1,1,1), (ni,nj,nk))
                elif stype == 'j-indices':
                    z1 = T.subzone(z, (1,1,1), (ni,j1,nk))
                    z2 = T.subzone(z, (1,j1,1), (ni,nj,nk))
                else:
                    z1 = T.subzone(z, (1,1,1), (ni,nj,k1))
                    z2 = T.subzone(z, (1,1,k1), (ni,nj,nk))
                CTK.replace(CTK.t, nob, noz, z1)
                CTK.add(CTK.t, nob, -1, z2)
            elif dims[0] == 'Unstructured' and dims[3] == 'BAR':
                B = T.splitBAR(z, ind[0])
                if len(B) == 1: CTK.replace(CTK.t, nob, noz, B[0])
                else:
                    CTK.replace(CTK.t, nob, noz, B[0])
                    CTK.add(CTK.t, nob, -1, B[1])
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if fail:
        CTK.TXT.insert('START', 'split fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    else:
        Panels.displayErrors(errors, header='Error: split')
        CTK.TXT.insert('START', 'split done.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# splitSize
# IN: t, cplot.selectedZones, N : taille max au-dessus de laquelle on split
# OUT: t modifie et affiche
#==============================================================================
def splitSize(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    try: size = int(VARS[0].get())
    except:
        CTK.TXT.insert('START', 'splitSize: invalid specified size.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    level = int(VARS[2].get())

    CTK.saveTree()
    fail = False; errors = []
    if nzs == []:
        try:
            CTK.t = T.splitSize(CTK.t, size, multigrid=level)
            CTK.display(CTK.t)
        except Exception as e:
            fail = True; errors += [0,str(e)]
    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            try:
                blocks = T.splitSize(CTK.t[2][nob][2][noz], size,
                                     multigrid=level)
                CTK.replace(CTK.t, nob, noz, blocks[0])
                for i in blocks[1:]: CTK.add(CTK.t, nob, -1, i)
            except Exception as e:
                fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'splitSize done.\n')
    else:
        Panels.displayErrors(errors, header='Error: splitSize')
        CTK.TXT.insert('START', 'Splitsize fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# splitNParts
# IN: t, cplot.selectedZones, N: nb de parties
# OUT: t modifie et affiche
#==============================================================================
def splitNParts(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    try: NParts = int(VARS[7].get())
    except:
        CTK.TXT.insert('START', 'splitNParts: invalid specified number.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    level = int(VARS[2].get())

    CTK.saveTree()
    fail = False; errors = []
    if nzs == []:
        try:
            CTK.t = T.splitNParts(CTK.t, NParts, multigrid=level)
            CTK.display(CTK.t)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    else:
        sel = []; nobs = {}; pname = {}
        # change le nom des zones dans la selection
        c = 0
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            zone = CTK.t[2][nob][2][noz]
            sel.append(zone)
            nobs['__z_%d__'%c] = (nob, noz)
            pname['__z_%d__'%c] = zone[0]
            c += 1
        sel = Internal.copyRef(sel)
        c = 0
        for s in sel:
            s[0] = '__z_%d__'%c; c += 1
        try:
            blocks = T.splitNParts(sel, NParts, multigrid=level)
            # match blocks by parent name
            rep = {}
            for b in blocks:
                name = b[0]
                name = name.rsplit('.', 1)[0]
                b[0] = C.getZoneName(pname[name])
                nob = nobs[name][0]; noz = nobs[name][1]
                if name not in rep: # remplace le parent
                    CTK.replace(CTK.t, nob, noz, b); rep[name] = 1
                else: CTK.add(CTK.t, nob, -1, b)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'splitNParts done.\n')
    else:
        Panels.displayErrors(errors, header='Error: splitNParts')
        CTK.TXT.insert('START', 'splitNParts fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# splitMultiplePoints: decoupe les blocs au niveau des raccords triples
# IN: t
# OUT: t modifie et affiche
#==============================================================================
def splitMultiplePoints():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: ndim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); ndim = 3

    CTK.saveTree()
    try:
        CTK.t = T.splitMultiplePts(CTK.t, dim=ndim)
        CTK.TXT.insert('START', 'splitMultiplePts done on all tree.\n')
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: splitMultiplePoints')
        CTK.TXT.insert('START', 'Split multiple points fails.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# splitConnexity
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche
#==============================================================================
def splitConnexity():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.setCursor(2, WIDGETS['splitConnexity'])
    CTK.saveTree()
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = G.close(z)
        try:
            splits = T.splitConnexity(z)
            CTK.replace(CTK.t, nob, noz, splits[0])
            for i in splits[1:]: CTK.add(CTK.t, nob, -1, i)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'splitConnexity done.\n')
    else:
        Panels.displayErrors(errors, header='Error: splitConnexity')
        CTK.TXT.insert('START', 'Split connexity fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['splitConnexity'])

#==============================================================================
# splitManifold
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche
#==============================================================================
def splitManifold():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.setCursor(2, WIDGETS['splitManifold'])
    CTK.saveTree()
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = G.close(z)
        try:
            splits = T.splitManifold(z)
            CTK.replace(CTK.t, nob, noz, splits[0])
            for i in splits[1:]: CTK.add(CTK.t, nob, -1, i)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'splitManifold done.\n')
    else:
        Panels.displayErrors(errors, header='Error: splitManifold')
        CTK.TXT.insert('START', 'Split manifold fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['splitManifold'])

# Split suivant les sharpEdges + distance
def splitSharpEdgesWithDelta__(a, alphaRef, deltaRef):
    try: import Post.PyTree as P; import Dist2Walls.PyTree as Dist2Walls
    except: raise ImportError("splitSharpEdges: requires Post, Dist2Walls modules.")
    e = P.sharpEdges(a, alphaRef)
    b = Dist2Walls.distance2Walls(a, e, loc='centers')
    a1 = P.selectCells(b, '{centers:TurbulentDistance}<%g'%deltaRef)
    a2 = P.selectCells(b, '{centers:TurbulentDistance}>=%g'%deltaRef)
    return [a1,a2]

#=========================================================================
# splitSharpAngles
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche
#=========================================================================
def splitSharpAngles():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    alphaRef = CTK.varsFromWidget(VARS[3].get(), type=1)
    if len(alphaRef) != 1:
        CTK.TXT.insert('START', 'Split angle is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    alphaRef = alphaRef[0]
    deltaRef = CTK.varsFromWidget(VARS[8].get(), type=1)
    if len(deltaRef) != 1:
        CTK.TXT.insert('START', 'Distance to sharp edges is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    deltaRef = deltaRef[0]

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = G.close(z)
        try:
            if deltaRef <= 0.: splits = T.splitSharpEdges(z, alphaRef)
            else: splits = splitSharpEdgesWithDelta__(z, alphaRef, deltaRef)
            CTK.replace(CTK.t, nob, noz, splits[0])
            for i in splits[1:]: CTK.add(CTK.t, nob, -1, i)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'splitSharpEdges done.\n')
    else:
        Panels.displayErrors(errors, header='Error: splitSharpEdges')
        CTK.TXT.insert('START', 'SplitSharpEdges fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkSplit  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Split meshes.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkSplit')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- splitSize -
    V = TK.StringVar(win); V.set('1000'); VARS.append(V)
    if 'tkSplitSize' in CTK.PREFS: V.set(CTK.PREFS['tkSplitSize'])
    # -1- direction pour subzone
    V = TK.StringVar(win); V.set('i-indices'); VARS.append(V)
    # -2- multigrid (niveau)
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkSplitSizeMultigrid' in CTK.PREFS:
        V.set(CTK.PREFS['tkSplitSizeMultigrid'])
    # -3- splitSharpEdges angle -
    V = TK.StringVar(win); V.set('30.'); VARS.append(V)
    if 'tkSplitSharpEdges' in CTK.PREFS:
        V.set(CTK.PREFS['tkSplitSharpEdges'])
    # -4- tol for join -
    V = TK.StringVar(win); V.set('1.e-8'); VARS.append(V)
    if 'tkSplitJoinTol' in CTK.PREFS:
        V.set(CTK.PREFS['tkSplitJoinTol'])
    # -5- tol for merge -
    V = TK.StringVar(win); V.set('1.e-8'); VARS.append(V)
    if 'tkSplitMergeTol' in CTK.PREFS:
        V.set(CTK.PREFS['tkSplitMergeTol'])
    # -6- constraint direction for merge -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -7- Nparts for splitNParts -
    V = TK.StringVar(win); V.set('2'); VARS.append(V)
    if 'tkSplitNParts' in CTK.PREFS: V.set(CTK.PREFS['tkSplitNParts'])
    # -8- Multigrid level for splitNParts -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkSplitNPartsMultigrid' in CTK.PREFS:
        V.set(CTK.PREFS['tkSplitNPartsMultigrid'])
    # -9- distance for splitSharpEdges
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)

    # - Buttons -
    # - split (subzone) -
    B = TTK.Button(Frame, text="Split", command=split)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split blocks following\n index direction (structured)\nor geometric direction.')
    B = TTK.OptionMenu(Frame, VARS[1], 'i-indices', 'j-indices', 'k-indices',
                       'X-coords', 'Y-coords', 'Z-coords',
                       'X-view')
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)

    # - splitSize -
    B = TTK.Button(Frame, text="SplitSize", command=splitSize)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split a mesh with respect to a maximum number of points.\nTree is modified.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=10)
    B.grid(row=1, column=1, sticky=TK.EW)
    B.bind('<Return>', splitSize)
    BB = CTK.infoBulle(parent=B, text='Max number of points.')
    B = TTK.OptionMenu(Frame, VARS[2], '0', '1', '2')
    B.grid(row=1, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Multigrid level.')

    # - splitNParts -
    B = TTK.Button(Frame, text="SplitNParts", command=splitNParts)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split a mesh in N Parts.\nTree is modified.')
    B = TTK.Entry(Frame, textvariable=VARS[7], background='White', width=10)
    B.grid(row=2, column=1, sticky=TK.EW)
    B.bind('<Return>', splitNParts)
    BB = CTK.infoBulle(parent=B, text='Number of parts.')
    B = TTK.OptionMenu(Frame, VARS[8], '0', '1', '2')
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Multigrid level.')

    # - splitManifold -
    B = TTK.Button(Frame, text="SplitManifold", command=splitManifold)
    WIDGETS['splitManifold'] = B
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split a mesh into manifold parts.\nTree is modified.')

    # - splitConnexity -
    B = TTK.Button(Frame, text="SplitConnexity", command=splitConnexity)
    WIDGETS['splitConnexity'] = B
    B.grid(row=3, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split a mesh into connex parts.\nTree is modified.')

    # - splitMultiplePoints -
    B = TTK.Button(Frame, text="MP", command=splitMultiplePoints)
    B.grid(row=3, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split a mesh at corners between more that 2 blocks.\nTree is modified.')

    # - splitSharpEdges -
    B = TTK.Button(Frame, text="SplitSharpEdges", command=splitSharpAngles)
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split a mesh into smooth parts.\nTree is modified.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=4)
    B.grid(row=4, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Angle of split.')
    B = TTK.Entry(Frame, textvariable=VARS[8], background='White', width=4)
    B.grid(row=4, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Distance to sharp edges.')

    # - Join selection -
    B = TTK.Button(Frame, text="Join", command=join)
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Join selection into a single zone.\nReplaced in tree.')
    B = TTK.Button(Frame, text="Zip", command=zip)
    B.grid(row=5, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Zip zones.\nTree is modified.')
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=4)
    B.grid(row=5, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Join/zip tolerance.')

    # - Merge surface grids by Rigby -
    B = TTK.Button(Frame, text="Merge", command=merge)
    B.grid(row=6, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Merge blocks if possible.\nTree is modified.')
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White', width=10)
    B.grid(row=6, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Merge tolerance.')
    B = TTK.OptionMenu(Frame, VARS[6], '0', '1', '2', '3', 'Names')
    B.grid(row=6, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Favorite merging direction.\n0: possibly merge in all directions,\n1: favor merging in i direction, ...')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['BlockNoteBook'].add(WIDGETS['frame'], text='tkSplit')
    except: pass
    CTK.WIDGETS['BlockNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['BlockNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkSplitSize'] = VARS[0].get()
    CTK.PREFS['tkSplitSizeMultigrid'] = VARS[2].get()
    CTK.PREFS['tkSplitSharpEdges'] = VARS[3].get()
    CTK.PREFS['tkSplitJoinTol'] = VARS[4].get()
    CTK.PREFS['tkSplitMergeTol'] = VARS[5].get()
    CTK.PREFS['tkSplitNParts'] = VARS[7].get()
    CTK.PREFS['tkSplitNPartsMultigrid'] = VARS[8].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('1000')
    VARS[2].set('0')
    VARS[3].set('30.')
    VARS[4].set('1.e-8')
    VARS[5].set('1.e-8')
    VARS[7].set('2')
    VARS[8].set('0')
    CTK.PREFS['tkSplitSize'] = VARS[0].get()
    CTK.PREFS['tkSplitSizeMultigrid'] = VARS[2].get()
    CTK.PREFS['tkSplitSharpEdges'] = VARS[3].get()
    CTK.PREFS['tkSplitJoinTol'] = VARS[4].get()
    CTK.PREFS['tkSplitMergeTol'] = VARS[5].get()
    CTK.PREFS['tkSplitNParts'] = VARS[7].get()
    CTK.PREFS['tkSplitNPartsMultigrid'] = VARS[8].get()
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
    (win, menu, file, tools) = CTK.minimal('tkSplit '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
