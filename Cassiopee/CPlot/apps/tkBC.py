# -- tkBCs --
"""Applet to view/set BCs in a pyTree."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.PyTree as X
import Post.PyTree as P

# local widgets list
WIDGETS = {}; VARS = []
__SPLITFACTOR__ = 180.

#==============================================================================
def viewl(v, l):
    v.set(l); view()

#==============================================================================
# Retourne les BCs definies dans t + les familles
def getAllDefinedBC(t):
    natives = set()
    zones = Internal.getZones(t)

    # FamilyBC
    FamilyBC = C.getFamilyBCNamesDict(t)

    # Defined BC
    for z in zones:
        nodes = Internal.getNodesFromType2(z, 'BC_t')
        for i in nodes:
            if Internal.getValue(i) != 'FamilySpecified':
                natives.add(Internal.getValue(i))
            else:
                f = Internal.getNodeFromType1(i, 'FamilyName_t')
                if f is not None:
                    name = Internal.getValue(f)
                    if name in FamilyBC: natives.add(FamilyBC[name])

    for z in zones:
        # BCMatch
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        for i in nodes:
            if i is not None: natives.add('BCMatch')
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for i in nodes:
            r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
            if r is not None:
                if Internal.getValue(r) == 'Abutting1to1': natives.add('BCMatch')

        # BCNearMatch
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for i in nodes:
            r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
            if r is not None:
                if Internal.getValue(r) == 'Abutting':
                    f = Internal.getNodeFromType1(i, 'FamilyName_t')
                    if f is not None and Internal.getValue(f).startswith('BCStage'):
                        natives.add(Internal.getValue(f))
                    else: natives.add('BCNearMatch')

        # BCOverlap
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for i in nodes:
            r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
            if r is not None:
                if Internal.getValue(r) == 'Overset': natives.add('BCOverlap')

    natives = list(natives)
    natives.sort(key=str.lower)

    # FamilyBC
    fams = []
    for i in FamilyBC: fams.append(i)
    fams.sort(key=str.lower)
    return natives + fams

#==============================================================================
# Pour view BC
def updateFamilyBCNameList1(event=None):
    if CTK.t == []: return
    m = WIDGETS['BCs'].children['menu']
    m.delete(0, TK.END)
    varsp = ['Mesh', '-All BC-']+getAllDefinedBC(CTK.t)
    for i in varsp:
        m.add_command(label=i, command=lambda v=VARS[0],l=i:viewl(v,l))

def updateFamilyBCNameList1_2(event=None):
    if CTK.t == []: return
    varsp = ['Mesh', '-All BC-']+getAllDefinedBC(CTK.t)
    if 'BCs' in WIDGETS: WIDGETS['BCs']['values'] = varsp

#================================================================================
# Pour list box
def updateBCNameList(event=None):
    if CTK.t == []: return

    lb = WIDGETS['BCLB']
    lb.focus_set()
    varsbc = ['-All BC-']+getAllDefinedBC(CTK.t)
    if len(varsbc) == lb.size(): return # un peu brutal

    # Remplace tous les elements
    lb.delete(0, TK.END)
    for i, value in enumerate(['-All BC-']+getAllDefinedBC(CTK.t)): lb.insert(i, value)
    return lb

#==============================================================================
# Pour rm BC
def updateFamilyBCNameList2(event=None):
    if CTK.t == []: return
    m = WIDGETS['BCs3'].children['menu']
    m.delete(0, TK.END)
    varsp = ['-All BC-']+getAllDefinedBC(CTK.t)
    for i in varsp:
        m.add_command(label=i, command=lambda v=VARS[5],l=i:v.set(l))

def updateFamilyBCNameList2_2(event=None):
    if CTK.t == []: return
    varsp = ['-All BC-']+getAllDefinedBC(CTK.t)
    if 'BCs3' in WIDGETS: WIDGETS['BCs3']['values'] = varsp

#==============================================================================
# Pour set BC
def updateFamilyBCNameList3(event=None):
    if CTK.t == []: return
    varsl = C.getFamilyBCNamesOfType(CTK.t)
    m = WIDGETS['BCs2'].children['menu']
    m.delete(0, TK.END)
    varsp = Internal.KNOWNBCS[:]
    if len(varsl) != 0: varsp += varsl
    for i in varsp:
        m.add_command(label=i, command=lambda v=VARS[6],l=i:v.set(l))

def updateFamilyBCNameList3_2(event=None):
    if CTK.t == []: return
    varsl = C.getFamilyBCNamesOfType(CTK.t)
    varsp = Internal.KNOWNBCS[:]
    if len(varsl) != 0:
        varsl.sort(key=str.lower)
        varsp += varsl
    if 'BCs2' in WIDGETS: WIDGETS['BCs2']['values'] = varsp
#==============================================================================
# Pour fillEmptyBC
def updateFamilyBCNameList4(event=None):
    if CTK.t == []: return
    varsl = C.getFamilyBCNamesOfType(CTK.t)
    m = WIDGETS['BCs4'].children['menu']
    m.delete(0, TK.END)
    varsp = Internal.KNOWNBCS[:]
    if len(varsl) != 0: varsp += varsl
    for i in varsp:
        m.add_command(label=i, command=lambda v=VARS[6],l=i:v.set(l))

def updateFamilyBCNameList4_2(event=None):
    if CTK.t == []: return
    varsl = C.getFamilyBCNamesOfType(CTK.t)
    varsp = Internal.KNOWNBCS[:]
    if len(varsl) != 0:
        varsl.sort(key=str.lower)
        varsp += varsl
    if 'BCs4' in WIDGETS: WIDGETS['BCs4']['values'] = varsp

#==============================================================================
# DisplayUndefinedBoundaries
#==============================================================================
def displayUndefinedBoundaries():
    CTK.TXT.insert('START', 'Display undefined boundaries.\n')
    check()

#==============================================================================
# Rebascule sur une vue dde CTK.t
#==============================================================================
def viewMesh(event=None):
    CTK.display(CTK.t)

#==============================================================================
# Construit dt en fonction du type de BC demandee
# N'agit que sur les zones CPlot-actives
# IN: t
# IN: VARS[0] state (type de BC demandee)
# OUT: dt construit et affiche
#==============================================================================
def view(event=None):
    if CTK.t == []: return

    BCTypes = []
    selection = WIDGETS['BCLB'].curselection()
    for s in selection:
        t = WIDGETS['BCLB'].get(s)
        if t not in Internal.KNOWNBCS: t = 'FamilySpecified:'+t
        BCTypes.append(t)
    if 'FamilySpecified:-All BC-' in BCTypes: BCTypes = ['*']

    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()

    tp = Internal.appendBaseName2ZoneName(CTK.t, updateRef=False,
                                          separator=Internal.SEP1)
    CTK.dt = C.newPyTree(['Base', 'Edges'])
    active = []
    for z in CTK.__MAINACTIVEZONES__: active.append(tp[2][CTK.Nb[z]+1][2][CTK.Nz[z]])

    Z = []
    for t in BCTypes:
        Z += C.extractBCOfType(active, t, topTree=tp)
        if t == 'BCWall': # Dans ce cas, affiche tous les types de BCWall
            Z += C.extractBCOfType(active, 'BCWallInviscid')
            Z += C.extractBCOfType(active, 'BCWallViscous')
            Z += C.extractBCOfType(active, 'BCWallViscousIsoThermal')
    CTK.dt[2][1][2] += Z

    if VARS[7].get() == '1': # display les edges des zones en +
        exts = []
        for z in active:
            ztype = Internal.getZoneType(z)
            if ztype == 1:
                zp = P.exteriorFacesStructured(z)
                exts += zp
            else:
                #zp = P.exteriorFaces(z)
                #zp = P.sharpEdges(zp)
                zp = []
                exts += zp

        CTK.dt[2][2][2] += exts
        C._fillMissingVariables(CTK.dt) # bug exteriorFaces

        # Activate
        lenZ = len(CTK.dt[2][1][2]); lenExts = len(CTK.dt[2][2][2])
        active = [(i,1) for i in range(lenZ+lenExts)]
        for i in range(lenZ): active[i] = (i,1)
        for i in range(lenExts): active[i+lenZ] = (i+lenZ,0)

        CTK.display(CTK.dt, mainTree=CTK.DEFINEDBC)
        CPlot.setActiveZones(active)
        CPlot.setState(edgifyDeactivatedZones=1)
    else:
        lenZ = len(CTK.dt[2][1][2])
        active = [(i,1) for i in range(lenZ)]
        C._fillMissingVariables(CTK.dt) # si BCDataSet != fields
        CTK.display(CTK.dt, mainTree=CTK.DEFINEDBC)
        CPlot.setActiveZones(active)
        CPlot.setState(edgifyDeactivatedZones=0)

#==============================================================================
# Construit dt contenant les BC non definies
# IN: t
# IN: dimension du probleme
# IN: splitFactor: slider
# OUT: dt construit et affiche
#==============================================================================
def check():
    if CTK.t == []: return
    if CTK.CADHOOK is not None: return

    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: ndim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); ndim = 3

    # Varie de 0 a 180 degres
    global __SPLITFACTOR__
    splitFactor = 180.-WIDGETS['splitFactor'].get()*180./100.
    __SPLITFACTOR__ = splitFactor

    wins = C.getEmptyBC(CTK.t, ndim, splitFactor)
    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()
    CTK.dt = C.newPyTree(['Base', 'Edges'])
    tp = Internal.appendBaseName2ZoneName(CTK.t, updateRef=False,
                                          separator=Internal.SEP1,
                                          trailing=Internal.SEP1)
    bases = Internal.getBases(tp)
    nb = 0
    for b in bases:
        nodes = Internal.getNodesFromType1(b, 'Zone_t')
        nz = 0
        for z in nodes:
            ztype = Internal.getZoneType(z)
            winz = wins[nb][nz]
            if ztype == 1: # structure
                for w in winz:
                    imin = w[0]; imax = w[1]
                    jmin = w[2]; jmax = w[3]
                    kmin = w[4]; kmax = w[5]
                    zp = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax))
                    CTK.dt[2][1][2].append(zp)
            else: # non structure
                for w in winz:
                    zp = T.subzone(z, w, type='faces')
                    CTK.dt[2][1][2].append(zp)
            nz += 1
        nb += 1

    if VARS[7].get() == '1': # display les edges des zones en +
        exts = []
        zones = Internal.getZones(tp)
        for z in zones:
            ztype = Internal.getZoneType(z)
            if ztype == 1:
                zp = P.exteriorFacesStructured(z)
                exts += zp
            else:
                #zp = P.exteriorFaces(z); zp = P.sharpEdges(zp)
                zp = []
                exts += zp
        CTK.dt[2][2][2] += exts
        #C._fillMissingVariables(CTK.dt) # bug exteriorFacesStruct

        # Activate
        lenZ = len(CTK.dt[2][1][2]); lenExts = len(exts)
        active = [(i,1) for i in range(lenZ+lenExts)]
        for i in range(lenZ): active[i] = (i,1)
        for i in range(lenExts): active[i+lenZ] = (i+lenZ,0)

        CTK.display(CTK.dt, mainTree=CTK.UNDEFINEDBC)
        CPlot.setActiveZones(active)
        CPlot.setState(edgifyDeactivatedZones=1)
    else:
        lenZ = len(CTK.dt[2][1][2])
        active = [(i,1) for i in range(lenZ)]
        CTK.display(CTK.dt, mainTree=CTK.UNDEFINEDBC)
        CPlot.setActiveZones(active)

    # modifie la couleur du bouton
    bases = Internal.getBases(CTK.dt)
    if len(bases) > 0: l = len(Internal.getZones(bases[0]))
    else: l = 0
    if l == 0: TTK.setButtonGreen(WIDGETS['undefinedBC'])
    else: TTK.setButtonRed(WIDGETS['undefinedBC'])
    WIDGETS['undefinedBC'].update()

#==============================================================================
# Builds the BCDegenerateLine and BCDegeneratePoint for selected zones
# IN: t
# IN: VARS[2] state (eps)
# OUT: t displayed per view
#==============================================================================
def setDegeneratedBC():
    if CTK.t == []: return
    eps = VARS[2].get()
    eps = float(eps)
    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: ndim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); ndim = 3

    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        try:
            CTK.t = X.setDegeneratedBC(CTK.t, tol=eps, dim=ndim)
            CTK.TKTREE.updateApp()
            CTK.TXT.insert('START', 'Degenerated BCs successfully set.\n')
        except Exception as e:
            Panels.displayErrors([0,str(e)], header='Error: setDegenratedBC')
            CTK.TXT.insert('START', 'Degenerated BCs failed.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    else:
        sel = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            sel.append(z)
        try:
            sel = X.setDegeneratedBC(sel, tol=eps, dim=ndim)
            c = 0
            for nz in nzs:
                nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
                CTK.t[2][nob][2][noz] = sel[c]; c += 1
            CTK.TKTREE.updateApp()
            CTK.TXT.insert('START', 'Degenerated BCs successfully set.\n')
        except Exception as e:
            Panels.displayErrors([0,str(e)], header='Error: setDegeneratedBC')
            CTK.TXT.insert('START', 'Degenerated BCs failed.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    check()

#==============================================================================
# Construit les BC Match dans t pour les zones selectionnees
# IN: t
# IN: VARS[2] state (eps)
# OUT: t affiche par view
#==============================================================================
def connectMatch():
    if CTK.t == []: return
    eps = VARS[2].get()
    eps = float(eps)
    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: ndim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); ndim = 3

    mode = VARS[9].get()
    translation = [0.,0.,0.]; rotationCenter = [0.,0.,0.]; rotationAngle = [0.,0.,0.]
    if mode == 'Translation':
        mode = 1; translation = CTK.varsFromWidget(VARS[10].get(), type=1)
    elif mode == 'Rotation (Degree)':
        mode = 1; f = CTK.varsFromWidget(VARS[10].get(), type=1)
        rotationCenter = f[0:3]; rotationAngle = f[3:6]
    else: mode = 0

    CTK.setCursor(2, WIDGETS['connectMatch'])
    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        try:
            if mode == 0: CTK.t = X.connectMatch(CTK.t, tol=eps, dim=ndim)
            else:
                CTK.t = X.connectMatchPeriodic(CTK.t, rotationCenter=rotationCenter,
                                               rotationAngle=rotationAngle,
                                               translation=translation, tol=eps, dim=ndim,
                                               unitAngle=None)
            CTK.TKTREE.updateApp()
            CTK.TXT.insert('START', 'Matching BCs successfully set.\n')
        except Exception as e:
            Panels.displayErrors([0,str(e)], header='Error: connectMatch')
            CTK.TXT.insert('START', 'Matching BCs failed.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    else:
        sel = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            sel.append(z)
        try:
            if mode == 0: sel = X.connectMatch(sel, tol=eps, dim=ndim)
            else:
                sel = X.connectMatchPeriodic(sel, rotationCenter=rotationCenter,
                                             rotationAngle=rotationAngle,
                                             translation=translation, tol=eps,
                                             dim=ndim, unitAngle=None)
            c = 0
            for nz in nzs:
                nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
                CTK.t[2][nob][2][noz] = sel[c]; c += 1
            CTK.TKTREE.updateApp()
            CTK.TXT.insert('START', 'Matching BCs successfully set.\n')
        except Exception as e:
            Panels.displayErrors([0,str(e)], header='Error: connectMatch')
            CTK.TXT.insert('START', 'Matching BCs failed.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    check()
    CTK.setCursor(0, WIDGETS['connectMatch'])

#==============================================================================
# Construit les BC NearMatch dans t pour les zones selectionnees
# IN: t
# IN: VARS[3] state (eps)
# OUT: t affiche par view
#==============================================================================
def connectNearMatch():
    if CTK.t == []: return
    eps = VARS[2].get(); eps = float(eps)
    ratio = VARS[3].get()
    ratio = CTK.varsFromWidget(ratio, type=2)

    if len(ratio) != 1 and len(ratio) != 3:
        CTK.TXT.insert('START', 'Invalid ratio for nearmatch.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if len(ratio) == 1: ratio = ratio[0]

    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: ndim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); ndim = 3

    CTK.setCursor(2, WIDGETS['connectNearMatch'])
    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        try:
            CTK.t = X.connectNearMatch(CTK.t, ratio=ratio, tol=eps, dim=ndim)
            CTK.TKTREE.updateApp()
            CTK.TXT.insert('START', 'n/m matching BCs successfully set.\n')
        except Exception as e:
            Panels.displayErrors([0,str(e)], header='Error: connectNearMatch')
            CTK.TXT.insert('START', 'n/m matching BCs failed.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    else:
        sel = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            sel.append(z)
        try:
            sel = X.connectNearMatch(sel, ratio=ratio, tol=eps, dim=ndim)
            c = 0
            for nz in nzs:
                nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
                CTK.t[2][nob][2][noz] = sel[c]; c += 1
            CTK.TKTREE.updateApp()
            CTK.TXT.insert('START', 'n/m matching BCs successfully set.\n')
        except Exception as e:
            Panels.displayErrors([0,str(e)], header='Error: connectNearMatch')
            CTK.TXT.insert('START', 'n/m matching BCs failed.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    check()
    CTK.setCursor(0, WIDGETS['connectNearMatch'])

#==============================================================================
# Remplit les BC vide de t avec le type demande
# IN: t
# IN: VARS[4] state (type de BC)
# OUT: t construit et affiche en fonction du filtre view BC
#==============================================================================
def fillEmptyBCWith():
    if CTK.t == []: return
    if CTK.CADHOOK is not None: return # not coded for now

    typeBC = VARS[4].get()
    if typeBC not in Internal.KNOWNBCS:
        nameBC = typeBC; typeBC = 'FamilySpecified:'+typeBC
    else: nameBC = typeBC

    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: ndim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); ndim = 3

    CTK.saveTree()
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        C._fillEmptyBCWith(CTK.t, nameBC+'Fill', typeBC, dim=ndim)
    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            C._fillEmptyBCWith(CTK.t[2][nob][2][noz],
                               nameBC+'Fill', typeBC, dim=ndim)
    CTK.TXT.insert('START', 'Empty BCs filled.\n')
    CTK.TKTREE.updateApp()
    check()

#==============================================================================
# Enleve les BC de t du type demande
# IN: t
# IN: VARS[5] state (type de BC)
# OUT: t construit et affiche en fonction du filtre view BC
#==============================================================================
def rmBCOfType():
    if CTK.t == []: return

    BCTypes = []
    selection = WIDGETS['BCLB'].curselection()
    for s in selection:
        t = WIDGETS['BCLB'].get(s)
        if t not in Internal.KNOWNBCS: t = 'FamilySpecified:'+t
        BCTypes.append(t)

    CTK.saveTree()
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        if 'FamilySpecified:-All BC-' in BCTypes:
            Internal._rmNodesByType(CTK.t, 'BC_t')
        else:
            for t in BCTypes: C._rmBCOfType(CTK.t, t)
    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            if 'FamilySpecified:-All BC-' in BCTypes:
                Internal._rmNodesByType(z, 'BC_t')
            else:
                for t in BCTypes: C._rmBCOfType(z, t)
    if len(BCTypes) > 0:
        CTK.TXT.insert('START', 'BCs of type %s have been removed.\n'%BCTypes[0])
    CTK.TKTREE.updateApp()
    check()

#==============================================================================
# set BC avec le type demande
# IN: t, BC selectionnee dans dt undef
# IN: VARS[6] state (type de BC)
# OUT: t construit et affiche en fonction du filtre view BC
#==============================================================================
def setBCWith():
    if CTK.t == []: return

    nzs = CPlot.getSelectedZones()
    if nzs == []: return
    typeBC = VARS[6].get()
    if typeBC not in Internal.KNOWNBCS:
        nameBC = typeBC; typeBC = 'FamilySpecified:'+typeBC
    else: nameBC = typeBC

    # On CAD tree, we set directly BC to full zone
    if CTK.CADHOOK is not None:
        CTK.saveTree()
        CTK.setCursor(2, WIDGETS['setBCWith'])
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            dim = Internal.getZoneDim(z)
            C._addBC2Zone(z, nameBC, typeBC, elementRange=[0,dim[2]])
        CTK.TXT.insert('START', 'BCs set to %s.\n'%typeBC)
        CTK.TKTREE.updateApp()
        CTK.setCursor(0, WIDGETS['setBCWith'])
        return

    # On volume tree, we set from undefinedbc
    if CTK.__MAINTREE__ != CTK.UNDEFINEDBC:
        CTK.TXT.insert('START', 'Fail on a this tree (view undefined BC before).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: ndim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); ndim = 3

    splitFactor = 180.-WIDGETS['splitFactor'].get()*180./100.
    if __SPLITFACTOR__ != splitFactor:
        CTK.TXT.insert('START', 'Split factor changed: view undefinedBC again.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['setBCWith'])

    wins = C.getEmptyBC(CTK.t, ndim, splitFactor)

    # calcul nzr : le nz dans la numerotation des fenetres emptys
    nzr = {}
    c = 0; cs = 0
    for wb in wins:
        for wz in wb:
            for w in wz:
                if isinstance(w, list): nzr[c] = cs; cs += 1
                c += 1
    c = 0; cu = 0
    for wb in wins:
        for wz in wb:
            for w in wz:
                if not isinstance(w, list): nzr[c] = cu+cs; cu += 1
                c += 1

    bases = Internal.getBases(CTK.t)
    for nz in nzs: # pour chaque zone selectionnee
        c = 0; nob = 0
        for wb in wins: # par base
            b = bases[nob]
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            noz = 0
            for wz in wb: # par zone
                z = zones[noz]
                dim = Internal.getZoneDim(z)
                for w in wz:
                    if nzr[c] == nz:
                        if dim[0] == 'Structured': # structure
                            C._addBC2Zone(z, nameBC, typeBC, w)
                        elif dim[0] == 'Unstructured' and dim[3] == 'NGON':
                            C._addBC2Zone(z, nameBC, typeBC, faceList=w)
                        else: # BE + BCC
                            zp = T.subzone(z, w, type='faces')
                            zp[0] = C.getZoneName(zp[0])
                            C._addBC2Zone(z, nameBC, typeBC, subzone=zp)
                    c += 1
                noz += 1
            nob += 1
        #print('BC is ', w, 'corresponding ', nob-1, noz-1)

    CTK.TXT.insert('START', 'BCs set to %s.\n'%typeBC)
    CTK.TKTREE.updateApp()
    check()
    CTK.setCursor(0, WIDGETS['setBCWith'])

#==============================================================================
def setSplitFactor(event=None):
    # Seulement pour l'info bulle
    val = WIDGETS['splitFactor'].get()
    VARS[8].set('Split more or less undefined BCs [%.2f deg.]. \nUsefull only for unstructured grids.'%(180.-val*180./100.))

#==============================================================================
def createBCFamily(event=None):
    if CTK.t == []: return
    name = VARS[11].get()
    if name == '':
        CTK.TXT.insert('START', 'FamilyBC name is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        bases = CTK.t[2]
        for b in bases:
            if b[3] == 'CGNSBase_t':
                C._addFamily2Base(b, name, VARS[12].get())
        CTK.TXT.insert('START', 'BC Family '+name+' added to all bases.\n')
    else:
        nob = CTK.Nb[nzs[0]]+1
        noz = CTK.Nz[nzs[0]]
        z = CTK.t[2][nob][2][noz]
        C._addFamily2Base(CTK.t[2][nob], name, VARS[12].get())
        CTK.TXT.insert('START', 'BC Family '+name+' added to base '+CTK.t[2][nob][0]+'.\n')
    CTK.TKTREE.updateApp()

#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkBC  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage boundary conditions.\nCtrl+w to close applet.', btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=4)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkBC')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # - 0 - Type de BC -
    V = TK.StringVar(win); V.set('Mesh'); VARS.append(V)
    # - 1 - Type de cas (2D/3D) -
    V = TK.StringVar(win); V.set('3D'); VARS.append(V)
    # - 2 - tol pour ConnectMatch -
    V = TK.StringVar(win); V.set('1.e-6'); VARS.append(V)
    if 'tkBCMatchTol' in CTK.PREFS: V.set(CTK.PREFS['tkBCMatchTol'])
    # - 3 - ratio pour ConnectNearMatch -
    V = TK.StringVar(win); V.set('2'); VARS.append(V)
    # - 4 - Type de BC pour fillEmptyBCWith -
    V = TK.StringVar(win); V.set('BCFarfield'); VARS.append(V)
    # - 5 - Type de BC pour rmBCOfType -
    V = TK.StringVar(win); V.set('-All BC-'); VARS.append(V)
    # - 6 - Type de BC pour setBCWith -
    V = TK.StringVar(win); V.set('BCWall'); VARS.append(V)
    # - 7 - Edges des zones du calcul
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkBCEdges' in CTK.PREFS: V.set(CTK.PREFS['tkBCEdges'])
    # -8- SplitFactor info bulle
    V = TK.StringVar(win); V.set('Split more or less undefined BCs. \nUsefull only for unstructured grids.'); VARS.append(V)
    # -9- Periodicity? in connectMatch
    V = TK.StringVar(win); V.set('Not periodic'); VARS.append(V)
    # -10- Periodicity field (0;0;0...)
    V = TK.StringVar(win); V.set('0.;0.;0.;0.;0.;0.'); VARS.append(V)
    if 'tkBCMatchPer' in CTK.PREFS: V.set(CTK.PREFS['tkBCMatchPer'])
    # -11- Nom de la FamilyBC (new) -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -12- Type de BC pour createFamilyBC -
    V = TK.StringVar(win); V.set('UserDefined'); VARS.append(V)

    # - View mesh -
    B = TTK.Button(Frame, text="View Mesh", command=viewMesh)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='View mesh.\nTree is NOT modified.')

    # - Edges -
    B = TTK.Checkbutton(Frame, text='Edges', variable=VARS[7])
    BB = CTK.infoBulle(parent=B, text='Show edges of zones of the tree.')
    B.grid(row=0, column=1, sticky=TK.EW)

    # - View type de BC -
    B = TTK.Button(Frame, text="View BC", command=view)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='View specified BC.\nTree is NOT modified.')

    # - Type of BC - ListBox Frame -
    LBFrame = TTK.Frame(Frame)
    LBFrame.grid(row=1, column=1, rowspan=4, sticky=TK.EW)
    LBFrame.rowconfigure(0, weight=1)
    LBFrame.columnconfigure(0, weight=1)
    LBFrame.columnconfigure(1, weight=0)
    SB = TTK.Scrollbar(LBFrame)
    LB  = TTK.Listbox(LBFrame, selectmode=TK.EXTENDED, height=6)
    LB.bind('<Double-1>', view)
    LB.bind('<Enter>', updateBCNameList)
    # LB.bind('<ButtonRelease-1>', view)
    for i, value in enumerate(['-All BC-']+getAllDefinedBC(CTK.t)): LB.insert(i, value)
    SB.config(command=LB.yview)
    LB.config(yscrollcommand=SB.set)
    LB.grid(row=0, column=0, sticky=TK.NSEW)
    SB.grid(row=0, column=1, sticky=TK.NSEW)
    LBFrame.bind('<Enter>', updateFamilyBCNameList1_2)
    WIDGETS['BCLB'] = LB

    # - View undefined BCs -
    B = TTK.Button(Frame, text="View undefined BC",
                   command=displayUndefinedBoundaries)
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    WIDGETS['undefinedBC'] = B
    BB = CTK.infoBulle(parent=B, text='View undefined BC in ALL tree.\nUse this to setBC or fillEmptyBC.')

    # - Slider for splitFactor -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL,
                  command=setSplitFactor, showvalue=0, borderwidth=1, value=0)
    WIDGETS['splitFactor'] = B
    B.grid(row=3, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[8])

    # - rmBCOfType -
    B = TTK.Button(Frame, text="rm BC", command=rmBCOfType)
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Remove the BCs of given type from tree.')

    # - setBCWith -
    B = TTK.Button(Frame, text="setBCWith", command=setBCWith)
    WIDGETS['setBCWith'] = B
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the BC with specified type.\nAdded to pyTree.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[6], *(Internal.KNOWNBCS))
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFamilyBCNameList3)
        F.grid(row=6, column=1, sticky=TK.EW)
        WIDGETS['BCs2'] = B
    else:
        B = TTK.Combobox(F, textvariable=VARS[6],
                         values=Internal.KNOWNBCS, state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFamilyBCNameList3_2)
        F.grid(row=5, column=1, sticky=TK.EW)
        WIDGETS['BCs2'] = B

    # - FillEmptyBCWith -
    B = TTK.Button(Frame, text="FillEmptyBCWith", command=fillEmptyBCWith)
    B.grid(row=6, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fill empty BCs with given type.\nAdded to pyTree.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[4], *(Internal.KNOWNBCS))
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFamilyBCNameList4)
        F.grid(row=6, column=1, sticky=TK.EW)
        WIDGETS['BCs4'] = B
    else:
        B = TTK.Combobox(F, textvariable=VARS[4],
                         values=Internal.KNOWNBCS, state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFamilyBCNameList4_2)
        F.grid(row=6, column=1, sticky=TK.EW)
        WIDGETS['BCs4'] = B

    # - Create BC family -
    B = TTK.Button(Frame, text="NewBCFamily", command=createBCFamily)
    B.grid(row=7, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create a new BC family tag.')
    B = TTK.Entry(Frame, textvariable=VARS[11], background='White', width=15)
    BB = CTK.infoBulle(parent=B, text='BC family name to be created.')
    B.grid(row=7, column=1, sticky=TK.EW)
    B.bind('<Return>', createBCFamily)

    if ttk is None:
        B = TK.OptionMenu(Frame, VARS[12], *(Internal.KNOWNBCS))
        BB = CTK.infoBulle(parent=B, text='BC family type to be created.')
        B.grid(row=8, column=1, sticky=TK.EW)
    else:
        B = TTK.Combobox(Frame, textvariable=VARS[12],
                         values=Internal.KNOWNBCS, state='readonly', width=10)
        BB = CTK.infoBulle(parent=B, text='BC family type to be created.')
        B.grid(row=8, column=1, sticky=TK.EW)

    # - setDegeneratedBC -
    B = TTK.Button(Frame, text="SetDegeneratedBC", command=setDegeneratedBC)
    B.grid(row=8, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Find the degenerated BCs,\nset them as "BCDegenerateLine"\nAdded to pyTree.')

    # - ConnectMatch -
    B = TTK.Button(Frame, text="ConnectMatch", command=connectMatch)
    WIDGETS['connectMatch'] = B
    B.grid(row=9, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Find the matching BCs\nAdded to pyTree.')
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White')
    B.grid(row=9, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Tolerance for matching.')

    # - Periodicity management -
    B = TTK.OptionMenu(Frame, VARS[9], 'Not periodic', 'Translation', 'Rotation (Degree)')
    B.grid(row=10, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[10], background='White')
    B.grid(row=10, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Periodic translation: tx;ty;tz\nPeriodic rotation: cx;cy;cz;ax;ay;az\nangles in degrees.')

    # - ConnectNearMatch -
    B = TTK.Button(Frame, text="ConnectNearMatch", command=connectNearMatch)
    WIDGETS['connectNearMatch'] = B
    B.grid(row=11, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Find the matching BCs\nAdded to pyTree.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White')
    B.grid(row=11, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Nearmatch ratio (2 or 2;1;2).')

#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['BCNoteBook'].add(WIDGETS['frame'], text='tkBC')
    except: pass
    CTK.WIDGETS['BCNoteBook'].select(WIDGETS['frame'])

#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['BCNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
def saveApp():
    CTK.PREFS['tkBCMatchTol'] = VARS[2].get()
    CTK.PREFS['tkBCMatchPer'] = VARS[10].get()
    CTK.PREFS['tkBCEdges'] = VARS[7].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[2].set('1.e-6')
    VARS[10].set('0;0;0;0;0;0')
    VARS[7].set('0')
    CTK.PREFS['tkBCMatchTol'] = VARS[2].get()
    CTK.PREFS['tkBCMatchPer'] = VARS[10].get()
    CTK.PREFS['tkBCEdges'] = VARS[7].get()
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
    (win, menu, file, tools) = CTK.minimal('tkBC'+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
