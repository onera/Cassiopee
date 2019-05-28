# - manipulate edges -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Post.PyTree as P
import Transform.PyTree as T
import Generator.PyTree as G
import Intersector.PyTree as XOR

try: range = xrange
except: pass

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def intersection():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if len(nzs) != 2:
        CTK.TXT.insert('START', 'Intersection requires exactely two surfaces.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    p = C.Internal.getNodeFromName1(CTK.t, 'CONTOURS')
    nob = CTK.Nb[0]+1
    noz = CTK.Nz[0]
    a1 = CTK.t[2][nob][2][noz]
    nob = CTK.Nb[1]+1
    noz = CTK.Nz[1]
    a2 = CTK.t[2][nob][2][noz]
    try:
        b = XOR.intersection(a1, a2)
        nob = C.getNobOfBase(p, CTK.t)
        CTK.add(CTK.t, nob, -1, b)
        CTK.TXT.insert('START', 'Intersection edge computed.\n')
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
        return
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: intersection')
        CTK.TXT.insert('START', 'Intersection failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
#==============================================================================
# Exterior faces: recupere les faces externes
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche. Les faces externes sont stockes dans une
# base CONTOURS. Ce sont des zones non-structurees.
#==============================================================================
def exteriorFaces():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    p = C.Internal.getNodeFromName1(CTK.t, 'CONTOURS')

    fail = False; exts = []; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Unstructured':
            z = G.close(z)
            try:
                ext = P.exteriorFaces(z)
                ext = G.close(ext)
                ext = T.splitConnexity(ext)
                exts += ext
            except Exception as e: 
                fail = True; errors += [0,str(e)]
        else:
            ext = P.exteriorFacesStructured(z)
            exts += ext

    if fail:
        Panels.displayErrors(errors, header='Error: exteriorFaces') 
        
    if exts == []:
        CTK.TXT.insert('START', 'External edges set is empty.\n')
    else:
        nob = C.getNobOfBase(p, CTK.t)
        for i in exts: CTK.add(CTK.t, nob, -1, i)
        if not fail:
            CTK.TXT.insert('START', 'External edges extracted.\n')
        else:
            print('Error: externalEdges: %s.'%str(e))
            CTK.TXT.insert('START', 'External edges fails for at least one zone.\n')
            CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
# silhouette: recupere les edges de silhouette
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche. Les edges sont stockes dans une
# base CONTOURS. Ce sont des zones BARS.
#==============================================================================
def silhouette():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    p = C.Internal.getNodeFromName1(CTK.t, 'CONTOURS')

    cam = CPlot.getState('posCam')
    eye = CPlot.getState('posEye')
    vector = [cam[0]-eye[0], cam[1]-eye[1], cam[2]-eye[2]]

    fail = False; exts = []; errors = []    
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = G.close(z, 1.e-6)
        try:
            ext = P.silhouette(z, vector)
            exts += ext
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if fail: 
        Panels.displayErrors(errors, header='Error: silhouette')
    if exts == []:
        CTK.TXT.insert('START', 'Silhouette set is empty.\n')
    else:
        nob = C.getNobOfBase(p, CTK.t)
        for i in exts: CTK.add(CTK.t, nob, -1, i)
        if not fail:
            CTK.TXT.insert('START', 'Silhouette edges extracted.\n')
        else:
            CTK.TXT.insert('START', 'Silhouette fails for at least one zone.\n')
            CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
# sharpEdges
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche
#==============================================================================
def sharpEdges():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    alphaRef = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(alphaRef) != 1:
        CTK.TXT.insert('START', 'Angle is incorrect.\n');
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    alphaRef = alphaRef[0]
    
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    
    CTK.saveTree()
    sharps = []
    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    p = C.Internal.getNodeFromName1(CTK.t, 'CONTOURS')
    
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = G.close(z)
        try: ext = P.sharpEdges(z, alphaRef)
        except Exception as e:
            ext = []; fail = True; errors += [0,str(e)]
        sharps += ext

    if fail: 
        Panels.displayErrors(errors, header='Error: sharpEdges')

    if sharps == []:
         CTK.TXT.insert('START', 'Sharp edges set is empty.\n');
    else:
        nob = C.getNobOfBase(p, CTK.t)
        for i in sharps: CTK.add(CTK.t, nob, -1, i)
        if not fail:
            CTK.TXT.insert('START', 'Sharp edges extracted.\n')
        else:
            CTK.TXT.insert('START', 'Sharp edges fails for at least one zone.\n')
            CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
# splitTBranches
#==============================================================================
def splitTBranches():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = CTK.varsFromWidget(VARS[1].get(), type=1)
    if (len(tol) != 1):
        CTK.TXT.insert('START', 'Split tolerance is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = tol[0]

    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    
    # Type de la premiere zone
    nob = CTK.Nb[0]+1
    noz = CTK.Nz[0]
    z = CTK.t[2][nob][2][noz]
    dim = Internal.getZoneDim(z)
    structured = 0
    if (dim[0] == 'Structured'): structured = 1

    CTK.saveTree()

    fail = False; errors = []
    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = C.convertArray2Tetra(z)
        zones.append(z)

    a = T.join(zones)
    a = G.close(a)
    try:
        splits = T.splitTBranches(a, tol)
        if (structured == 1): splits = C.convertBAR2Struct(splits)
        n = len(nzs); ns = len(splits)

        for c in range(n):
            nob = CTK.Nb[nzs[c]]+1
            noz = CTK.Nz[nzs[c]]
            if c < ns:
                CTK.replace(CTK.t, nob, noz, splits[c])
            else:
                baseName = CTK.t[2][nob][0]
                zoneName = CTK.t[2][nob][2][noz][0]
                i = CPlot.getCPlotNumber(CTK.t, baseName, zoneName)
                Internal._rmNode(CTK.t, CTK.t[2][nob][2][noz])
                CPlot.delete([i])
        for i in splits[n:]: CTK.add(CTK.t, nob, -1, i)

    except Exception as e:
        fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'splitTBranches done.\n')
    else:
        Panels.displayErrors(errors, header='Error: splitTBranches')
        CTK.TXT.insert('START', 'Split T branch fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def convertBAR2Struct():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
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
            zp = C.convertBAR2Struct(z)
            CTK.replace(CTK.t, nob, noz, zp)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'BARS converted to STRUCT.\n')
    else:
        Panels.displayErrors(errors, header='Error: convertBAR2Struct')
        CTK.TXT.insert('START', 'SRUCT conversion fails for at least one zone.\n')
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
                           text='tkExtractEdges', font=CTK.FRAMEFONT, 
                           takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Extract edges.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkExtractEdges')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Edge split angle -
    V = TK.StringVar(win); V.set('30.'); VARS.append(V)
    if 'tkExtractEdgesSplitAngle' in CTK.PREFS:
        V.set(CTK.PREFS['tkExtractEdgesSplitAngle'])
    # -1- Tolerance for T-branches splitting -
    V = TK.StringVar(win); V.set('1.e-8'); VARS.append(V)
    
    # - External edges -
    B = TTK.Button(Frame, text="External edges", command=exteriorFaces)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract the exterior edges of a surface.\nAdded to tree.')

    # - Surface intersection -
    B = TTK.Button(Frame, text="Intersection", command=intersection)
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get intersection of exactly two surfaces.')

    # - sharpEdges -
    B = TTK.Button(Frame, text="sharpEdges", command=sharpEdges)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract sharp edges of a surface.\nTree is modified.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=10)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Threshold angle.')

    # - silhouette -
    #B = TTK.Button(Frame, text="Silhouette", command=silhouette)
    #B.grid(row=2, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Extract silhouette edges of a surface.\nTree is modified.')

    # - splitTBranches -
    B = TTK.Button(Frame, text="splitTBranches", command=splitTBranches)
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split BAR edges at T-branch points.\nTree is modified.')
    #B = TK.Entry(Frame, textvariable=VARS[1], background='White', width=10)
    #B.grid(row=2, column=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Tolerance.')

    # - convertBAR2Struct -
    B = TTK.Button(Frame, text="BAR2Struct", command=convertBAR2Struct)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Convert BARs to i-arrays.')
    
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
    CTK.PREFS['tkExtractEdgesSplitAngle'] = VARS[0].get()
    CTK.savePrefFile()
    
#==============================================================================
def resetApp():
    VARS[0].set('30.')
    CTK.PREFS['tkExtractEdgesSplitAngle'] = VARS[0].get()
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
    (win, menu, file, tools) = CTK.minimal('tkExtractEdges '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
