# - tkBlocks -
"""Block operations in a pyTree."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Post.PyTree as P
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Delete le block selectionne de t
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche
#==============================================================================
def rmBlock():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    deletedZoneNames = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        deletedZoneNames.append(CTK.t[2][nob][0]+Internal.SEP1+CTK.t[2][nob][2][noz][0])

    CTK.saveTree()
    CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)

    CTK.TXT.insert('START', 'Selected zones deleted.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.delete(deletedZoneNames) 
    CPlot.render()

#==============================================================================
# Copie le block selectionne dans t
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche
#==============================================================================
def cpBlock():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'COPY', 3)
    base = Internal.getNodesFromName1(CTK.t, 'COPY')[0]
    gnob = C.getNobOfBase(base, CTK.t)
    newFamilyZoneNames = set()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = Internal.copyRef(CTK.t[2][nob][2][noz])
        z[0] = C.getZoneName(z[0]+'.dup')
        CTK.add(CTK.t, gnob, -1, z)
        # Create new Family Name for Family Zones
        nbdup = z[0].split('.dup')[-1]
        nodes = Internal.getNodesFromType1(z, 'FamilyName_t')
        for f in nodes:
            newFamilyZoneName = Internal.getValue(f)+'.dup'+nbdup
            f[1] = newFamilyZoneName
            if newFamilyZoneName not in newFamilyZoneNames:
                newFamilyZoneNames.add(newFamilyZoneName)
        CTK.TXT.insert('START', CTK.t[2][nob][0]+'/'+CTK.t[2][nob][2][noz][0]
                       +' duplicated.\n')
    
    # Create new Family_t node for each Tag Zone
    for f in newFamilyZoneNames:
        Internal.newFamily(name=f, parent=base)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
# Converti un bloc ou tous les blocs en tetra
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche
#==============================================================================
def convert2Tetra():
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
        try:
            a = C.convertArray2Tetra(CTK.t[2][nob][2][noz])
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail: CTK.TXT.insert('START', 'Zones converted to tetra.\n')
    else:
        Panels.displayErrors(errors, header='Error: convert2Tetra')
        CTK.TXT.insert('START', 'Tetra conversion fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Converti un bloc ou tous les blocs en hexa
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche. 
#==============================================================================
def convert2Hexa():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    list = []
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dim = Internal.getZoneDim(z)
        try:
            if dim[0] == 'Unstructured' and dim[3] == 'TRI':
                a, b = C.convertTri2Quad(z)
                CTK.replace(CTK.t, nob, noz, a)
                CTK.add(CTK.t, nob, -1, b)
            else:
                a = C.convertArray2Hexa(z)
                CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True; errors += [0,str(e)]
    if not fail: CTK.TXT.insert('START', 'Zones converted to hexa.\n')
    else:
        Panels.displayErrors(errors, header='Error: convert2Hexa')
        CTK.TXT.insert('START', 'Hexa conversion fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)   
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Converti un bloc en node
# IN: t, cplot.selectedZones
# OUT: t modifie et affiche. 
#==============================================================================
def convert2Node():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]        
        a = C.convertArray2Node(CTK.t[2][nob][2][noz])
        CTK.replace(CTK.t, nob, noz, a)
        
    CTK.TXT.insert('START', 'Zones converted to node.\n')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
# Exterior faces : recupere les faces externes
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
    p = Internal.getNodesFromName1(CTK.t, 'CONTOURS')
    gnob = C.getNobOfBase(p[0], CTK.t)

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        try:
            ext = P.exteriorFaces(CTK.t[2][nob][2][noz])
            ext = T.splitConnexity(ext)
            for i in ext: CTK.add(CTK.t, gnob, -1, i)
        except TypeError as e: # type d'element non reconnu
            fail = True; errors += [0,str(e)]
        except ValueError: # empty set
            pass
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    if not fail: CTK.TXT.insert('START', 'Exterior faces done.\n')
    else:
        Panels.displayErrors(errors, header='Error: exteriorFaces')
        CTK.TXT.insert('START', 'Exterior faces fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CPlot.render()
    
#=========================================================================
# oneovern: 1 point sur n
# IN: t, cplot.selectedZones
# OUT: t avec zones remplacees et affiche
#=========================================================================
def oneovern():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return 
    args = VARS[0].get()
    args = args.split(';')
    if (len(args) != 3):
        CTK.TXT.insert('START', 'oneovern requires 3 steps.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    ni = int(args[0]); nj = int(args[1]); nk = int(args[2])

    CTK.saveTree()

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        try:
            a = T.oneovern(CTK.t[2][nob][2][noz], (ni,nj,nk))
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail: CTK.TXT.insert('START', 'oneovern done.\n')
    else:
        Panels.displayErrors(errors, header='Error: oneovern')
        CTK.TXT.insert('START', 'oneovern fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# close: ferme les zones
# IN: t, cplot.selectedZones, eps
# OUT: t avec zones remplacees et affiche
#==============================================================================
def close():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    eps = float(VARS[1].get())
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
    try: zones = G.close(zones, eps)
    except Exception as e:
        fail = True; errors += [0,str(e)]
        
    if not fail:
        c = 0
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            a = zones[c]; c += 1
            CTK.replace(CTK.t, nob, noz, a)
            CTK.TXT.insert('START', 'Zones closed.\n')
    else:
        Panels.displayErrors(errors, header='Error: close')
        CTK.TXT.insert('START', 'Close fails at least for one zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE, 
                           text='tkBlock  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='General block operations.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkBlock')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- oneovern -
    V = TK.StringVar(win); V.set('2;2;2'); VARS.append(V)
    if 'tkBlockOneovern' in CTK.PREFS: V.set(CTK.PREFS['tkBlockOneovern'])
    # -1- eps pour close
    V = TK.StringVar(win); V.set('1.e-6'); VARS.append(V)
    if 'tkBlockClose' in CTK.PREFS: V.set(CTK.PREFS['tkBlockClose'])

    # - Rm block from model -
    B = TTK.Button(Frame, text="Rm block", command=rmBlock)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Remove selection from tree.')

    # - Copy block -
    B = TTK.Button(Frame, text="Copy block", command=cpBlock)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Duplicate block in tree.')

    # - Convert a block to tetra -
    B = TTK.Button(Frame, text="convert2Tetra", command=convert2Tetra)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Convert selection to tetra.\nTree is modified.')

    # - Convert a block to Hexa -
    B = TTK.Button(Frame, text="convert2Hexa", command=convert2Hexa)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Convert selection to hexa.\nTree is modified.')

    # - Convert a block to node -
    B = TTK.Button(Frame, text="convert2Node", command=convert2Node)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Convert selection to nodes.\nTree is modified.')

    # - Get exteriorFaces -
    B = TTK.Button(Frame, text="Exterior faces", command=exteriorFaces)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get the exterior faces of selection.\nAdded to tree.')

    # - Close -
    B = TTK.Button(Frame, text="Close", command=close)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Close mesh.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White')
    BB = CTK.infoBulle(parent=B, text='Close accuracy.')
    B.grid(row=3, column=1, sticky=TK.EW)

    # - oneovern -
    B = TTK.Button(Frame, text="Oneovern", command=oneovern)
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract one point over n points.\nTree is modified.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    BB = CTK.infoBulle(parent=B, text='Steps in each direction.')
    B.grid(row=4, column=1, sticky=TK.EW)
    
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['BlockNoteBook'].add(WIDGETS['frame'], text='tkBlock')
    except: pass
    CTK.WIDGETS['BlockNoteBook'].select(WIDGETS['frame'])

#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['BlockNoteBook'].hide(WIDGETS['frame'])
    
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkBlockOneovern'] = VARS[0].get()
    CTK.PREFS['tkBlockClose'] = VARS[1].get()
    CTK.savePrefFile()
    
#==============================================================================
def resetApp():
    VARS[0].set('2;2;2')
    VARS[1].set('1.e-6')
    CTK.PREFS['tkBlockOneovern'] = VARS[0].get()
    CTK.PREFS['tkBlockClose'] = VARS[1].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if (__name__ == "__main__"):
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkBlock '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
