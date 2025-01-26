# - projection apps -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Transform.PyTree as T
import Converter.Internal as Internal
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setSurface():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    selected = ''
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        selected += CTK.t[2][nob][0]+'/'+z[0]+';'
    selected = selected[0:-1]
    VARS[0].set(selected)

#==============================================================================
def orthoProject():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # surfaces
    name = VARS[0].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): surfaces.append(z)
    if surfaces == []:
        CTK.TXT.insert('START', 'Projection surface is empty.\n')
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
        try:
            a = T.projectOrtho(z, surfaces)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Zones projected.\n')
    else:
        Panels.displayErrors(errors, header='Error: projectOrtho')
        CTK.TXT.insert('START', 'Projection fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def dirProject():
    if CTK.t == []: return
    if (CTK.__MAINTREE__ <= 0):
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    # surfaces
    name = VARS[0].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if (bases != []):
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): surfaces.append(z)
    if (surfaces == []):
        CTK.TXT.insert('START', 'Projection surface is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    eye = CPlot.getState('posEye')
    cam = CPlot.getState('posCam')
    dir = (eye[0]-cam[0], eye[1]-cam[1], eye[2]-cam[2])
    CTK.saveTree()
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            a = T.projectDir(z, surfaces, dir)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if (fail == False):
        CTK.TXT.insert('START', 'Zones projected.\n')
    else:
        Panels.displayErrors(errors, header='Error: projectDir')
        CTK.TXT.insert('START', 'Projection fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkProjection  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Project meshes on surface meshes.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=2)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkProjection')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Projection surface -
    V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Surface -
    B = TTK.Button(Frame, text="Proj. surf", command=setSurface,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set projection surfaces.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Projection surfaces.')

    # - Orthogonal projection -
    B = TTK.Button(Frame, text="Ortho project", command=orthoProject)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Project selection orthogonally on a surface.')

    # - Dir projection -
    B = TTK.Button(Frame, text="Dir project", command=dirProject)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Project selection following eye direction on a surface.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkProjection')
    except: pass
    CTK.WIDGETS['SurfNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['SurfNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

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
    (win, menu, file, tools) = CTK.minimal('tkProjection '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
