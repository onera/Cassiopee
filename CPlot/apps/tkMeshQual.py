# - mesh quality -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Generator.PyTree as G
import Converter.Internal as Internal
import Post.PyTree as P

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Filters
#==============================================================================
def F1(v):
    if v <= 0: return 1
    else: return 0
    
#==============================================================================
# La variable var existe t'elle dans la premiere zone de l'arbre?
# Retourne 0: non
# Retourne 1: oui, en noeuds
# Retourne 2: oui, en centres
#==============================================================================
def findVar(var):
    v = C.getVarNames(CTK.t)
    if len(v) > 0: vars = v[0]
    else: vars = []
    for i in vars:
        if (var == i): return 1
        if ('centers:'+var == i): return 2
    return 0

#==============================================================================
def computeQual():
    if CTK.t == []: return
    CTK.saveTree()
    qtype = VARS[0].get()
    if qtype == 'Volume map':
        CTK.t = G.getVolumeMap(CTK.t)
        CTK.TXT.insert('START', 'Volume map computed.\n')
    elif qtype == 'Orthogonality map':
        CTK.t = G.getOrthogonalityMap(CTK.t)
        CTK.TXT.insert('START', 'Orthogonality map computed.\n')
    elif qtype == 'Regularity map':
        CTK.t = G.getRegularityMap(CTK.t)
        CTK.TXT.insert('START', 'Regularity map computed.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def viewQual():
    if CTK.t == []: return
    qtype = VARS[1].get()
    if qtype == 'Neg. volume cells':
        res = findVar('vol')
        if res == 0: CTK.t = G.getVolumeMap(CTK.t)
        if CTK.__MAINTREE__ == 1:
            CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()
        active = []
        zones = Internal.getZones(CTK.t)
        for z in CTK.__MAINACTIVEZONES__: active.append(CTK.t[2][CTK.Nb[z]+1][2][CTK.Nz[z]])

        temp = C.newPyTree(['Base']); temp[2][1][2] += active
        Z = C.initVars(temp, 'centers:__tag__', F1, ['centers:vol'])
        Z = P.selectCells2(Z, 'centers:__tag__')
        if Z is not None:
            CTK.TXT.insert('START', 'Viewing '+qtype+'.\n')
            CTK.dt = Z
            CTK.display(CTK.dt, mainTree=CTK.MESHQUAL)
    elif qtype == 'Mesh':
        CTK.display(CTK.t)
                
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkMeshQual', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Analyse mesh quality.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkMeshQual')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Qual to compute -
    V = TK.StringVar(win); V.set('Volume map'); VARS.append(V)
    # -0- Filter for view -
    V = TK.StringVar(win); V.set('Mesh'); VARS.append(V)
    
    # - Compute -
    B = TTK.Button(Frame, text="Compute", command=computeQual)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Compute this field. Tree is modified.')
    B = TTK.OptionMenu(Frame, VARS[0], 'Volume map', 'Orthogonality map', 'Regularity map')
    B.grid(row=0, column=1, sticky=TK.EW)

    # - View -
    B = TTK.Button(Frame, text="View", command=viewQual)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='View this filter.')
    B = TTK.OptionMenu(Frame, VARS[1], 'Mesh', 'Neg. volume cells')
    B.grid(row=1, column=1, sticky=TK.EW)
    
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
    (win, menu, file, tools) = CTK.minimal('tkMeshQual '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
