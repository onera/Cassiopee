# - Reorder des blocs d'un pyTree -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Transform.PyTree as T
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Fait un reorder de t ou d'une partie de t
# IN: t, cplot.selectedZone
# OUT: t modifie et affiche
#==============================================================================
def reorder():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    ii = VARS[0].get()
    if ii == 'I -> I': i1 = 1
    elif ii == 'I -> -I': i1 = -1
    elif ii == 'I -> J': i1 = 2
    elif ii == 'I -> -J': i1 = -2
    elif ii == 'I -> K': i1 = 3
    elif ii == 'I -> -K': i1 = -3
    jj = VARS[1].get()
    if jj == 'J -> I': j1 = 1
    elif jj == 'J -> -I': j1 = -1
    elif jj == 'J -> J': j1 = 2
    elif jj == 'J -> -J': j1 = -2
    elif jj == 'J -> K': j1 = 3
    elif jj == 'J -> -K': j1 = -3
    kk = VARS[2].get()
    if kk == 'K -> I': k1 = 1
    elif kk == 'K -> -I': k1 = -1
    elif kk == 'K -> J': k1 = 2
    elif kk == 'K -> -J': k1 = -2
    elif kk == 'K -> K': k1 = 3
    elif kk == 'K -> -K': k1 = -3
    if abs(i1)+abs(j1)+abs(k1) != 6:
        CTK.TXT.insert('START', 'Reordering settings is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dim = Internal.getZoneDim(z)
        try:
            if dim[0] == 'Unstructured':
                T._reorder(z, (-1,), topTree=CTK.t)
            else:
                #ni = dim[1]; nj = dim[2]; nk = dim[3]
                #if (nk == 1): a = T.reorder(z, (-1,2,3), topTree=CTK.t)
                #elif (nj == 1): a = T.reorder(z, (-1,2,3), topTree=CTK.t)
                #elif (ni == 1): a = T.reorder(z, (1,2,-3), topTree=CTK.t)
                #else: a = T.reorder(z, (-1,2,3), topTree=CTK.t)
                T._reorder(z, (i1,j1,k1), topTree=CTK.t)
            CTK.replace(CTK.t, nob, noz, z)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Zones reordered.\n')
    else:
        Panels.displayErrors(errors, header='Error: reorder')
        CTK.TXT.insert('START', 'Reorder fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Fait un reorderAll de t
# IN: t
# OUT: t modifie et affiche
#==============================================================================
def reorderAll():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    ii = VARS[0].get()
    if ii == 'I -> I': i1 = 1
    elif ii == 'I -> -I': i1 = -1
    try:
        CTK.saveTree()
        CTK.t = T.reorderAll(CTK.t, i1)
        CTK.TKTREE.updateApp()
        CTK.display(CTK.t)
        CTK.TXT.insert('START', 'All blocks reordered.\n')
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: reorderAll')
        CTK.TXT.insert('START', 'Reorder all fails.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# Fait un makeDirect sur la selection
# IN: t, cplot.selectedZone
# OUT: t modifie et affiche
#==============================================================================
def makeDirect():
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
        try:
            a = T.makeDirect(z)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Zones made direct.\n')
    else:
        Panels.displayErrors(errors, header='Error: makeDirect')
        CTK.TXT.insert('START', 'MakeDirect fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkReorder  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Reorder blocks.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkReorder')
    WIDGETS['frameMenu'] = FrameMenu

    # -0- Transformation de i -
    V = TK.StringVar(win); V.set('I -> I'); VARS.append(V)
    # -1- Transformation de j -
    V = TK.StringVar(win); V.set('J -> J'); VARS.append(V)
    # -2- Transformation de k -
    V = TK.StringVar(win); V.set('K -> K'); VARS.append(V)

    # - Index switch for structured grids 
    B = TTK.OptionMenu(Frame, VARS[0], 'I -> I', 'I -> -I', 'I -> J', 'I -> -J', 'I -> K', 'I -> -K')
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Transformation of i index.')
    B = TTK.OptionMenu(Frame, VARS[1], 'J -> I', 'J -> -I', 'J -> J', 'J -> -J', 'J -> K', 'J -> -K')
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Transformation of j index.')
    B = TTK.OptionMenu(Frame, VARS[2], 'K -> I', 'K -> -I', 'K -> J', 'K -> -J', 'K -> K', 'K -> -K')
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Transformation of k index.')

    # - Bouton reorder -
    B = TTK.Button(Frame, text="Reorder", command=reorder)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reorder the numerotation of block.')

    # - Bouton reorderAll -
    B = TTK.Button(Frame, text="ReorderAll", command=reorderAll)
    BB = CTK.infoBulle(parent=B, text='Reorder all blocks even with overlaps.')
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)

    # - Bouton makeDirect -
    B = TTK.Button(Frame, text="makeDirect", command=makeDirect)
    BB = CTK.infoBulle(parent=B, text='Make structured blocks direct.')
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)

#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['BlockNoteBook'].add(WIDGETS['frame'], text='tkReorder')
    except: pass
    CTK.WIDGETS['BlockNoteBook'].select(WIDGETS['frame'])

#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['BlockNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
def updateApp(): return

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t, mode=1)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkReorder'+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
