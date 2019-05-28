# - slice -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Converter.Internal as Internal
import Post.PyTree as P
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Generator.PyTree as G

# local widgets list
WIDGETS = {}; VARS = []

VALUE = 0.
DELTA = 0.1
    
#==============================================================================
def fit():
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []: bb = G.bbox(CTK.t)
    else:
        sel = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            sel.append(z)
        bb = G.bbox(sel)
        
    xmin = bb[0]; ymin = bb[1]; zmin = bb[2]
    xmax = bb[3]; ymax = bb[4]; zmax = bb[5]
    
    plane = VARS[0].get()

    if plane == 'Mesh':
        delta = max(xmax-xmin, ymax-ymin)
        delta = max(delta, zmax-zmin)
        delta = delta / 50.
        VARS[4].set(str(delta))
        x = 0.5*(xmax+xmin)
        VARS[1].set(str(x))
    elif plane == 'X':
        delta = (xmax-xmin) / 50.
        VARS[4].set(str(delta))
        x = 0.5*(xmax+xmin)
        VARS[1].set(str(x))
    elif plane == 'Y':
        delta = (ymax-ymin) / 50.
        VARS[4].set(str(delta))
        x = 0.5*(ymax+ymin)
        VARS[1].set(str(x))
    elif plane == 'Z':
        delta = (zmax-zmin) / 50.
        VARS[4].set(str(delta))
        x = 0.5*(zmax+zmin)
        VARS[1].set(str(x))
    
#==============================================================================
def movePlus():
    pos = float(VARS[1].get())
    delta = float(VARS[4].get())
    pos = pos + delta
    VARS[1].set(str(pos))
    view()

#==============================================================================
def moveMoins():
    pos = float(VARS[1].get())
    delta = float(VARS[4].get())
    pos = pos - delta
    VARS[1].set(str(pos))
    view()
    
#==============================================================================
def unselect(event=None):
    CPlot.unselectAllZones()

#==============================================================================
def view(event=None):
    if CTK.t == []: return
    pos = float(VARS[1].get())
    global VALUE; VALUE = pos
    delta = float(VARS[4].get())
    global DELTA; DELTA = delta
    plane = VARS[0].get()
    order = int(VARS[3].get())
    eps = float(VARS[2].get())
    algo = VARS[5].get()

    nzs = CPlot.getSelectedZones()
    if nzs != []:
        point = CPlot.getActivePoint()
        if len(point) == 3:
            if plane == 'X': pos = point[0]
            elif plane == 'Y': pos = point[1]
            elif plane == 'Z': pos = point[2]
            VARS[1].set(str(pos)); VALUE = pos
        
    if plane == 'Mesh': CTK.display(CTK.t); return
    try:
        if CTK.__MAINTREE__ == 1:
            CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()
        active = []
        tp = Internal.appendBaseName2ZoneName(CTK.t, updateRef=False, 
                                              separator=Internal.SEP1)
        for z in CTK.__MAINACTIVEZONES__: active.append(tp[2][CTK.Nb[z]+1][2][CTK.Nz[z]])

        temp = C.newPyTree(['Base']); temp[2][1][2] += active
        if plane == 'X' and algo == 'Slice1':
            p = P.isoSurfMC(active, 'CoordinateX', pos)
        elif plane == 'Y' and algo == 'Slice1':
            p = P.isoSurfMC(active, 'CoordinateY', pos)
        elif plane == 'Z' and algo == 'Slice1':
            p = P.isoSurfMC(active, 'CoordinateZ', pos)
        elif plane == 'X' and algo == 'Slice2':
            p = P.extractPlane(active, (1,0,0,-pos), order=order, tol=eps)
        elif plane == 'Y' and algo == 'Slice2':
            p = P.extractPlane(active, (0,1,0,-pos), order=order, tol=eps)
        elif plane == 'Z' and algo == 'Slice2':
            p = P.extractPlane(active, (0,0,1,-pos), order=order, tol=eps)
        elif plane == 'X' and algo == 'Select+':
            p = P.selectCells(temp, '{CoordinateX}>='+str(VALUE))
        elif plane == 'Y' and algo == 'Select+':
            p = P.selectCells(temp, '{CoordinateY}>='+str(VALUE))
        elif plane == 'Z' and algo == 'Select+':
            p = P.selectCells(temp, '{CoordinateZ}>='+str(VALUE))
        elif plane == 'X' and algo == 'Select-':
            p = P.selectCells(temp, '{CoordinateX}<='+str(VALUE))
        elif plane == 'Y' and algo == 'Select-':
            p = P.selectCells(temp, '{CoordinateY}<='+str(VALUE))
        elif plane == 'Z' and algo == 'Select-':
            p = P.selectCells(temp, '{CoordinateZ}<='+str(VALUE))
        elif plane == 'X' and algo == 'Select=':
            p = P.selectCells(temp, '({CoordinateX}>='+str(VALUE-DELTA)+') & ({CoordinateX}<='+str(VALUE+DELTA)+')')
        elif plane == 'Y' and algo == 'Select=':
            p = P.selectCells(temp, '({CoordinateY}>='+str(VALUE-DELTA)+') & ({CoordinateY}<='+str(VALUE+DELTA)+')')
        elif plane == 'Z' and algo == 'Select=':
            p = P.selectCells(temp, '({CoordinateZ}>='+str(VALUE-DELTA)+') & ({CoordinateZ}<='+str(VALUE+DELTA)+')')
        CTK.dt = C.newPyTree(['Base'])
        if algo == 'Slice1': CTK.dt[2][1][2] += p
        elif algo == 'Slice2': CTK.dt[2][1][2] += [p]
        else: CTK.dt[2][1][2] += p[2][1][2]
        CTK.display(CTK.dt, mainTree=CTK.SLICE)
        if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()
    except ValueError:
        CTK.TXT.insert('START', 'Intersection is empty.\n'); return
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: slice')
        CTK.TXT.insert('START', 'Slice failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return    

#==============================================================================
def extract(event=None):
    if CTK.t == []: return
    pos = float(VARS[1].get())
    global VALUE; VALUE = pos
    delta = float(VARS[4].get())
    global DELTA; DELTA = delta
    plane = VARS[0].get()
    order = int(VARS[3].get())
    eps = float(VARS[2].get())
    algo = VARS[5].get()

    nzs = CPlot.getSelectedZones()
    if nzs != []:
        point = CPlot.getActivePoint()
        if plane == 'X': pos = point[0]
        elif plane == 'Y': pos = point[1]
        elif plane == 'Z': pos = point[2]
        VARS[1].set(str(pos)); VALUE = pos
         
    if plane == 'Mesh': return
    try:
        CTK.saveTree()
        if CTK.__MAINTREE__ == 1:
            CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()
        active = []
        zones = Internal.getZones(CTK.t)
        for z in CTK.__MAINACTIVEZONES__: active.append(zones[z])
        temp = C.newPyTree(['Base']); temp[2][1][2] += active
        if plane == 'X' and algo == 'Slice1':
            p = P.isoSurfMC(active, 'CoordinateX', pos)
        elif plane == 'Y' and algo == 'Slice1':
            p = P.isoSurfMC(active, 'CoordinateY', pos)
        elif plane == 'Z' and algo == 'Slice1':
            p = P.isoSurfMC(active, 'CoordinateZ', pos)
        elif plane == 'X' and algo == 'Slice2':
            p = P.extractPlane(active, (1,0,0,-pos), order=order, tol=eps)
        elif plane == 'Y' and algo == 'Slice2':
            p = P.extractPlane(active, (0,1,0,-pos), order=order, tol=eps)
        elif plane == 'Z' and algo == 'Slice2':
            p = P.extractPlane(active, (0,0,1,-pos), order=order, tol=eps)
        elif plane == 'X' and algo == 'Select+':
            p = P.selectCells(temp, '{CoordinateX}>='+str(VALUE))
        elif plane == 'Y' and algo == 'Select+':
            p = P.selectCells(temp, '{CoordinateY}>='+str(VALUE))
        elif plane == 'Z' and algo == 'Select+':
            p = P.selectCells(temp, '{CoordinateZ}>='+str(VALUE))
        elif plane == 'X' and algo == 'Select-':
            p = P.selectCells(temp, '{CoordinateX}<='+str(VALUE))
        elif plane == 'Y' and algo == 'Select-':
            p = P.selectCells(temp, '{CoordinateY}<='+str(VALUE))
        elif plane == 'Z' and algo == 'Select-':
            p = P.selectCells(temp, '{CoordinateZ}<='+str(VALUE))
        elif plane == 'X' and algo == 'Select=':
            p = P.selectCells(temp, '({CoordinateX}>='+str(VALUE-DELTA)+') & ({CoordinateX}<='+str(VALUE+DELTA)+')')
        elif plane == 'Y' and algo == 'Select=':
            p = P.selectCells(temp, '({CoordinateY}>='+str(VALUE-DELTA)+') & ({CoordinateY}<='+str(VALUE+DELTA)+')')
        elif plane == 'Z' and algo == 'Select=':
            p = P.selectCells(temp, '({CoordinateZ}>='+str(VALUE-DELTA)+') & ({CoordinateZ}<='+str(VALUE+DELTA)+')')
        CTK.t = C.addBase2PyTree(CTK.t, 'EXTRACT', 2)
        base = Internal.getNodeFromName1(CTK.t, 'EXTRACT')
        if algo == 'Slice1':
            for i in p: i[0] = C.getZoneName(i[0])
            base[2] += p
        elif algo == 'Slice2':
            p[0] = C.getZoneName(p[0])
            base[2] += [p]
        else:
            p = C.deleteEmptyZones(p)
            for i in p[2][1][2]: i[0] = C.getZoneName(i[0])
            base[2] += p[2][1][2]
        #C._fillMissingVariables(CTK.t)
        CTK.TXT.insert('START', 'Slice extracted.\n')
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CTK.display(CTK.t)
        if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()
    except ValueError:
        CTK.TXT.insert('START', 'Intersection is empty.\n'); return
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: slice')
        CTK.TXT.insert('START', 'Slice failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkSlice', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Visualize/extract slices.\nCtrl+c to close applet.', btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=4)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkSlice')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Slice dir
    V = TK.StringVar(win); V.set('X'); VARS.append(V)
    # -1- position
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    #V.trace_add("write", unselect)
    V.trace("w", lambda name, index, mode, V=V: unselect(V))
    # -2- epsilon for 2D slices
    V = TK.StringVar(win); V.set('1.e-6'); VARS.append(V)
    # -3- Order
    V = TK.StringVar(win); V.set('2'); VARS.append(V)
    # -4- Delta pour le move
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)
    if 'tkSliceStep' in CTK.PREFS: V.set(CTK.PREFS['tkSliceStep'])
    # -5- slice algorithm
    V = TK.StringVar(win); V.set('Slice1'); VARS.append(V)

    # - Settings -
    #B = TK.Entry(Frame, textvariable=VARS[2], background='White', width=3)
    #B.grid(row=0, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Tolerance for surface slice.')
    #B = TK.Entry(Frame, textvariable=VARS[3], background='White', width=3)
    #B.grid(row=0, column=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Order for interpolation.')
    #B = TK.Button(Frame, text="Fit", command=fit)
    #B.grid(row=0, column=2, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Fit slice to selection.')

    # Move
    B = TTK.Button(Frame, text="+", command=movePlus)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Move slice +.')
    B = TTK.Button(Frame, text="-", command=moveMoins)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Move slice -.')
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=3)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Move step.')

    # Slice algorithm
    B = TTK.OptionMenu(Frame, VARS[5], 'Slice1', 'Slice2', 'Select+', 'Select-', 'Select=')
    BB = CTK.infoBulle(parent=B, text='Type of slice.')
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    
    # - Position / type -
    B = TTK.OptionMenu(Frame, VARS[0], 'Mesh', 'X', 'Y', 'Z')
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Slice direction.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=3)
    B.grid(row=1, column=2, sticky=TK.EW)
    B.bind('<Return>', unselect)
    BB = CTK.infoBulle(parent=B, text='Plane position.\nTaken from selection or set it here with no selection.')

    # - Extract / view -
    B = TTK.Button(Frame, text="View", command=view)
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='View a slice.')
    B = TTK.Button(Frame, text="Extract", command=extract)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract a slice to pyTree.')
    B = TTK.Button(Frame, text="Fit", command=fit)
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fit slice to selection.')
    
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
    CTK.PREFS['tkSliceStep'] = VARS[4].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[4].set('0.1')
    CTK.PREFS['tkSliceStep'] = VARS[4].get()
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
    (win, menu, file, tools) = CTK.minimal('tkPersonal '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
