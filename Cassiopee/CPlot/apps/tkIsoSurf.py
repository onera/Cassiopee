# - tkIsoSurf -
"""Compute isosurfaces."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Post.PyTree as P
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def updateVarNameList(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        zvars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        zvars = C.getVarNames(CTK.t[2][nob][2][noz])
    m = WIDGETS['field'].children['menu']
    m.delete(0, TK.END)
    if len(zvars) == 0: return
    for i in zvars[0]:
        m.add_command(label=i, command=lambda v=VARS[0],l=i:v.set(l))

def updateVarNameList2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        zvars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        zvars = C.getVarNames(CTK.t[2][nob][2][noz])

    if len(zvars) == 0: return
    if 'field' in WIDGETS:
        WIDGETS['field']['values'] = zvars[0]

#==============================================================================
def extractIsoSurf(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    field = VARS[0].get()
    value = VARS[1].get()

    try: value = float(value)
    except: value = 1.

    nzs = CPlot.getSelectedZones()
    CTK.saveTree()

    if nzs == []:
        z = Internal.getZones(CTK.t)
    else:
        z = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z.append(CTK.t[2][nob][2][noz])

    isos = []
    CTK.setCursor(2, WIDGETS['frame'])
    try:
        iso = P.isoSurfMC(z, field, value)
        isos += iso
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: isoSurf')
    if isos == []:
        CTK.TXT.insert('START', 'isoSurf failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    else:
        CTK.TXT.insert('START', 'isoSurf of '+field+'='
                       +str(value)+' computed.\n')
    CTK.setCursor(0, WIDGETS['frame'])
    for i in isos: i[0] = C.getZoneName(i[0]) # unique name
    CTK.t = C.addBase2PyTree(CTK.t, 'SURFACES', 2)
    base = Internal.getNodeFromName1(CTK.t, 'SURFACES')
    nob = C.getNobOfBase(base, CTK.t)
    for i in isos: CTK.add(CTK.t, nob, -1, i)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def getValueFromMouse():
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    var = VARS[0].get()
    if var == 'CoordinateX':
        point = CPlot.getActivePoint()
        val = point[0]
    elif var == 'CoordinateY':
        point = CPlot.getActivePoint()
        val = point[1]
    elif var == 'CoordinateZ':
        point = CPlot.getActivePoint()
        val = point[2]
    else:
        val = None; c = 0
        for i in WIDGETS['field']['values']:
            if var == i: break
            c += 1
        c = c-3 # a cause des coord
        values = CPlot.getActivePointF()
        if values != []: val = values[c]
    if val is not None: VARS[1].set(str(val))

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkIsoSurf  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Compute iso-surfaces.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkIsoSurf')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Field name -
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -1- value -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkIsoSurfValue' in CTK.PREFS: 
        V.set(CTK.PREFS['tkIsoSurfValue'])

    # - field name -
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[0], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList)
        F.grid(row=0, column=0, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Extracted field.')
        WIDGETS['field'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[0], 
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=0, column=0, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Extracted field.')
        WIDGETS['field'] = B

    # - Value -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=7)
    BB = CTK.infoBulle(parent=B, text='Extracted field value.')
    B.bind('<Return>', extractIsoSurf)
    B.grid(row=0, column=1, sticky=TK.EW)

    # - Get value from mouse -
    B = TTK.Button(Frame, image=iconics.PHOTO[8],
                   command=getValueFromMouse, padx=0)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get value from mouse.')

    # - Extract -
    B = TTK.Button(Frame, text="Extract isosurf", command=extractIsoSurf)
    B.grid(row=1, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract isosurf to SURFACES.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['PostNoteBook'].add(WIDGETS['frame'], text='tkIsoSurf')
    except: pass
    CTK.WIDGETS['PostNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['PostNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkIsoSurfValue'] = VARS[1].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[1].set('1.')
    CTK.PREFS['tkIsoSurfValue'] = VARS[1].get()
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
    (win, menu, file, tools) = CTK.minimal('tkIsoSurf '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
