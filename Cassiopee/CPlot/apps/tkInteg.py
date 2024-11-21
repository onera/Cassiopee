# - tkInteg -
"""Perform field integration."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Post.PyTree as P

WIDGETS = {}; VARS = []

#==============================================================================
def updateVarNameList__(no):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        vars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz])
    m = WIDGETS['variable'+str(no)].children['menu']
    m.delete(0, TK.END)
    allvars = []
    if len(vars) > 0:
        for v in vars[0]: allvars.append(v)
    for i in allvars:
        m.add_command(label=i, command=lambda v=VARS[1+no], l=i:v.set(l))

def updateVarNameList2__(no):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        vars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz])

    allvars = []
    if len(vars) > 0:
        for v in vars[0]: allvars.append(v)
    if 'variable'+str(no) in WIDGETS:
        WIDGETS['variable'+str(no)]['values'] = allvars

#==============================================================================
def updateVarNameList1(event=None):
    updateVarNameList__(1)
def updateVarNameList2(event=None):
    updateVarNameList__(2)
def updateVarNameList3(event=None):
    updateVarNameList__(3)
def updateVarNameList1_2(event=None):
    updateVarNameList2__(1)
def updateVarNameList2_2(event=None):
    updateVarNameList2__(2)
def updateVarNameList3_2(event=None):
    updateVarNameList2__(3)

#==============================================================================
def compute():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    type = VARS[0].get(); var1 = VARS[2].get()

    if type == 'INT((v1,v2,v3).ndS)' or type == 'INT((v1,v2,v3)^OMdS)':
        var2 = VARS[3].get(); var3 = VARS[4].get()
        vector= [var1,var2,var3]
    if type == 'INT((v1,v2,v3)^OMdS)' or type == 'INT(v1n^OMdS)':
        center = CTK.varsFromWidget(VARS[1].get(), type=1)
        if len(center) != 3:
            CTK.TXT.insert('START', 'Center for moment integration is incorrect.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    res1 = 0.; res2 = 0.; res3 = 0.
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        if type == 'INT(v1dS)':
            res1 += P.integ(z, var1)[0]
        elif type == 'INT(v1.ndS)':
            res = P.integNorm(z,var1)[0]
            res1 += res[0]; res2 += res[1]; res3 += res[2]
        elif type == 'INT((v1,v2,v3).ndS)':
            res1 += P.integNormProduct(z,vector)
        elif type == 'INT((v1,v2,v3)^OMdS)':
            res = P.integMoment(z,center,vector)
            res1 += res[0]; res2 += res[1]; res3 += res[2]
        elif type == 'INT(v1n^OMdS)':
            res = P.integMomentNorm(z,center,var1)[0]
            res1 += res[0]; res2 += res[1]; res3 += res[2]
    if type == 'INT((v1,v2,v3)^OMdS)' or type == 'INT(v1n^OMdS)' or type == 'INT(v1.ndS)':
        res = [res1,res2,res3]
    else: res = res1
    CTK.TXT.insert('START', 'Res='+str(res)+'.\n')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkInteg  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Integrate fields.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkInteg')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Integration type -
    V = TK.StringVar(win); V.set('INT(v1dS)'); VARS.append(V)
    # -1- Moment center  -
    V = TK.StringVar(win); V.set('0;0;0'); VARS.append(V)
    # -2- Var0 for integration -
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -3- Var1 for integration -
    V = TK.StringVar(win); V.set('CoordinateY'); VARS.append(V)
    # -4- Var2 for integration -
    V = TK.StringVar(win); V.set('CoordinateZ'); VARS.append(V)

    # - Menu des integrations -
    B = TTK.OptionMenu(Frame, VARS[0], \
                       'INT(v1dS)', 'INT(v1.ndS)', \
                       'INT((v1,v2,v3).ndS)', \
                       'INT((v1,v2,v3)^OMdS)', 'INT(v1n^OMdS)')
    B.grid(row=0, column=0, sticky=TK.EW)

    # - Centre des moments -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Moment center.')

    # - Bouton compute -
    B = TTK.Button(Frame, text="Compute", command=compute)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Compute integration on contour/surface.')

    # - Menu des variables -
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[2], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList1)
        F.grid(row=1, column=0, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 1 (v1).')
        WIDGETS['variable1'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[2], 
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList1_2)
        F.grid(row=1, column=0, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 1 (v1).')
        WIDGETS['variable1'] = B

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[3], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=1, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 2 (v2).')
        WIDGETS['variable2'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[3], 
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2_2)
        F.grid(row=1, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 2 (v2).')
        WIDGETS['variable2'] = B

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[4], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList3)
        F.grid(row=1, column=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 3 (v3).')
        WIDGETS['variable3'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[4], 
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList3_2)
        F.grid(row=1, column=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 3 (v3).')
        WIDGETS['variable3'] = B

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['PostNoteBook'].add(WIDGETS['frame'], text='tkInteg')
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
    (win, menu, file, tools) = CTK.minimal('tkInteg '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
