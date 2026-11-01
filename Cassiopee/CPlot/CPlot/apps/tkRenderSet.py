# - tkRenderSet -
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Converter.Internal as Internal
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

MATERIALS = ['Solid', 'Flat', 'Glass', 'Chrome',
             'Metal', 'Wood', 'Marble', 'Granite', 'Brick', 'XRay',
             'Cloud', 'Gooch', 'Sphere', 'Texmat']

COLORS = ['White', 'Black', 'Grey', 'Blue', 'Red', 'Green', 'Yellow',
          'Orange', 'Brown', 'Magenta', 'Custom>']

#==============================================================================
# Appele quand une zone color est selectionnee (optionMenu)
def setColorVar(l):
    if l == 'Custom>':
        import tkinter.colorchooser as tkColorChooser
        ret = tkColorChooser.askcolor()
        l = ret[1]
    VARS[1].set(l)

#==============================================================================
# Appele quand une zone color est selectionnee (combobox)
def setColorVar2(event=None):
    l = VARS[1].get()
    if l == 'Custom>':
        import tkinter.colorchooser as tkColorChooser
        ret = tkColorChooser.askcolor()
        l = ret[1]
        VARS[1].set(l)

#==============================================================================
# Appele quand une mesh color est selectionnee (optionMenu)
def setMeshColorVar(l):
    if l == 'Custom>':
        import tkinter.colorchooser as tkColorChooser
        ret = tkColorChooser.askcolor()
        l = ret[1]
    VARS[7].set(l)

#==============================================================================
# Appele quand une mesh color est selectionnee (combobox)
def setMeshColorVar2(event=None):
    l = VARS[7].get()
    if l == 'Custom>':
        import tkinter.colorchooser as tkColorChooser
        ret = tkColorChooser.askcolor()
        l = ret[1]
        VARS[7].set(l)

#==============================================================================
# Cree la liste des zone colors + variables (optionMenu)
#==============================================================================
def updateVarNameList(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        zvars = C.getVarNames(CTK.dt, excludeXYZ=True, mode=1)
    else:
        zvars = C.getVarNames(CTK.t, excludeXYZ=True, mode=1)
    m = WIDGETS['colors'].children['menu']
    m.delete(0, TK.END)
    allvars = COLORS
    if len(zvars) > 0:
        for v in zvars[0]: allvars.append('Iso:'+v)
    for i in allvars:
        m.add_command(label=i, command=lambda v=VARS[1],l=i:setColorVar(l))

#==============================================================================
# Cree la liste des zone colors + variables (combobox)
#==============================================================================
def updateVarNameList2(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        zvars = C.getVarNames(CTK.dt, excludeXYZ=True, mode=1)
    else:
        zvars = C.getVarNames(CTK.t, excludeXYZ=True, mode=1)
    allvars = COLORS
    if len(zvars) > 0:
        for v in zvars[0]: allvars.append('Iso:'+v)
    if 'colors' in WIDGETS:
        WIDGETS['colors']['values'] = allvars

#==============================================================================
# Cree la liste des mesh colors (optionMenu)
#==============================================================================
def updateMeshColorList(event=None):
    m = WIDGETS['meshColors'].children['menu']
    m.delete(0, TK.END)
    allvars = COLORS
    for i in allvars:
        m.add_command(label=i, command=lambda v=VARS[1],l=i:setMeshColorVar(l))

#==============================================================================
# Cree la liste des mesh colors (combobox)
#==============================================================================
def updateMeshColorList2(event=None):
    allvars = COLORS
    if 'colors' in WIDGETS:
        WIDGETS['meshColors']['values'] = allvars

#==============================================================================
# set zone material
#==============================================================================
def setMaterial():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    material = VARS[0].get()
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], material=material)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    Panels.updateRenderPanel()
    CPlot.render()

#==============================================================================
# set zone color
#==============================================================================
def setColor():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    color = VARS[1].get()
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], color=color)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    Panels.updateRenderPanel()
    CPlot.render()

#==============================================================================
# set mesh color
#==============================================================================
def setMeshColor():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    color = VARS[7].get()
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], meshColor=color)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    Panels.updateRenderPanel()
    CPlot.render()

#==============================================================================
# set mesh width
#==============================================================================
def setMeshWidth(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    width = VARS[8].get()
    width = float(width)
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], meshWidth=width)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    Panels.updateRenderPanel()
    CPlot.render()

#==============================================================================
def setAll():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    material = VARS[0].get()
    color = VARS[1].get()
    blending = WIDGETS['blending'].get() / 100.
    VARS[6].set('Blending [%.2f].'%blending)
    shaderParameter2 = (WIDGETS['param2'].get()) / 50.
    shaderParameter1 = (WIDGETS['param1'].get()) / 50.
    meshOverlay = VARS[3].get()
    meshColor = VARS[7].get()
    meshWidth = VARS[8].get()
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz],
                                 material=material,
                                 color=color,
                                 blending=blending,
                                 shaderParameters=[shaderParameter1,
                                                   shaderParameter2])
        if meshOverlay == "1":
            CPlot._addRender2Zone(a,
                                  meshOverlay=meshOverlay,
                                  meshColor=meshColor,
                                  meshWidth=meshWidth)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    Panels.updateRenderPanel()
    CPlot.render()

#==============================================================================
def setBlending(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    blending = WIDGETS['blending'].get() / 100.
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    VARS[6].set('Blending [%.2f]'%(WIDGETS['blending'].get() / 100.))

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], blending=blending)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    Panels.updateRenderPanel()
    CPlot.render()

#==============================================================================
def setMeshOverlay():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    meshOverlay = VARS[3].get()
    if meshOverlay == "1":
        WIDGETS['meshFrame'].grid(row=1, column=0, columnspan=2)
    else:
        WIDGETS['meshFrame'].grid_remove()

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz],
                                 meshOverlay=meshOverlay)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    Panels.updateRenderPanel()
    CPlot.render()

#==============================================================================
def setShaderParameter(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    shaderParameter2 = (WIDGETS['param2'].get()) / 50.
    shaderParameter1 = (WIDGETS['param1'].get()) / 50.
    VARS[4].set('Shader parameter 1 [%.2f].'%shaderParameter1)
    VARS[5].set('Shader parameter 2 [%.2f].'%shaderParameter2)

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz],
                                 shaderParameters=[shaderParameter1,
                                                   shaderParameter2])
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    Panels.updateRenderPanel()
    CPlot.render()

#==============================================================================
def getData():
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []: return
    # get first of selection
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    zone = CTK.t[2][nob][2][noz]
    if zone is not None:
        n = Internal.getNodeFromPath(zone, '.RenderInfo/Material')
        if n is not None:
            val = Internal.getValue(n)
            VARS[0].set(val)
        n = Internal.getNodeFromPath(zone, '.RenderInfo/Color')
        if n is not None:
            val = Internal.getValue(n)
            VARS[1].set(val)
        n = Internal.getNodeFromPath(zone, '.RenderInfo/Blending')
        if n is not None:
            val = Internal.getValue(n)
            WIDGETS['blending'].set(val*100)
            VARS[6].set('Blending [%.2f].'%val)
        n = Internal.getNodeFromPath(zone, '.RenderInfo/MeshOverlay')
        if n is not None:
            val = Internal.getValue(n)
            VARS[3].set(val)
            if val == "1": WIDGETS['meshFrame'].grid(row=1, column=0, columnspan=2)
            else: WIDGETS['meshFrame'].grid_remove()
        n = Internal.getNodeFromPath(zone, '.RenderInfo/MeshColor')
        if n is not None:
            val = Internal.getValue(n)
            VARS[7].set(val)
        n = Internal.getNodeFromPath(zone, '.RenderInfo/MeshWidth')
        if n is not None:
            val = Internal.getValue(n)
            val = str(val)
            VARS[8].set(val)
        n = Internal.getNodeFromPath(zone, '.RenderInfo/ShaderParameters')
        if n is not None:
            val = Internal.getValue(n)
            WIDGETS['param1'].set(val[0]*50)
            VARS[4].set('Shader parameter 1 [%.2f].'%val[0])
            WIDGETS['param2'].set(val[1]*50)
            VARS[4].set('Shader parameter 2 [%.2f].'%val[1])

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkRenderSet  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Customize block rendering\n(render mode).\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkRenderSet')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Material -
    V = TK.StringVar(win); V.set('Solid'); VARS.append(V)
    # -1- Color -
    V = TK.StringVar(win); V.set('White'); VARS.append(V)
    # -2- Blending
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -3- Mesh overlay
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -4- Shader parameter1 info bulle
    V = TK.StringVar(win); V.set('Shader parameter 1.'); VARS.append(V)
    # -5- Shader parameter2 info bulle
    V = TK.StringVar(win); V.set('Shader parameter 2.'); VARS.append(V)
    # -6- Blending info bulle
    V = TK.StringVar(win); V.set('Blending.'); VARS.append(V)
    # -7- Mesh color
    V = TK.StringVar(win); V.set('Black'); VARS.append(V)
    # -8- Mesh width
    V = TK.StringVar(win); V.set('1'); VARS.append(V)

    # - set Blending -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL,
                  command=setBlending, showvalue=0, borderwidth=1, value=100)
    WIDGETS['blending'] = B
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[6])

    # - set mesh overlay -
    B = TTK.Checkbutton(Frame, text='Mesh', variable=VARS[3],
                        command=setMeshOverlay)
    BB = CTK.infoBulle(parent=B, text='Add underlaying mesh.')
    B.grid(row=0, column=1, sticky=TK.EW)

    # - set mesh color -
    F2 = TTK.Frame(Frame, borderwidth=0)
    F2.columnconfigure(0, weight=1)
    F2.columnconfigure(1, weight=1)
    B = TTK.Button(F2, text="Mesh Color", command=setMeshColor)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the color of mesh lines.')

    F = TTK.Frame(F2, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[7], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateMeshColorList)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['meshColors'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[7],
                         values=[], state='readonly', height=11)
        B.bind('<<ComboboxSelected>>', setMeshColorVar2)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateMeshColorList2)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['meshColors'] = B
    #F2.grid(row=1, column=0, columnspan=2)
    WIDGETS['meshFrame'] = F2

    # - set mesh width -
    B = TTK.Button(F2, text="Mesh Width", command=setMeshWidth)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the width of mesh lines.')

    B = TTK.Entry(F2, textvariable=VARS[8], background='White', width=10)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Width of mesh lines.')
    B.bind('<Return>', setMeshWidth)
    WIDGETS['meshWidth'] = B

    # - set shader parameter 1 & 2-
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL,
                  command=setShaderParameter, showvalue=0, borderwidth=1, value=50)
    WIDGETS['param1'] = B
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[4])

    # - set shader parameter 2 -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL,
                  command=setShaderParameter, showvalue=0, borderwidth=1, value=50)
    WIDGETS['param2'] = B
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[5])

    # - set Material -
    B = TTK.Button(Frame, text="Set Material", command=setMaterial)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the material render property.')
    B = TTK.OptionMenu(Frame, VARS[0], *MATERIALS)
    B.grid(row=3, column=1, sticky=TK.EW)

    # - set Color -
    B = TTK.Button(Frame, text="Set Color", command=setColor)
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the color render property.')

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[1], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList)
        F.grid(row=4, column=1, sticky=TK.EW)
        WIDGETS['colors'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[1],
                         values=[], state='readonly', height=11)
        B.bind('<<ComboboxSelected>>', setColorVar2)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=4, column=1, sticky=TK.EW)
        WIDGETS['colors'] = B

    # - set all properties -
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=0)

    B = TTK.Button(F, text="Set all properties", command=setAll)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set all render properties.')
    B = TTK.Button(F, command=getData,
                   image=iconics.PHOTO[8], padx=0, pady=0)
    BB = CTK.infoBulle(parent=B, text='Get properties from selected zone.')
    B.grid(row=0, column=1, sticky=TK.EW)
    F.grid(row=5, column=0, columnspan=2, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['RenderNoteBook'].add(WIDGETS['frame'], text='tkRenderSet')
    except: pass
    CTK.WIDGETS['RenderNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['RenderNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkRenderSet '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
