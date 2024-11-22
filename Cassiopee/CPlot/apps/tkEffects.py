# - tkEffects -
"""Setup render effects (shadow, dof, ...)"""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Converter.Internal as Internal
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setShadow(event=None):
    if VARS[1].get() == '1': shadow= 1 
    else: shadow = 0 
    CPlot.setState(shadow=shadow)
    CPlot._addRender2PyTree(CTK.t, shadow=shadow)

#==============================================================================
def setDOF(event=None):
    if VARS[2].get() == '1': dof = 1
    else: dof = 0 
    CPlot.setState(dof=dof)
    CPlot._addRender2PyTree(CTK.t, dof=dof)

#==============================================================================
def setDofPower(event=None):
    off = WIDGETS['dofPower'].get()*4./100.
    VARS[5].set('Depth of field power [%.2f].'%off)
    CPlot.setState(dofPower=off)
    CPlot._addRender2PyTree(CTK.t, dofPower=off)

#==============================================================================
def setSharpenPower(event=None):
    off = WIDGETS['sharpenPower'].get()*2./100.
    VARS[8].set('Sharpen power [%.2f].'%off)
    CPlot.setState(sharpenPower=off)
    CPlot._addRender2PyTree(CTK.t, sharpenPower=off)

#==============================================================================
def setLightOffsetX(event=None):
    off = WIDGETS['lightOffsetX'].get() / 50. - 1.
    off = -off
    VARS[3].set('Light offset in x [%.2f %%].'%off)
    CPlot.setState(lightOffset=(off, -999))
    CPlot._addRender2PyTree(CTK.t, lightOffsetX=off)

#==============================================================================
def setLightOffsetY(event=None):
    off = WIDGETS['lightOffsetY'].get() / 50.
    VARS[4].set('Light offset in y [%.2f %%].'%off)
    CPlot.setState(lightOffset=(-999, off))
    CPlot._addRender2PyTree(CTK.t, lightOffsetY=off)

#==============================================================================
def setViewAngle(event=None):
    angle = VARS[0].get()
    angle = float(angle)
    CPlot.setState(viewAngle=angle)

#==============================================================================
def setGammaCorrection(event=None):
    off = WIDGETS['gammaCorrection'].get()*2.5/100.
    VARS[6].set('Gamma correction [%.2f].'%off)
    CPlot.setState(gamma=off)
    CPlot._addRender2PyTree(CTK.t, gamma=off)

#==============================================================================
def setToneMapping(event=None):
    ntype = VARS[7].get()
    rtype = 0
    if ntype == 'None': rtype = 0
    elif ntype == 'ACE': rtype = 1
    elif ntype == 'Filmic': rtype = 2
    elif ntype == 'Uchimura': rtype = 3
    CPlot.setState(toneMapping=rtype)
    CPlot._addRender2PyTree(CTK.t, tone=rtype)

#==============================================================================
def getData(event=None):
    renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
    if renderInfo is None: return
    slot = Internal.getNodeFromName1(renderInfo, 'Slot0')
    if slot is None: return

    pos = Internal.getNodeFromName1(renderInfo, 'materials')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i))
        CPlot.setState(materials=out)
    pos = Internal.getNodeFromName1(renderInfo, 'bumpMaps')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i))
        CPlot.setState(bumpMaps=out)  
    pos = Internal.getNodeFromName1(renderInfo, 'billBoards')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i)) #; out += [1,1]
        CPlot.setState(billBoards=out, billBoardSize=0.8)

    pos = Internal.getNodeFromName1(slot, 'dof')
    if pos is not None:
        dof = Internal.getValue(pos)
        VARS[2].set(dof)
        CPlot.setState(dof=dof)

    pos = Internal.getNodeFromName1(slot, 'tone')
    if pos is not None:
        tone = Internal.getValue(pos)
        if tone == 0: VARS[7].set('None')
        elif tone == 1: VARS[7].set('ACE')
        elif tone == 2: VARS[7].set('Filmic')
        elif tone == 3: VARS[7].set('Uchimura')
        CPlot.setState(toneMapping=tone)

    pos = Internal.getNodeFromName1(slot, 'gamma')
    if pos is not None:
        gamma = Internal.getValue(pos)
        WIDGETS['gammaCorrection'].set(100/2.5)
        CPlot.setState(gamma=gamma)

    pos = Internal.getNodeFromName1(slot, 'shadow')
    if pos is not None:
        shadow = Internal.getValue(pos)
        VARS[1].set(shadow)
        CPlot.setState(shadow=shadow)

    pos = Internal.getNodeFromName1(slot, 'lightOffsetX')
    if pos is not None:
        val = Internal.getValue(pos)
        WIDGETS['lightOffsetX'].set(-val * 50. + 1)
        CPlot.setState(lightOffset=(val,-999))

    pos = Internal.getNodeFromName1(slot, 'lightOffsetY')
    if pos is not None:
        val = Internal.getValue(pos)
        WIDGETS['lightOffsetY'].set(val * 50.)
        CPlot.setState(lightOffset=(-999,val))

    pos = Internal.getNodeFromName1(slot, 'dofPower')
    if pos is not None:
        dofPower = Internal.getValue(pos)
        WIDGETS['dofPower'].set(100./4.*dofPower)
        CPlot.setState(dofPower=dofPower)

    pos = Internal.getNodeFromName1(slot, 'sharpenPower')
    if pos is not None:
        sharpenPower = Internal.getValue(pos)
        WIDGETS['sharpenPower'].set(100./2.*sharpenPower)

        CPlot.setState(sharpenPower=sharpenPower)

    CTK.TXT.insert('START', 'Effects data loaded.\n')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkEffects  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Set anaglyph mode.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkEffects')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- cam view angle -
    V = TK.StringVar(win); V.set('50.'); VARS.append(V)
    if 'tkEffectsAngle' in CTK.PREFS: V.set(CTK.PREFS['tkEffectsAngle'])
    # -1- Shadow
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkEffectsShadow' in CTK.PREFS:
        V.set(CTK.PREFS['tkEffectsShadow'])
    # -2- DOF
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkEffectsDOF' in CTK.PREFS: V.set(CTK.PREFS['tkEffectsDOF'])
    # -3- Light offset X info bulle
    V = TK.StringVar(win); V.set('Light offset in x.'); VARS.append(V)
    # -4- Light offset Y info bulle
    V = TK.StringVar(win); V.set('Light offset in y.'); VARS.append(V)
    # -5- Depth of field info bulle
    V = TK.StringVar(win); V.set('Depth of field power.'); VARS.append(V)
    # -6- Gamma correction info bulle
    V = TK.StringVar(win); V.set('Gamma correction.'); VARS.append(V)
    # -7- Type of tone mapping
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    # -8- Sharpen info bulle
    V = TK.StringVar(win); V.set('Sharpen power.'); VARS.append(V)

    # - Camera angle -
    B = TTK.Button(Frame, text="Cam angle", command=setViewAngle)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the camera view angle.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
    B.bind('<Return>', setViewAngle)
    BB = CTK.infoBulle(parent=B, text='Cam view angle (deg).')

    # - Shadow -
    B = TTK.Checkbutton(Frame, text='Shadow', variable=VARS[1],
                        command=setShadow)
    BB = CTK.infoBulle(parent=B, text='Toggle shadows.')
    B.grid(row=1, column=0, sticky=TK.EW)

    # - Light offset X -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setLightOffsetX, value=50)
    WIDGETS['lightOffsetX'] = B
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[3])

    # - Light offset Y -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setLightOffsetY, value=0)
    WIDGETS['lightOffsetY'] = B
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[4])

    # - Post -
    B = TTK.Checkbutton(Frame, text='Post', variable=VARS[2], command=setDOF)
    BB = CTK.infoBulle(parent=B, text='Toggle image post processing (depth of field blur/gamma correction/tone mapping).')
    B.grid(row=2, column=0, sticky=TK.EW)

    # - Tone mapping -
    B = TTK.OptionMenu(Frame, VARS[7], 'None', 'ACE', 'Filmic', 'Uchimura', command=setToneMapping)
    B.grid(row=2, column=1, sticky=TK.EW)

    # - Gamma correction -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setGammaCorrection, value=50)
    WIDGETS['gammaCorrection'] = B
    B.grid(row=2, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[6])

    # - DOF -
    B = TTK.Label(Frame, text='DOF')
    BB = CTK.infoBulle(parent=B, text='Depth of field blur (around cliked point).')
    B.grid(row=3, column=0, sticky=TK.EW)

    # - DOF power -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setDofPower, value=0)
    WIDGETS['dofPower'] = B
    B.grid(row=3, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[5])

    # - Sharpen -
    B = TTK.Label(Frame, text='Sharp')
    BB = CTK.infoBulle(parent=B, text='Sharpen image.')
    B.grid(row=4, column=0, sticky=TK.EW)

    # - Sharpen power -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setSharpenPower, value=0)
    WIDGETS['sharpenPower'] = B
    B.grid(row=4, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[8])

    # - Get data from tree -
    B = TTK.Button(Frame, command=getData,
                   image=iconics.PHOTO[8], padx=0, pady=0)
    BB = CTK.infoBulle(parent=B, text='Get data from tree.')
    B.grid(row=4, column=2, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['RenderNoteBook'].add(WIDGETS['frame'], text='tkEffects')
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
def saveApp():
    CTK.PREFS['tkEffectsAngle'] = VARS[0].get()
    CTK.PREFS['tkEffectsShadow'] = VARS[1].get()
    CTK.PREFS['tkEffectsDOF'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('50.')
    VARS[1].set('0')
    VARS[2].set('0')
    CTK.PREFS['tkEffectsAngle'] = VARS[0].get()
    CTK.PREFS['tkEffectsShadow'] = VARS[1].get()
    CTK.PREFS['tkEffectsDOF'] = VARS[2].get()
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
    (win, menu, file, tools) = CTK.minimal('tkEffects '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
