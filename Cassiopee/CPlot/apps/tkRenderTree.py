# - tkRenderTree -
# Add global rendering options to tree
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Converter.Internal as Internal
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def addBillboardFile(event=None):
    if CTK.t == []: return
    v = VARS[0].get()
    v.replace(' ', '')
    v = v.split(';')
    vo = []
    for i in v:
        if i != '': vo.append(i)
    CPlot._addRender2PyTree(CTK.t, billBoards=vo)
    CTK.TKTREE.updateApp()

    renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
    pos = Internal.getNodeFromName1(renderInfo, 'billBoards')
    out = []
    for i in pos[2]: out.append(Internal.getValue(i))
    CPlot.setState(billBoards=out, billBoardSize=0.8)
    CTK.TXT.insert('START', 'New billboard texture added.\n')

#==============================================================================
def addBumpMapFile(event=None):
    if CTK.t == []: return
    v = VARS[0].get()
    v.replace(' ', '')
    v = v.split(';')
    vo = []
    for i in v:
        if i != '': vo.append(i)
    CPlot._addRender2PyTree(CTK.t, bumpMaps=vo)
    CTK.TKTREE.updateApp()

    renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
    pos = Internal.getNodeFromName1(renderInfo, 'bumpMaps')
    out = []
    for i in pos[2]: out.append(Internal.getValue(i))
    CPlot.setState(bumpMaps=out)
    CTK.TXT.insert('START', 'New bumpmap texture added.\n')

#==============================================================================
def addTextureFile(event=None):
    if CTK.t == []: return
    # Get files to add
    v = VARS[0].get()
    v.replace(' ', '')
    v = v.split(';')
    vo = []
    for i in v:
        if i != '': vo.append(i)
    # Get previous materials
    renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
    pos = Internal.getNodeFromName1(renderInfo, 'materials')
    if pos is not None: nprev = len(pos[2])
    else: nprev = 0

    # Update shaderParameters
    for z in Internal.getZones(CTK.t):
        info = Internal.getNodeFromName1(z, '.RenderInfo')
        if info is not None:
            material = Internal.getNodeFromName1(info, 'Material')
            if material is not None and Internal.getValue(material) == 'Texmat':
                params = Internal.getNodeFromName1(info, 'ShaderParameters')[1]
                params[1] = params[1]*nprev/(nprev+len(vo))

    # Add new materials
    CPlot._addRender2PyTree(CTK.t, materials=vo)
    CTK.TKTREE.updateApp()

    # Update CPlot textures
    renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
    pos = Internal.getNodeFromName1(renderInfo, 'materials')
    out = []
    for i in pos[2]: out.append(Internal.getValue(i))
    CPlot.setState(materials=out)
    CPlot.display(CTK.t)
    CTK.TXT.insert('START', 'New color texture added.\n')

#==============================================================================
def chooseFile():
    try: import tkFileDialog
    except: import tkinter.filedialog as tkFileDialog
    import os.path
    init = VARS[0].get()
    init = init.split(';')[0]
    files = tkFileDialog.askopenfilenames(
        filetypes=CTK.fileTypes, initialfile=init, multiple=1)
    if files == '' or files == None or files == (): # user cancel
        return
    # strangely, initfile is part of the return
    files = CTK.fixFileString__(files, init)
    s = ''
    for f in files:
        bname = os.path.basename(f)
        if os.path.exists(bname): s += bname+';' # short
        else: s += f+';' # full name
    s = s[:-1]
    VARS[0].set(s)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkRenderTree  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Add global rendering options to tree.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkRenderTree')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Name of file to add -
    V = TK.StringVar(win); V.set('file.png'); VARS.append(V)

    # - Buttons -
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=0)
    B = TTK.Entry(F, textvariable=VARS[0], background='White')
    BB = CTK.infoBulle(parent=B, text='Name of image file to add as billboard or texture.')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text="...", padx=0, command=chooseFile)
    BB = CTK.infoBulle(parent=B, text='Select an image file (png).')
    B.grid(row=0, column=1, sticky=TK.EW)
    F.grid(row=0, column=0, sticky=TK.EW)

    B = TTK.Button(Frame, text="Add Base Color textures", command=addTextureFile)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Add given files as color texture file to be referenced\n in Texmat material of tkRenderSet.')
    B = TTK.Button(Frame, text="Add Normal maps", command=addBumpMapFile)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Add given files as normal/bumpmap file to be referenced\n in Texmat material of tkRenderSet.')
    B = TTK.Button(Frame, text="Add Billboards textures", command=addBillboardFile)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Add file names as image billboard file to be referenced\n in Sphere material of tkRenderSet.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['RenderNoteBook'].add(WIDGETS['frame'], text='tkRenderTree')
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
    (win, menu, file, tools) = CTK.minimal('tkRenderTree '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
