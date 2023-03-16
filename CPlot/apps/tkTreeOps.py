# - tkTreeOps -
"""General operations on pyTree."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Converter.Filter as Filter
import numpy

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Cree une liste des bases de l'arbre (optionMenu)
#==============================================================================
def updateBaseNameList(event=None):
    bases = Internal.getBases(CTK.t)
    bs = ['New Base']
    for b in bases: bs.append(b[0])
    m = WIDGETS['base'].children['menu']
    m.delete(0, TK.END)
    for i in bs:
        m.add_command(label=i, command=lambda v=VARS[0],l=i:v.set(l))

#==============================================================================
# Cree une liste des bases de l'arbre (comboBox)
#==============================================================================
def updateBaseNameList2(event=None):
    bases = Internal.getBases(CTK.t)
    bs = []
    for b in bases: bs.append(b[0])
    if 'base' in WIDGETS:
        WIDGETS['base']['values'] = bs

#==============================================================================
# Move selected zone in a new base
#==============================================================================
def moveSelection():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    baseName = VARS[0].get()

    base = Internal.getNodeFromName1(CTK.t, baseName)
    if base is None: # if base doesnt exist, create it
        C._addBase2PyTree(CTK.t, baseName, 3)
        base = Internal.getNodeFromName1(CTK.t, baseName)

    CTK.saveTree()
    Z = []
    deletedZoneNames = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        zone = CTK.t[2][nob][2][noz]
        deletedZoneNames.append(CTK.t[2][nob][0]+Internal.SEP1+CTK.t[2][nob][2][noz][0])
        Z.append(zone)
    CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
    CPlot.delete(deletedZoneNames)
    nob = C.getNobOfBase(base, CTK.t)
    for i in Z: CTK.add(CTK.t, nob, -1, i)
    
    CTK.TXT.insert('START', 'Selection moved to %s.\n'%baseName)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def moveNodeUp():
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node is None or node[3] == 'CGNSTree_t': return # Tree node can not move
    (p, c) = Internal.getParentOfNode(CTK.t, node)
    if c == 0: return # already first
    
    if node[3] == 'Zone_t': # optimise
        z1 = p[2][c-1]
        z2 = p[2][c]
        if z1[3] != 'Zone_t':
            temp = p[2][c-1]
            p[2][c-1] = p[2][c]
            p[2][c] = temp
            CTK.TKTREE.updateApp()
            return
        (nob1, noz1) = C.getNobNozOfZone(z1, CTK.t)
        (nob2, noz2) = C.getNobNozOfZone(z2, CTK.t)
        CTK.replace(CTK.t, nob1, noz1, z2)
        CTK.replace(CTK.t, nob2, noz2, z1)
        CTK.TKTREE.updateApp()
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CPlot.render()
    elif node[3] == 'CGNSBase_t': # non optimise
        temp = p[2][c-1]
        p[2][c-1] = p[2][c]
        p[2][c] = temp
        CTK.TKTREE.updateApp()
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.display(CTK.t)
    else: # n'impacte pas CPlot
        temp = p[2][c-1]
        p[2][c-1] = p[2][c]
        p[2][c] = temp
        CTK.TKTREE.updateApp()

#==============================================================================
def moveNodeDown():
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node is None or node[3] == 'CGNSTree_t': return # Tree node can not move
    (p, c) = Internal.getParentOfNode(CTK.t, node)
    if c == len(p[2])-1: return # already last

    if node[3] == 'Zone_t': # optimise
        z1 = p[2][c+1]
        z2 = p[2][c]
        if z1[3] != 'Zone_t':
            temp = p[2][c+1]
            p[2][c+1] = p[2][c]
            p[2][c] = temp
            CTK.TKTREE.updateApp()
            return
        (nob1, noz1) = C.getNobNozOfZone(z1, CTK.t)
        (nob2, noz2) = C.getNobNozOfZone(z2, CTK.t)
        CTK.replace(CTK.t, nob2, noz2, z1)
        CTK.replace(CTK.t, nob1, noz1, z2)
        CTK.TKTREE.updateApp()
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CPlot.render()
    elif node[3] == 'CGNSBase_t': # non optimise
        temp = p[2][c+1]
        p[2][c+1] = p[2][c]
        p[2][c] = temp
        CTK.TKTREE.updateApp()
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.display(CTK.t)
    else: # n'impacte pas CPlot
        temp = p[2][c+1]
        p[2][c+1] = p[2][c]
        p[2][c] = temp
        CTK.TKTREE.updateApp()

# format string pour ecriture casual (6)
def strFormat(v):
    if isinstance(v, int): return "%ld"%v
    elif isinstance(v, float): return "%.5g"%v
    else: return "%s"%v
    
# format string pour ecriture precise (12)
def strFormat2(v):
    if isinstance(v, int): return "%ld"%v
    elif isinstance(v, float): return "%.12g"%v
    else: return "%s"%v

#==============================================================================
# Affiche la valeur de noeud et le type du noeud
#==============================================================================
def showNodeValue(event=None):
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []:
        index = 0; VARS[1].set(str(index))
        VARS[2].set(''); VARS[3].set(''); VARS[4].set(''); return

    VARS[3].set(node[3]); VARS[4].set(node[0])
    
    index = VARS[1].get()
    try: index = int(index)
    except: index = 0; VARS[1].set('0')
    if index < 0: index = 0; VARS[1].set(str(index))
    v = node[1]
    if isinstance(v, float): v = strFormat2(v); index = 0
    elif isinstance(v, int): v = strFormat2(v); index = 0
    elif isinstance(v, numpy.ndarray):
        if v.dtype == 'c': v = Internal.getValue(node); index = 0
        else:
            sh = v.shape
            if len(sh) == 3:
                dim = '('+str(sh[0])+';'+str(sh[1])+';'+str(sh[2])+'): '
            elif len(sh) == 2:
                dim = '('+str(sh[0])+';'+str(sh[1])+'): '
            elif len(sh) == 1:
                dim = '('+str(sh[0])+'): '
            else:
                dim = str(sh)+': '
            
            try: vf = v.ravel(order='F')
            except: vf = v.flat
            if index >= v.size: index = v.size-1; VARS[1].set(str(index))
            # Construit une chaine representant le tableau a plat
            if index > 3:
                if v.size <= index+1: flatView = ''
                elif v.size <= index+2:
                    flatView = ' '+strFormat(vf[index+1])+''
                elif v.size <= index+3:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+''
                else:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+'...'
                flatView += '\n'
                CTK.TXT.insert('START', flatView)
                flatView = strFormat(vf[index])
                CTK.TXT.insert('START', flatView, 'Warning')
                flatView = dim+'...'+strFormat(vf[index-2])+' '+strFormat(vf[index-1])+' '
                CTK.TXT.insert('START', flatView)
                
            elif index == 2:
                if v.size <= index+1: flatView = ''
                elif v.size <= index+2:
                    flatView = ' '+strFormat(vf[index+1])+''
                elif v.size <= index+3:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+''
                else:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+'...'
                flatView += '\n'
                CTK.TXT.insert('START', flatView)
                flatView = strFormat(vf[index])
                CTK.TXT.insert('START', flatView, 'Warning')
                flatView = dim+strFormat(vf[index-2])+' '+strFormat(vf[index-1])+' '
                CTK.TXT.insert('START', flatView)
                  
            elif index == 1:
                if v.size <= index+1: flatView = ''
                elif v.size <= index+2:
                    flatView = ' '+strFormat(vf[index+1])+''
                elif v.size <= index+3:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+''
                else:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+'...'
                flatView += '\n'
                CTK.TXT.insert('START', flatView)
                flatView = strFormat(vf[index])
                CTK.TXT.insert('START', flatView, 'Warning')
                flatView = dim+strFormat(vf[index-1])+' '
                CTK.TXT.insert('START', flatView)
                 
            elif index == 0:
                if v.size <= index+1: flatView = ''
                elif v.size <= index+2:
                    flatView = ' '+strFormat(vf[index+1])+''
                elif v.size <= index+3:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+''
                else:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+'...'
                flatView += '\n'
                CTK.TXT.insert('START', flatView)
                flatView = strFormat(vf[index])
                CTK.TXT.insert('START', flatView, 'Warning')
                flatView = dim
                CTK.TXT.insert('START', flatView)                
            v = strFormat2(vf[index])
    VARS[1].set(str(index))
    VARS[2].set(v)

#==============================================================================
# increase index in inspector
#==============================================================================
def increaseIndex(event=None):
    if CTK.t == []: return
    index = VARS[1].get()
    index = int(index)
    index += 1
    VARS[1].set(str(index))
    showNodeValue()

#==============================================================================
# decrease index in inspector
#==============================================================================
def decreaseIndex(event=None):
    if CTK.t == []: return
    index = VARS[1].get()
    index = int(index)
    index -= 1
    VARS[1].set(str(index))
    showNodeValue()
    
#==============================================================================
# Changement de la valeur du noeud
#==============================================================================
def setNodeValue(event=None):
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    index = VARS[1].get()
    try: index = int(index)
    except: index = 0
    val = VARS[2].get()
    v = node[1]
    if isinstance(v, float):
        Internal.setValue(node, val)
        CTK.TXT.insert('START', 'Value "%d" set in node.\n.'%val) # not done by showNodeValue
    elif isinstance(v, int):
        Internal.setValue(node, val)
        CTK.TXT.insert('START', 'Value "%g" set in node.\n.'%val) # not done by showNodeValue
    elif isinstance(v, numpy.ndarray):
        if v.dtype == 'c':
            Internal.setValue(node, val)
            CTK.TXT.insert('START', 'Value "%s" set in node.\n.'%val) # not done by showNodeValue
        elif v.dtype == 'int32':
            try: vf = v.ravel(order='F')
            except: vf = v.flat
            vf[index] = int(val)
            #Internal.setValue(node, node[1])
        else:
            try: vf = v.ravel(order='F')
            except: vf = v.flat
            vf[index] = float(val)
            #Internal.setValue(node, node[1])
            CTK.display(CTK.t) # view may be impacted
    showNodeValue()

#==============================================================================
# Change le type du noeud. Ca peut etre tres dangereux.
#==============================================================================
def setNodeType(event=None):
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    node[3] = VARS[3].get()
    
#==============================================================================
# Change le nom du noeud.
#==============================================================================
def setNodeName(event=None):
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    node[0] = VARS[4].get()
    CTK.TKTREE.updateApp()

#==============================================================================
# Appele de tkTree, pour informer qu'un noeud a ete selectionne
#==============================================================================
def updateNode(node):
    showNodeValue()

#==============================================================================
def rmNodes():
    if CTK.t == []: return

#==============================================================================
def setByLevel(node1, node2, depth, maxDepth):
    #print("setting", node1[0])
    node1[1] = node2[1]
    if depth < maxDepth or maxDepth == -1:
        for c, n in enumerate(node1[2]):
            setByLevel(n, node2[2][c], depth+1, maxDepth)

#==============================================================================
def freeByLevel(node1, depth, maxDepth):
    #print("freeing", node1[0])
    node1[1] = None
    if depth < maxDepth or maxDepth == -1:
        for n in node1[2]:
            freeByLevel(n, depth+1, maxDepth)

#==============================================================================
def loadNode():
    if CTK.t == []: return
    if CTK.HANDLE is None: fileName = CTK.FILE
    else: fileName = CTK.HANDLE.fileName
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    path = Internal.getPath(CTK.t, node)
    depth = VARS[7].get()
    depth = int(depth)
    nodes = Filter.readNodesFromPaths(fileName, [path], maxDepth=depth)
    
    # depth replace
    if depth == 0: node[1] = nodes[0][1]
    elif depth == -1: 
        node[1] = nodes[0][1]
        node[2] = nodes[0][2]
    else: 
        setByLevel(node, nodes[0], 0, depth)
    updateNode(node)

#==============================================================================
def freeNode():
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    depth = VARS[7].get()
    depth = int(depth)
    freeByLevel(node, 0, depth)
    updateNode(node)

#==============================================================================
def saveNode():
    if CTK.t == []: return
    if CTK.HANDLE is None: fileName = CTK.FILE
    else: fileName = CTK.HANDLE.fileName
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    path = Internal.getPath(CTK.t, node)
    depth = VARS[7].get()
    depth = int(depth)
    Filter.writeNodesFromPaths(fileName, [path], [node], maxDepth=depth, mode=1)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkTreeOps', font=CTK.FRAMEFONT, takefocus=1)
    
    #BB = CTK.infoBulle(parent=Frame, text='Operations on pyTree.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    Frame.columnconfigure(3, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkTreeOps')
    WIDGETS['frameMenu'] = FrameMenu
    
    # - VARS -
    # -0- Base to move to -
    V = TK.StringVar(win); V.set('New Base'); VARS.append(V)
    # -1- Displayed index
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -2- Displayed node value
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -3- Displayed node type
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -4- Displayed node name
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -5- Name/Type to be deleted
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -6- Type of deletion
    V = TK.StringVar(win); V.set('Name'); VARS.append(V)
    # -7- Maxdepth for load
    V = TK.StringVar(win); V.set('0'); VARS.append(V)

    # - Move selection to base -
    B = TTK.Button(Frame, text="Move sel. to", command=moveSelection)
    BB = CTK.infoBulle(parent=B, text='Move selected zones to another base.')
    B.grid(row=0, column=0, sticky=TK.EW)

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[0], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateBaseNameList)
    else:
        VARS[0].set('newBase')
        B = ttk.Combobox(F, textvariable=VARS[0], 
                         values=[], state='normal')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateBaseNameList2)
    F.grid(row=0, column=1, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Destination Base name.')
    WIDGETS['base'] = B

    # - Move nodes up/down/import selection in tkTree
    B = TTK.Button(Frame, text="Node Up", command=moveNodeUp)
    BB = CTK.infoBulle(parent=B, text='Move selected node in tkTree up.')
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Button(Frame, text="Node Down", command=moveNodeDown)
    BB = CTK.infoBulle(parent=B, text='Move selected node in tkTree down.')
    B.grid(row=1, column=1, columnspan=3, sticky=TK.EW)

    # - Rm Nodes -
    #B = TK.Button(Frame, text="Rm node(s)", command=rmNodes)
    #BB = CTK.infoBulle(parent=B, text='Remove nodes from tree.\nYou can also use the Suppr key.')
    #B.grid(row=2, column=0, sticky=TK.EW)
    #B = TK.OptionMenu(Frame, VARS[6], 'Name', 'Type', 'Value')
    #B.grid(row=2, column=1, sticky=TK.EW)
    #B = TK.Entry(Frame, textvariable=VARS[5], background='White', width=5)
    #BB = CTK.infoBulle(parent=B, text='Enter node name or types. Wildcards are accepted.')
    #B.grid(row=2, column=2, sticky=TK.EW)

    # - Small node editor -
    F = TTK.LabelFrame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE,
                       text='Edit node', takefocus=1)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=4)
    F.grid(row=3, column=0, columnspan=4, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[4], background='White', width=5)
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Selected node name.\nYou can change this by pressing <Enter>.')
    B.bind('<Return>', setNodeName)

    B = TTK.Entry(F, textvariable=VARS[3], background='White', width=5)
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Selected node type.\nYou can change this by pressing <Enter>.')
    B.bind('<Return>', setNodeType)
    
    FrameA = TK.Frame(F)
    FrameA.columnconfigure(0, weight=0)
    FrameA.columnconfigure(1, weight=1)
    FrameA.columnconfigure(2, weight=0)
    FrameA.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(FrameA, textvariable=VARS[1], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Displayed index (start: 0).')
    B.bind('<Return>', showNodeValue)
    B = TTK.Button(FrameA, text="<", command=decreaseIndex, width=2)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(FrameA, text=">", command=increaseIndex, width=2)
    B.grid(row=0, column=2, sticky=TK.EW)
    
    B = TTK.Entry(F, textvariable=VARS[2], background='White', width=5)
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Selected node value at given index.\nYou can change this by pressing <Enter>.')
    B.bind('<Return>', setNodeValue)

    # - Load/free Nodes -
    B = TTK.Button(Frame, text="Load node", command=loadNode)
    BB = CTK.infoBulle(parent=B, text='Load selected node from file.')
    B.grid(row=4, column=0, sticky=TK.EW)
    #B = TTK.Button(Frame, text="Save node", command=saveNode)
    #BB = CTK.infoBulle(parent=B, text='Save selected node.')
    #B.grid(row=4, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Button(Frame, text="Free", command=freeNode)
    BB = CTK.infoBulle(parent=B, text='Free selected node.')
    B.grid(row=4, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[7], background='White', width=2)
    B.grid(row=4, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Max depth for load, save and free (-1: full).')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.NSEW)

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
    (win, menu, file, tools) = CTK.minimal('tkTreeOps '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
