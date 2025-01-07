# Search BAR
try: import tkinter as TK
except: import Tkinter as TK
from . import Ttk as TTK
from . import Tk as CTK

# Dictionaire Action->applet
applet = {
    # tkNodeEdir
    'tkNodeEdit':'tkNodeEdit', 'Edit/change node value':'tkNodeEdit', 'Load a node':'tkNodeEdit',
    # tkTreeOps
    'tkTreeOps':'tkTreeOps', 'Move node':'tkTreeOps', 'Delete/remove node':'tkTreeOps',
    'Move zone to another base':'tkTreeOps',
    # tkCheckPyTree
    'tkCheckPyTree':'tkCheckPyTree', 'Correct pyTree':'tkCheckPyTree',
    'Correct nodes':'tkCheckPyTree', 'Correct BCs':'tkCheckPyTree',
    # tkFilter
    'tkFilter':'tkFilter', 'Find multigrid zones':'tkFilter',
    'Find zone of given name':'tkFilter', 'Find zone of given proc':'tkFilter',
    # tkFamily
    'tkFamily':'tkFamily', 'Create/tag zone family':'tkFamily',
    'Create BCFamily':'tkFamily',
    # tkCADFix
    'tkCADFix':'tkCADFix', 'Read/load CAD file':'tkCADFix',
    'Write/save CAD file':'tkCADFix', 'Repair CAD':'tkCADFix',
    # tkState
    'tkState':'tkState', 'Set problem dimension':'tkState',
    'Set reference state':'tkState',
    # tkPrefs
    'tkPrefs':'tkPrefs', 'Change background image/color':'tkPrefs',
    'Change export resolution':'tkPrefs', 'Customize':'tkPrefs',
    # tkPerfo
    'tkPerfo':'tkPerfo', 'Improve performance':'tkPerfo',
    'Set number of threads':'tkPerfo', 'Approximate mesh display':'tkPerfo',
    # tkContainers
    'tkContainers':'tkContainers', 'My fields dont display!':'tkContainers',
    'Change FlowSolution': 'tkContainers',
    # tkCamera
    'tkCamera':'tkCamera', 'Camera position':'tkCamera', 'View position':'tkCamera',
    # tkRuler
    'tkRuler':'tkRuler', 'Measure distance': 'tkRuler',
    # tkFind
    'tkFind':'tkFind', 'Find index/cell':'tkFind', 'Extract one/some cells':'tkFind',
    # tkProbe
    'tkProbe':'tkProbe', 'Probe value in cell or vertex':'tkProbe', 'Set a field value in a cell':'tkProbe',
    # tkCanvas
    'tkCanvas':'tkCanvas', 'Create Canvas for drawing':'tkCanvas',
    'Enlarge canvas':'tkCanvas', 'Reduce canvas':'tkCanvas',
    # tkPoint
    'tkPoint':'tkPoint', 'Create a point':'tkPoint',
    'Modify point coordinates':'tkPoint',
    # tkDraw
    'tkDraw':'tkDraw', 'Draw lines/beziers/splines':'tkDraw',
    'Draw circle/rectangle':'tkDraw',
    'Draw freely':'tkDraw', 'Draw on surfaces':'tkDraw',
    # tkExtractEdges
    'tkExtractEdges':'tkExtractEdges', 'Get external edges':'tkExtractEdges',
    'Get two edges intersection':'tkExtractEdges', 'Get sharp edges':'tkExtractEdges',
    'Convert BAR to Struct':'tkExtractEdges', 'Split TBranches':'tkExtractEdges',
    'BAR2Struct':'tkExtractEdges',
    # tkMapEdge
    'tkMapEdges':'tkMapEdge', 'Uniformize edge distrib':'tkMapEdge',
    'Refine edge distrib':'tkMapEdge', 'Smooth edge distrib':'tkMapEdge',
    'Copy a distrib to another edge':'tkMapEdge',
    'Enforce step in edge':'tkMapEdge',
    # tkBasicSurfs
    'tkBasicSurfs':'tkBasicSurfs', 'Create Sphere/Cube/Cylinder surface mesh':'tkBasicSurfs',
    'Create Tetra/Hexa/Pyramid/Torus surface mesh':'tkBasicSurfs',
    'Create Plane mesh':'tkBasicSurfs',
    'Create Classical surface meshes':'tkBasicSurfs', 'Create mesh from Geometry data base':'tkBasicSurfs',
    # tkText
    'tkText':'tkText', 'Create mesh of text':'tkText',
    # tkCADMesh
    'tkCADMesh':'tkCADMesh', 'Generate CAD surface mesh':'tkCADMesh',
    'Mesh CAD edges':'tkCADMesh',
    # tkFixer
    'tkFixer2':'tkFixer2', 'Fill hole in surface':'tkFixer2',
    'ConformUnstr':'tkFixer2', 'Conformize a TRI surface':'tkFixer2',
    # tkBoolean
    'tkBoolean':'tkBoolean', 'Union two surfaces':'tkBoolean',
    'Intersection of rwo surfaces':'tkBoolean',
    'Difference of two surfaces':'tkBoolean',
    # tkMapUV
    'tkMapUV':'tkMapUV', 'UV map of surface':'tkMapUV',
    # tkCartWrap
    'tkCartWrap':'tkCartWrap', 'Remesh a surface with cartesian wrapper':'tkCartWrap',
    'Wrap surface -watertight-':'tkCartWrap', 'Cartesian wrapper':'tkCartWrap',
    # tkOffset
    'tkOffset':'tkOffset', 'Offset a surface of a given distance':'tkOffset',
    # tkMMGs
    'tkMMGs':'tkMMGs', 'Remesh a TRI surface':'tkMMGs',
    'Refine a TRI surface':'tkMMGs',
    # tkSurfaceWalk
    'tkSurfaceWalk':'tkSurfaceWalk',
    'Mesh a surface -othogonal/structured-':'tkSurfaceWalk',
    # tkProjection
    'tkProjection':'tkProjection', 'Project a mesh on a surface':'tkProjection',
    'Project ortho':'tkProjection', 'Project following view':'tkProjection',
    # tkCells
    'tkCells':'tkCells', 'Suppress cells':'tkCells', 'Refine cells':'tkCells',
    'Select some cells':'tkCells',
    # tkStretch
    'tkStretch':'tkStretch', 'Refine edge/mesh':'tkStretch',
    'Remesh edge/mesh':'tkStretch', 'Enforce step in mesh':'tkStretch',
    # tkExtrusion
    'tkExtrusion':'tkExtrusion', 'Add k-planes':'tkExtrusion',
    'Extrude with normals':'tkExtrusion', 'Revolve mesh':'tkExtrusion',
    'Axisym mesh':'tkExtrusion',
    # tkTetraMesher
    'tkTetraMesher':'tkTetraMesher', 'Fill with tetras':'tkTetraMesher',
    'Fill with triangles':'tkTetraMesher', 'Triangulate mesh':'tkTetraMesher',
    'Generate tetra/triangle mesh':'tkTetraMesher',
    # tkTFI
    'TkTFI':'tkTFI', 'TFI Mesh -structured-':'tkTFI',
    'O Mesh -structured-':'tkTFI',
    # tkSmooth
    'tkSmooth':'tkSmooth', 'Smooth mesh':'tkSmooth',
    # tkOctree
    'tkOctree':'tkOctree', 'Create Octree -unstructured-':'tkOctree',
    'Create Octree -structured-':'tkOctree',
    # tkMeshQual
    'tkMeshQual':'tkMeshQual', 'Check mesh quality':'tkMeshQual',
    'Check mesh regularity':'tkMeshQual',
    'Check mesh orthogonality':'tkMeshQual',
    'Compute cell volume':'tkMeshQual',
    'View negative volume cells':'tkMeshQual',
    # tkMeshInfo
    'tkMeshInfo':'tkMeshInfo', 'Mesh number of points/cells':'tkMeshInfo',
    'Display variables min/max':'tkMeshInfo',
    # tkBlock
    'tkBlock':'tkBlock', 'Delete/remove block':'tkBlock', 'Convert to tetra':'tkBlock',
    'Convert to hexa':'tkBlock', 'Convert to node':'tkBlock',
    'Exterior faces':'tkBlock', 'Close block':'tkBlock',
    'Suppress multiple points':'tkBlock', 'Take one over n points':'tkBlock',
    'Element type conversion':'tkBlock',
    # tkTransform
    'tkTransform':'tkTransform', 'Translate block':'tkTransform',
    'Rotate block':'tkTransform', 'Scale block':'tkTransform',
    'Mirror block':'tkTransform',
    'Stretch block':'tkTransform',
    # tkNGon
    'tkNGon':'tkNGon', 'Convert to NGON':'tkNGon', 'Convert to elements':'tkNGon',
    'Dual of mesh':'tkNGon', 'Conformize faces (NGON)':'tkNGon',
    'Break elements':'tkNGon',
    # tkGhostCells
    'tkGhostCells':'tkGhostCells', 'Add/remove ghost cells':'tkGhostCells',
    # tkSplit
    'tkSplit':'tkSplit', 'Split blocks':'tkSplit', 'Join blocks':'tkSplit',
    'Subzone block':'tkSplit', 'Split in connex parts':'tkSplit',
    'Split in manifold parts':'tkSplit',
    # tkReorder
    'tkReorder':'tkReorder', 'Reorder block':'tkReorder', 'Make direct':'tkReorder',
    # tkBC
    'tkBC':'tkBC', 'Set BC':'tkBC', 'View BC':'tkBC', 'Connect match':'tkBC',
    'Set Boundary conditions':'tkBC', 'Fill empty BCs':'tkBC',
    'View undefined BC':'tkBC', 'Remove BC':'tkBC',
    # tkIBC
    'tkIBC':'tkIBC', 'Set snear on surface':'tkIBC', 'Immersed boundaries':'tkIBC', 'Set data for IBM':'tkIBC',
    # tkChimera
    'tkChimera':'tkChimera', 'Blank cells':'tkChimera',
    'Optimize overlap':'tkChimera',
    # tkExtractBC
    'tkExtractBC':'tkExtractBC', 'Extract BC':'tkExtractBC',
    # tkInit
    'tkInit':'tkInit', 'Init field from reference state':'tkInit',
    # tkDistributor
    'tkDistributor':'tkDistributor', 'Distribute over processors':'tkDistributor',
    # tkDist2Walls
    'tkDist2Walls':'tkDist2Walls', 'Compute wall distance':'tkDist2Walls',
    # tkTime
    'tkTime':'tkTime', 'View time motion':'tkTime',
    # tkRigidMotion
    'tkRigidMotion':'tkRigidMotion', 'Set rigid motion in tree':'tkRigidMotion',
    # tkElsaSolver
    'tkElsaSolver':'tkElsaSolver', 'Create elsAHybrid':'tkElsaSolver',
    'Adapt tree for elsA':'tkElsaSolver',
    # tkFastSolver
    'tkFastSolver':'tkFastSolver', 'Compute CFD with IBM':'tkFastSolver',
    # tkVariables
    'tkVariables':'tkVariables', 'Compute variables/fields':'tkVariables',
    'Center2Node or Node2Center variables/fields':'tkVariables',
    'Grad/curl variables/fields':'tkVariables', 'Remove/delete variables/fields':'tkVariables',
    'Compute Pressure/Mach/Vorticity':'tkVariables',
    'Init variables/fields -with formula-':'tkVariables',
    # tkExtractMesh
    'tkExtractMesh':'tkExtractMesh', 'Interpolate fields on mesh':'tkExtractMesh',
    # tkStream
    'tkStream':'tkStream', 'Compute stream/ribbon/surface':'tkStream',
    # tkIsoSurf
    'tkIsoSurf':'tkIsoSurf', 'Compute iso surface':'tkIsoSurf',
    # tkInteg
    'tkInteg':'tkInteg', 'Integrate field on mesh':'tkInteg',
    # tkView
    'tkView':'tkView', 'View mesh':'tkView', 'View field/variable':'tkView',
    'Save view point':'tkView', 'View normal orientation':'tkView',
    'Display legend':'tkView', 'Switch 2D/3D view':'tkView',
    'Swtich to render mode':'tkView',
    # tkPlotXY
    'tkPlotXY':'tkPlotXY', 'Plot curves':'tkPlotXY', 'Graph':'tkPlotXY',
    # tkSlice
    'tkSlice':'tkSlice', 'Slice mesh':'tkSlice', 'View inside mesh':'tkSlice',
    'Cut mesh':'tkSlice',
    # tkIJK
    'tkIJK':'tkIJK', 'View IJK planes':'tkIJK',
    # tkCellN
    'tkCellN':'tkCellN', 'View blanking':'tkCellN', 'View chimera data':'tkCellN',
    'View Orphan points -chimera-':'tkCellN',
    # tkBackground
    'tkBackground':'tkBackground', 'Add a background mesh':'tkBackground',
    # tkRender
    'tkRenderTree':'tkRenderTree', 'Load textures':'tkRenderTree',
    'tkRenderSet':'tkRenderSet', 'Set surface material': 'tkRenderSet',
    'Set Chrome/Wood/Glass/Stone effect on surface':'tkRenderSet',
    'Set XRay/Metal/Gooch/Smoke effect on surface':'tkRenderSet',
    'Set Iso-color+effect on surface':'tkRenderSet',
    'Set color on surface':'tkRenderSet',
    'Set Surface texture':'tkRenderSet', 'Set Transparency/blending':'tkRenderSet',
    # tkStereo
    'tkStereo':'tkStereo', 'View anaglyph':'tkStereo', '3D effect':'tkStereo',
    'Use Red/Blue glasses':'tkStereo',
    # tkEffects
    'tkEffects':'tkEffects', 'Add shadow':'tkEffects', 'Change camera angle':'tkEffects',
    'Add depth of field': 'tkEffects', 'Set gamma': 'tkEffects', 'Set camera angle': 'tkEffects',
    # tkDemo
    'tkDemo':'tkDemo', 'Automatic camera motion':'tkDemo'
}
# Liste des actions
lista = applet.keys()

# Class entry avec auto-completion
class AutocompleteEntry(TK.Entry):
    def __init__(self, lista, *args, **kwargs):
        TK.Entry.__init__(self, *args, **kwargs)
        TK.Entry.config(self, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
        self.lista = lista
        self.var = self["textvariable"]
        if self.var == '':
            self.var = self["textvariable"] = TK.StringVar()
            self.var.set("Ask me...")
        self.var.trace('w', self.changed)
        self.bind("<FocusIn>", self.focusin)
        self.bind("<Right>", self.selection)
        self.bind("<Up>", self.up)
        self.bind("<Down>", self.down)
        self.bind("<Return>", self.selection)
        self.bind("<Control-c>", self.clearVar)
        self.bind("<Control-u>", self.clearVar)
        self.bind("<Escape>", self.clearVar)
        self.lb_up = False

    def clearVar(self, event):
        self.var.set('')

    def focusin(self, event):
        if self.var.get() == 'Ask me...': self.var.set('')

    def changed(self, name, index, mode):
        if self.var.get() == '':
            if self.lb_up: self.lb.destroy()
            self.lb_up = False
        else:
            words = self.comparison()
            if words:
                if not self.lb_up: # listbox exists
                    self.lb = TTK.Listbox(width=self.winfo_width())
                    self.lb.bind("<Double-Button-1>", self.selection)
                    self.lb.bind("<Right>", self.selection)
                    self.lb.bind("<Return>", self.selection)
                    self.lb.place(x=self.winfo_x(), y=self.winfo_y()+2*self.winfo_height()+11)
                    self.lb_up = True
                    #self.sb = TTK.Scrollbar()
                    #self.lb.config(yscrollcommand=self.sb.set)
                    #self.sb.config(command=self.lb.yview)
                    #self.sb.place(x=self.winfo_x()-10, y=self.winfo_y()+self.winfo_height())

                self.lb.delete(0, TK.END)
                for w in words:
                    self.lb.insert(TK.END,w)
            else:
                if self.lb_up:
                    self.lb.destroy()
                    self.lb_up = False

    def selection(self, event):
        if self.lb_up:
            self.var.set(self.lb.get(TK.ACTIVE))
            self.lb.destroy()
            self.lb_up = False
            self.icursor(TK.END)
        self.enter()

    def up(self, event):
        if self.lb_up:
            if self.lb.curselection() == (): index = '0'
            else: index = self.lb.curselection()[0]
            if index != '0':
                self.lb.selection_clear(first=index)
                index = str(int(index)-1)
                self.lb.selection_set(first=index)
                self.lb.activate(index)

    def down(self, event):
        if self.lb_up:
            if self.lb.curselection() == (): index = '-1'
            else: index = self.lb.curselection()[0]
            if index != TK.END:
                self.lb.selection_clear(first=index)
                index = str(int(index)+1)
                self.lb.selection_set(first=index)
                self.lb.activate(index)

    def comparison(self):
        import re
        askString = self.var.get()
        askString = askString.split(' ')
        pattern = []; la = len(askString)
        for a in askString:
            pattern.append(re.compile('.*' + a + '.*', re.IGNORECASE))

        sol = [[]]*la

        for w in self.lista:
            ma = 0
            for p in pattern:
                if re.match(p, w): ma += 1
            if ma > 0: sol[la-ma].append(w)

        ret = []
        for s in sol:
            lr = len(ret); ls = len(s)
            if lr < 8:
                if lr + ls < 8: ret += s
                else: ret += s[0:8-lr]
        #ret = [w for w in self.lista if re.match(pattern, w)]
        #import difflib
        #word =  self.var.get()
        #ret = difflib.get_close_matches(word, self.lista, n=8, cutoff=0.4)
        return ret

    def enter(self, event=None):
        word = self.var.get()
        if word in applet:
            # Get applet name
            app = applet[word]
            # activate the APP
            module = CTK.getModule(app)
            module.showApp()
            frame = None; menu = None
            from tkCassiopee import TREEAPPS, STATEAPPS, EDGEAPPS, SURFAPPS, MESHAPPS, BLOCKAPPS, BCAPPS, MOTIONAPPS, SOLVERAPPS, POSTAPPS, VISUAPPS, RENDERAPPS
            frames = CTK.WIDGETS['noteBookFrames']
            menus = CTK.WIDGETS['noteBookMenus']
            buttons = CTK.WIDGETS['noteBookButtons']
            if app in TREEAPPS:
                frame = frames[0]; menu = menus[0]; bt = buttons[0]
            elif app in STATEAPPS:
                frame = frames[1]; menu = menus[1]; bt = buttons[1]
            elif app in EDGEAPPS:
                frame = frames[2]; menu = menus[2]; bt = buttons[2]
            elif app in SURFAPPS:
                frame = frames[3]; menu = menus[3]; bt = buttons[3]
            elif app in MESHAPPS:
                frame = frames[4]; menu = menus[4]; bt = buttons[4]
            elif app in BLOCKAPPS:
                frame = frames[5]; menu = menus[5]; bt = buttons[5]
            elif app in BCAPPS:
                frame = frames[6]; menu = menus[6]; bt = buttons[6]
            elif app in MOTIONAPPS:
                frame = frames[7]; menu = menus[7]; bt = buttons[7]
            elif app in SOLVERAPPS:
                frame = frames[8]; menu = menus[8]; bt = buttons[8]
            elif app in POSTAPPS:
                frame = frames[9]; menu = menus[9]; bt = buttons[9]
            elif app in VISUAPPS:
                frame = frames[10]; menu = menus[10]; bt = buttons[10]
            elif app in RENDERAPPS:
                frame = frames[11]; menu = menus[11]; bt = buttons[11]
            if frame is not None:
                CTK.WIDGETS['noteBook'].display(frame, menu)
                for i in buttons: TTK.deselectRadioButton(i)
                TTK.selectRadioButton(bt)

# Cree la search bar (hidden)
def createSearchBar():
    T = TK.Toplevel(border=0)
    T.withdraw()
    F = TTK.Frame(T, takefocus=1)
    F.bind('<Enter>', lambda event : F.focus_set())
    F.columnconfigure(0, weight=1)
    entry = AutocompleteEntry(lista, F)
    entry.grid(row=0, column=0, sticky=TK.EW)
    F.grid(row=0, column=0, sticky=TK.EW)
    return T

def createSearchBar2(Frame):
    entry = AutocompleteEntry(lista, Frame)
    return entry

if __name__ == '__main__':
    root = TK.Tk()
    root.columnconfigure(0, weight=1)
    F = TTK.Frame(root, takefocus=1)
    F.bind('<Enter>', lambda event : F.focus_set())
    F.columnconfigure(0, weight=1)
    F.grid(row=0, column=0, sticky=TK.EW)
    entry = AutocompleteEntry(lista, F)
    entry.grid(row=0, column=0, sticky=TK.EW)
    root.mainloop()
