# -- Cassiopee main app --
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import os, os.path, sys

# Liste des apps par sous menu et perso
TREEAPPS = ['tkNodeEdit', 'tkTreeOps', 'tkCheckPyTree', '---',
            'tkFamily', '---',
            'tkCADFix']
STATEAPPS = ['tkState', '---',
             'tkPrefs', 'tkPerfo', 'tkContainers', 'tkCamera', '---',
             #'tkLogFile', '---',
             'tkFilter', 'tkFind', 'tkProbe', '---',
             'tkRuler']
EDGEAPPS = ['tkCanvas', 'tkPoint', 'tkDraw','---',
            'tkExtractEdges', 'tkMapEdge']
SURFAPPS = ['tkBasicSurfs', 'tkText', '---',
            'tkCADMesh', 
            'tkFixer2', 'tkBoolean', '---',
            'tkMapUV', 'tkSculpt', '---',
            'tkMMGs', 'tkCartWrap', 'tkOffset', 'tkSurfaceWalk', '---',
            'tkProjection']
MESHAPPS = ['tkCells', 'tkStretch', '---',
            #tkMirabelle,
            'tkExtrusion', 'tkTetraMesher', 'tkTFI', 'tkSmooth', '---',
            'tkOctree', 'tkCollarMesh', 'tkBlader', '---',
            'tkMeshQual', 'tkMeshInfo']
BLOCKAPPS = ['tkBlock', '---',
             'tkTransform', 'tkNGon', 'tkGhostCells', '---',
             'tkSplit', 'tkReorder']
BCAPPS = ['tkBC', '---',
          'tkChimera', 'tkIBC',
          #'tkIBC2',
          '---',
          'tkExtractBC']
MOTIONAPPS = ['tkRigidMotion', 'tkTime']
SOLVERAPPS = ['tkInit', 'tkDistributor', 'tkDist2Walls', '---',
              #tkCassiopeeSolver,
              'tkElsaSolver', 'tkFastSolver']
POSTAPPS = ['tkVariables', '---',
            'tkExtractMesh', '---',
            'tkStream', 'tkIsoLine', 'tkIsoSurf', '---',
            'tkInteg']
VISUAPPS = ['tkView', #'tkPlot', 
            'tkPlotXY', '---',
            'tkSlice', 'tkIJK', 'tkCellN', '---',
            'tkBackground']
RENDERAPPS = ['tkRenderTree', 'tkRenderSet', '---',
              'tkStereo', 'tkEffects', 'tkDemo', '---',
              'tkPovRay', 'tkLuxRender']

ALLAPPS = TREEAPPS + STATEAPPS + EDGEAPPS + SURFAPPS + MESHAPPS + \
          BLOCKAPPS + BCAPPS + MOTIONAPPS + SOLVERAPPS + POSTAPPS + \
          VISUAPPS + RENDERAPPS
PERSOAPPS = []

#==============================================================================
# Add a personal app to pref file
#==============================================================================
def addPersonalApp():
    try: import tkinter.filedialog as tkFileDialog 
    except: import tkFileDialog
    file = tkFileDialog.askopenfilename(filetypes=[('python', '*.py')])
    a = os.access(file, os.F_OK)
    if not a: return
    CTK.loadPrefFile()
    if 'module' in CTK.PREFS: CTK.PREFS['module'] += ' ;'+file
    else: CTK.PREFS['module'] = file
    CTK.savePrefFile()
    file = os.path.split(file)
    moduleName = os.path.splitext(file[1])[0]
    pathName = file[0]
    try:
        orig = sys.path; local = orig
        local.append(pathName)
        sys.path = local
        module = __import__(moduleName)
        CTK.TKMODULES[moduleName] = module
        sys.path = orig
        #tools.add_command(label=moduleName,
        #                  command=module.showApp)
        frame = CTK.TKMODULEFRAMES['tkTreeOps']
        module.createApp(frame)
        PERSOAPPS.append(moduleName)
        CTK.TXT.insert('START', 'Module %s added in tools menu (next restart)\n'%moduleName)
    except:
        CTK.TXT.insert('START', 'Can not import '+moduleName+'.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
def notImplemented():
    CTK.TXT.insert('START', 'This functionality is not implemented.\n')
    CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# To be called when CTK.t is set
#==============================================================================
def run(t=None):

    if t is not None:
        if Internal.isTopTree(t): CTK.t = t
        else: CTK.t, ntype = Internal.node2PyTree(t)

    if CTK.t != []:
        # upgrade tree
        CTK.t = CTK.upgradeTree(CTK.t)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        fileName = os.path.split(CTK.FILE)
        CPlot.setWindowTitle(fileName[1], fileName[0])

    # - Verifie l'arbre -
    errors = []
    if CTK.t != []:
        errors = Internal.checkPyTree(CTK.t, level=5)
        if errors == []: CTK.display(CTK.t)

    # Load and set prefs for interface
    CTK.loadPrefFile(); CTK.setPrefs()

    # Modules a ouvrir automatiquement
    auto = {}
    for app in ALLAPPS:
        app = app.split('/')
        if len(app) == 2: app = app[1]
        else: app = app[0]
        auto[app] = 0
    if 'auto' in CTK.PREFS:
        p = CTK.PREFS['auto']; p = p.split(';')
        for i in p:
            i = i.strip()
            auto[i] = 1

    # Main window
    (win, frames, menu, menus, file, tools) = CTK.minimal2('Cassiopee '+C.__version__,
                                                           show=False)

    fileName = os.path.split(CTK.FILE)
    CTK.changeWindowTitle(fileName[1], fileName[0])

    # - Apps -
    submenus = {}
    for app in TREEAPPS: CTK.addMenuItem(app, menus[0], frames[0], submenus, auto)
    submenus = {}
    for app in STATEAPPS: CTK.addMenuItem(app, menus[1], frames[1], submenus, auto)
    submenus = {}
    for app in EDGEAPPS: CTK.addMenuItem(app, menus[2], frames[2], submenus, auto)
    submenus = {}
    for app in SURFAPPS: CTK.addMenuItem(app, menus[3], frames[3], submenus, auto)
    submenus = {}
    for app in MESHAPPS: CTK.addMenuItem(app, menus[4], frames[4], submenus, auto)
    submenus = {}
    for app in BLOCKAPPS: CTK.addMenuItem(app, menus[5], frames[5], submenus, auto)
    submenus = {}
    for app in BCAPPS: CTK.addMenuItem(app, menus[6], frames[6], submenus, auto)
    submenus = {}
    for app in MOTIONAPPS: CTK.addMenuItem(app, menus[7], frames[7], submenus, auto)
    submenus = {}
    for app in SOLVERAPPS: CTK.addMenuItem(app, menus[8], frames[8], submenus, auto)
    submenus = {}
    for app in POSTAPPS: CTK.addMenuItem(app, menus[9], frames[9], submenus, auto)
    submenus = {}
    for app in VISUAPPS: CTK.addMenuItem(app, menus[10], frames[10], submenus, auto)
    submenus = {}
    for app in RENDERAPPS: CTK.addMenuItem(app, menus[11], frames[11], submenus, auto)

    # Updated Apps from tree (containers from tree containers)
    module = CTK.getModule('tkContainers'); module.updateApp()

    # Get tkPlotXY if any
    try:
        module = CTK.getModule('tkPlotXY')
        if module is not None: CTK.TKPLOTXY = module
    except: pass

    # - Personal apps  -
    tools.add_command(label='Add a personal app',
                      command=addPersonalApp)

    if 'module' in CTK.PREFS:
        mod = CTK.PREFS['module']
        mod = mod.split(';')
        for i in mod:
            i = i.strip()
            file = os.path.split(i)
            moduleName = os.path.splitext(file[1])[0]
            pathName = file[0]
            try:
                orig = sys.path; local = orig
                local.append(pathName)
                sys.path = local
                module = __import__(moduleName)
                CTK.TKMODULES[moduleName] = module
                sys.path = orig
                tools.add_command(label=moduleName,
                                  command=module.showApp)
                module.createApp(frames[0]); module.hideApp()
                PERSOAPPS.append(moduleName)
            except:
                CTK.TXT.insert('START', 'can not import '+moduleName+'.\n')
                CTK.TXT.insert('START', 'Error: ', 'Error')

    # Place win devant les autres fenetres
    win.deiconify(); win.focus_set()

    # - Erreur dans l'arbre -
    if errors != []:
        Panels.displayErrors(errors, header='Checking pyTree')
        CTK.t = Internal.correctPyTree(CTK.t, level=5)
        CTK.display(CTK.t)

    # - Update apps -    
    CTK.TKTREE.updateApp()
    if CTK.TKMODULES['tkContainers'] is not None: CTK.TKMODULES['tkContainers'].updateApp()
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

    # - open load panel if partial load -
    if CTK.t != []:
        zones = Internal.getZones(CTK.t)
        if len(zones) == 0: # all skeletons certainely
            Panels.openLoadPanel()

    # Load textures, billboards file names into CPlot
    CPlot.loadImageFiles(CTK.t)

    # Automatically adjust view if any slot
    if CTK.t != []:
        renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
        if renderInfo is not None:
            CTK.getModule("tkView")
            import tkView; tkView.loadSlot()

    # - Main loop -
    win.mainloop()

    # Del photos
    CTK.PHOTOS = []    

#==============================================================================
if __name__ == "__main__":
    # Ouverture du fichier de la ligne de commande
    import sys
    if len(sys.argv) >= 2:
        files = sys.argv[1:]
        CTK.tkLoadFile(files, mode='auto')
    run()
