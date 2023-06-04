# -- cassiopee cgnsview main app --
import Converter.PyTree as C
import CPlot.Tk as CTK
import os, os.path, sys

#==============================================================================
# To be called when CTK.t is set
def run(a, q):
    
    # Set cassiopee prefs
    CTK.loadPrefFile(); CTK.setPrefs()

    # Passe la dimension par une pseudo pref
    CTK.PREFS['tkTreeWidth'] = 380
    CTK.PREFS['tkTreeHeight'] = 290

    # Main window
    (win, frames, menu, menus, file, tools) = CTK.minimal2('kcgnsview '+C.__version__,
                                                           show=False)
    fileName = os.path.split(CTK.FILE)
    CTK.changeWindowTitle(fileName[1], fileName[0])

    # Add some apps
    auto = {}
    auto['tkTreeOps'] = True
    auto['tkMeshInfo'] = True
    auto['tkPlotXY'] = False
    submenus = {}
    for app in ['tkTreeOps']: CTK.addMenuItem(app, menus[0], frames[0], submenus, auto)
    submenus = {}
    for app in ['tkMeshInfo']: CTK.addMenuItem(app, menus[4], frames[4], submenus, auto)
    try:
        import tkPlotXY
        submenus = {}
        for app in ['tkPlotXY']: CTK.addMenuItem(app, menus[10], frames[10], submenus, auto)
    except: pass
    
    # Place win devant les autres fenetres
    win.deiconify(); win.focus_set()

    # get data from reading process
    if a is not None and q is not None:
        CTK.t = q.get()
        a.join()

    # - Update apps -    
    CTK.TKTREE.updateApp()
    #if 'tkContainers' in CTK.TKMODULES: CTK.TKMODULES['tkContainers'].updateApp()
    
    # - Main loop -
    win.mainloop()

    # Del photos
    CTK.PHOTOS = []


#==============================================================================
# load data skeleton in a separate process
# IN: q: queue for data export
# IN: h: file handle
def loadSkeleton(q, h):
    t2 = h.loadSkeleton(maxDepth=-1)
    q.put(t2)
    
#==============================================================================
if __name__ == "__main__":
    # Ouverture du fichier de la ligne de commande
    import sys
    a = None; q = None
    if len(sys.argv) >= 2:
        files = sys.argv[1:]
        import Converter.Filter as Filter
        CTK.HANDLE = Filter.Handle(files[0])
        CTK.FILE = files[0]
        
        # Threading is limited to one core because of GIL
        #import threading
        #a = threading.Thread(target=loadFullSkeleton)
        #a.start()
        
        # with multiprocessing, we must communicate CTK.t through Queue
        import multiprocessing
        q = multiprocessing.Queue()
        a = multiprocessing.Process(target=loadSkeleton, args=(q,CTK.HANDLE))
        a.start()

    run(a, q)