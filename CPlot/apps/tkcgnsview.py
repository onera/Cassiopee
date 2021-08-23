# -- cassiopee cgnsview main app --
try: import Tkinter as TK
except: import tkinter as TK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import os
import os.path, sys

#==============================================================================
# To be called when CTK.t is set
def run():
    
    # Main window
    (win, frames, menu, menus, file, tools) = CTK.minimal2('kcgnsview '+C.__version__,
                                                           show=False)
    fileName = os.path.split(CTK.FILE)
    CTK.changeWindowTitle(fileName[1], fileName[0])

    # Place win devant les autres fenetres
    win.deiconify(); win.focus_set()

    # open tkTreeOps
    app = 'tkTreeOps'    
    CTK.TKMODULES[app] = None
    CTK.TKMODULEFRAMES[app] = frames[0]
    #name = app; name = '  '+name
    #menu[0].add_command(label=name, command=lambda x=app:CTK.openApp(x))
    CTK.openApp(app)

    # - Update apps -    
    CTK.TKTREE.updateApp()
    #if 'tkContainers' in CTK.TKMODULES: CTK.TKMODULES['tkContainers'].updateApp()
    
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
        import Converter.Filter as Filter
        CTK.HANDLE = Filter.Handle(files[0])
        CTK.t = CTK.HANDLE.loadSkeleton(maxDepth=-1)
        CTK.FILE = files[0]
    run()    
