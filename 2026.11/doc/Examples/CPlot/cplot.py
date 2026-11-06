#!/usr/bin/env python
# Pg principal
import time
def infLoop():
    a = False
    while (a == False): time.sleep(10000)

def CPlotDisplay__(file):
    import os.path
    ext = os.path.splitext(file)
    extension = ext[len(ext)-1]
    if extension == '.plt' or extension == '.PLT':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'bin_tp')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.dat' or extension == '.tp' or extension == '.DAT' or extension == '.TP':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'fmt_tp')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.v3d' or extension == '.V3D':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'bin_v3d')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.mesh' or extension == '.MESH':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'fmt_mesh')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.stl' or extension == '.STL':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'bin_stl')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.fig' or extension == '.FIG':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'fmt_xfig')
        CPlot.display(a, dim=2); CPlot.setFileName__(file); infLoop();
    elif extension == '.svg' or extension == '.SVG':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'fmt_svg')
        CPlot.display(a, dim=2); CPlot.setFileName__(file); infLoop();
    elif extension == '.pov' or extension == '.POV':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'fmt_pov')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.3ds' or extension == '.3DS':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'bin_3ds')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.obj' or extension == '.OBJ':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'fmt_obj')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.gts' or extension == '.GTS':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'fmt_gts')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.cgns' or extension == '.CGNS':
        import Converter.PyTree as CP; import CPlot.PyTree
        a = CP.convertFile2PyTree(file, 'bin_cgns')
        CPlot.PyTree.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.adf' or extension == '.ADF':
        import Converter.PyTree as CP; import CPlot.PyTree
        a = CP.convertFile2PyTree(file, 'bin_cgns')
        CPlot.PyTree.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.hdf' or extension == '.HDF':
        import Converter.PyTree as CP; import CPlot.PyTree
        a = CP.convertFile2PyTree(file, 'bin_hdf')
        CPlot.PyTree.display(a); CPlot.setFileName__(file); infLoop();
    elif extension == '.pickle' or extension == '.PICKLE' or extension[0:4] == '.ref' or extension[0:4] == '.REF':
        import Converter as C; import CPlot
        a = C.convertFile2Arrays(file, 'bin_pickle')
        CPlot.display(a); CPlot.setFileName__(file); infLoop();

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print ("cplot: cplot <file>")
    else:
        file = sys.argv[1]
        CPlotDisplay__(file)
