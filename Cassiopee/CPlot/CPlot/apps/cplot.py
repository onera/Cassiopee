#!/usr/bin/env python
import sys
import os.path
import time
import CPlot

# FILE_FORMATS: dictionary mapping each file extension to
#   - its format
#   - whether it is a PyTree (True, import Converter.PyTree) or an array (False, import Converter)
#   - the display dimension.
FILE_FORMATS = {
    ".plt": ("bin_tp", False, 3),
    ".dat": ("fmt_tp", False, 3),
    ".tp":  ("fmt_tp", False, 3),
    ".v3d": ("bin_v3d", False, 3),
    ".mesh": ("fmt_mesh", False, 3),
    ".stl": ("bin_stl", False, 3),
    ".fig": ("fmt_xfig", False, 2),
    ".svg": ("fmt_svg", False, 2),
    ".pov": ("fmt_pov", False, 3),
    ".3ds": ("bin_3ds", False, 3),
    ".obj": ("fmt_obj", False, 3),
    ".gts": ("fmt_gts", False, 3),
    ".cgns": ("bin_cgns", True, 3),
    ".adf": ("bin_cgns", True, 3),
    ".hdf": ("bin_hdf", True, 3),
    ".pickle":  ("bin_pickle", False, 3),
}

def infLoop():
    while True: time.sleep(10000)

def CPlotDisplay__(filename):
    ext = os.path.splitext(filename)[-1].lower()
    extInfo = FILE_FORMATS.get(ext, None)
    if extInfo is None:
        if ext[:4].lower() == ".ref":
            fmt, isPyTree, dim = FILE_FORMATS[".pickle"]
        else:
            formats = ", ".join(sorted(FILE_FORMATS))
            raise ValueError(
                f"Format '{ext}' not supported."
                f"\nAvailable formats are:\n{formats}"
            )
    else:
        fmt, isPyTree, dim = extInfo

    if isPyTree:
        import Converter.PyTree as C
        a = C.convertFile2PyTree(filename, fmt)
        CPlot.PyTree.display(a, dim=dim)
    else:
        import Converter
        a = Converter.convertFile2Arrays(filename, fmt)
        CPlot.display(a, dim=dim)

    CPlot.setFileName__(filename)
    infLoop()

def main():
    if len(sys.argv) < 2:
        print("Usage: cplot <filename>")
        sys.exit(1)
    else:
        filename = sys.argv[1]
        CPlotDisplay__(filename)

if __name__ == "__main__":
    main()
