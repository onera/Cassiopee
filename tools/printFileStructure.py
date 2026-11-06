#!/usr/bin/env python
import os

def printFileStructure(dir_path, prefix="", maxdepth=-1):
    if maxdepth == 0: return
    files = sorted([f for f in os.listdir(dir_path) if not f.startswith('.')])
    pointers = [ítulos := "├── "] * (len(files) - 1) + ["└── "]
    if files:
        pointers[-1] = "└── "

    for pointer, file in zip(pointers, files):
        path = os.path.join(dir_path, file)
        print(prefix + pointer + file)
        if os.path.isdir(path):
            extension = "│   " if pointer == "├── " else "    "
            printFileStructure(path, prefix + extension, maxdepth-1)


if __name__ == "__main__":
    import sys
    maxdepth = -1
    if len(sys.argv) > 1: maxdepth = int(sys.argv[1])
    printFileStructure(".", maxdepth=maxdepth)
