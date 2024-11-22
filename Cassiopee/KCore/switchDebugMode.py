# Switch debug mode in KCore/Dist.py, ie, toggles DEBUG
# Usage: python switchDebugMode.py --activate
#        python switchDebugMode.py --deactivate
import os
import sys
import argparse

# Parse command-line arguments
def parseArgs():
    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--activate", action="store_true",
                        help="Activate debug mode.")
    # Parse arguments
    return parser.parse_args()

# Get installation path of Cassiopee
def getInstallPath():
    try:
        try: import KCore.installPath as installPath
        except: import installPath
        cassiopeeIncDir = installPath.includePath
        cassiopeeIncDir = os.path.dirname(cassiopeeIncDir)
        return cassiopeeIncDir
    except ImportError:
        print("Warning: KCore library is required to use this script. Skipping...")
        sys.exit()

# Read KCore/Dist
def readDist():
    cassiopeeIncDir = getInstallPath()
    filename = os.path.join(cassiopeeIncDir, "KCore/Dist.py")
    if not os.access(filename, os.R_OK):
        raise Exception("Dist.py can't be read at: {}".format(filename))
    contents = []
    with open (filename, 'r') as f:
        for line in f: contents.append(line)
    return contents

# Edit Dist
def editDist(contents, dbgMode):
    for i, line in enumerate(contents):
        if "DEBUG =" in line:
            contents[i] = "DEBUG = " + dbgMode + '\n'
            return
    return

# Write KCore/Dist
def writeDist(contents):
    cassiopeeIncDir = getInstallPath()
    filename = os.path.join(cassiopeeIncDir, "KCore/Dist.py")
    if not os.access(filename, os.W_OK):
        raise Exception("Modified Dist.py can't be written at: {}".format(filename))
    with open (filename, 'w') as f:
        for line in contents: f.write(line)
    return

# Check ELSAPROD
def check_elsaprod(dbgMode):
    elsaprod = os.getenv("ELSAPROD")
    if elsaprod is not None:
        if dbgMode and not "_DBG" in elsaprod:
            print("Add '_DBG' suffix to $ELSAPROD: {}_DBG".format(elsaprod))
        elif not dbgMode and "_DBG" in elsaprod:
            print("Remove '_DBG' suffix from $ELSAPROD: {}".format(elsaprod[:-4]))
    return

if __name__ == '__main__':
    args = parseArgs()
    dbgMode = "True" if args.activate else "False"
    contents = readDist()
    editDist(contents, dbgMode)
    writeDist(contents)
    check_elsaprod(args.activate)
