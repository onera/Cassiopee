# Switch size of Int in KCore/Dist.py, ie, toggles EDOUBLEINT and GDOUBLEINT
# Usage: python switchIntSize.py -i=4
#        python switchIntSize.py --int=8
import os
import sys
import argparse

# Parse command-line arguments
def parseArgs():
    def _checkInt(i):
        def _throwError():
            raise argparse.ArgumentTypeError("Int size must be 4 or 8.")
            sys.exit()
        try: i = int(i)
        except: _throwError()
        if i in [4, 8]: return i
        elif i == 32: return 4
        elif i == 64: return 8
        else: _throwError()

    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--int", type=_checkInt,
                        help="Size of Int: 4 or 8. Default: 4.")
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
def editDist(contents, intState):
    sw1 = True
    sw2 = True
    for i, line in enumerate(contents):
        if sw1 and "EDOUBLEINT =" in line:
            contents[i] = "EDOUBLEINT = " + intState + '\n'
            sw1 = False
        if "GDOUBLEINT =" in line:
            contents[i] = "GDOUBLEINT = " + intState + '\n'
            sw2 = False
        if not (sw1 or sw2): return
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
def check_elsaprod(intSize):
    elsaprod = os.getenv("ELSAPROD")
    if elsaprod is not None:
        if intSize == 4 and "_i8" in elsaprod:
            print("Remove '_i8' suffix from $ELSAPROD: {}".format(
                elsaprod.replace('_i8', '')))
        elif intSize == 8 and "_i8" not in elsaprod:
            if '_DBG' in elsaprod:
                print("Add '_i8_DBG' suffix to $ELSAPROD: {}_i8_DBG".format(
                    elsaprod[:-4]))
            else:
                print("Add '_i8' suffix to $ELSAPROD: {}_i8".format(elsaprod))
    return

if __name__ == '__main__':
    args = parseArgs()
    intState = "True" if args.int == 8 else "False"
    contents = readDist()
    editDist(contents, intState)
    writeDist(contents)
    check_elsaprod(args.int)
