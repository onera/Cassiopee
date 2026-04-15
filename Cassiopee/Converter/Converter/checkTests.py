# check unitary tests
# for every cassiopee functions
import os, ast, subprocess

MODULEBLACKLIST = []
FILEBLACKLIST = ['__init__', 'cgnslib', 'cgnserrors', 'cgnskeywords', 'cgnstypes',
                 'cgnsutils', 'vera', 'text1', 'chancery', 'courier', 'nimbus']
FUNCTIONBLACKLIST = ['send']

def getCassiopeeSourceDir():
    try:
        # Check installPath
        import KCore.installPath
        cassiopeeIncDir = KCore.installPath.includePath
        cassiopeeIncDir = os.path.dirname(cassiopeeIncDir)
    except ImportError:
        raise SystemError("Error: KCore module is required to use this script.")
    return cassiopeeIncDir

def getModules():
    cassiopeeIncDir = getCassiopeeSourceDir()
    try: mods = os.listdir(cassiopeeIncDir)
    except: mods = []
    out = []
    for mod in mods:
        a = os.access(os.path.join(cassiopeeIncDir, mod, 'test'), os.F_OK)
        if a:
            if mod not in MODULEBLACKLIST:
                out.append(mod)
    return out

# get the functions of module
# return a dict for each module python file containing
# (function name, function doc string)
def getFunctions(module):
    cassiopeeIncDir = getCassiopeeSourceDir()
    path = os.path.join(cassiopeeIncDir, module, module)
    files = os.listdir(path)
    modFuncs = {}
    for f in files: # python files in module

        if not f.endswith(".py"): continue

        name = f.replace(".py", "")

        # filter by file name
        if name in FILEBLACKLIST: continue

        modFuncs[name] = []

        file = os.path.join(path, f)
        with open(file, "r", encoding="utf-8") as f:
            tree = ast.parse(f.read())

        functions = []
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                functions.append( (node.name, ast.get_docstring(node)) ) 
            if isinstance(node, ast.ClassDef):
                for sub in node.body:
                    if isinstance(sub, ast.FunctionDef):
                        functions.append( (f"{node.name}.{sub.name}", ast.get_docstring(sub)) )

        # filter by function names
        for fu in functions:
            if fu[0].endswith('__'): continue
            if fu[0][0] == '_': continue 
            modFuncs[name].append(fu)

    return modFuncs

# IN: runTest: if true, run the tests of functionName
# IN: runHeader: if true, check test headers
# IN: runMissing: if true, check if tests are missing for functionName
# IN: runDoc: if true, check doc string and rst for functionName
def checkUnitaryTests(module, functionName, docString, runTest=False, runHeader=False, 
                      runMissing=False, runDoc=False):

    if functionName in FUNCTIONBLACKLIST: return 0

    ret = 0
    cassiopeeIncDir = getCassiopeeSourceDir()
    path = os.path.join(cassiopeeIncDir, module, 'test')
    
    # check doc string
    if runDoc:
        if docString is None: 
            print("Warning: %s docstring is missing."%(functionName))
            ret += 1

    # check name.py (doc test)
    file = os.path.join(path, functionName+'.py')
    a = os.path.exists(file)
    if not a:
        if runMissing:
            print("Warning: file %s is missing."%(functionName+'.py'))
            ret += 1
    else:
        # check test
        if runTest:
            try:
                subprocess.check_output(["python", file], stderr=subprocess.STDOUT, text=True)
            except Exception:
                print("Error: test %s fails."%file)
                ret += 1

        # check header
        if runHeader:
            with open(file, "r", encoding="utf-8") as f:
                l = f.readline()
                ref = "# - %s (array) -"%functionName
                l = l[0:len(ref)]
                if l != ref:
                    print("Warning: bad header in file: %s."%file)
                    print("%s %s"%(ref, l))
                    ret += 1

    # check namePT.py (doc test)
    file = os.path.join(path, functionName+'PT.py')
    a = os.path.exists(file)
    if not a:
        if runMissing:
            print("Warning: file %s is missing."%(functionName+'PT.py'))
            ret += 1
    else:
        # check test
        if runTest:
            try:
                subprocess.check_output(["python", file], stderr=subprocess.STDOUT, text=True)
            except Exception:
                print("Error: test %s fails."%file)
                ret += 1
        # check header
        if runHeader:
            with open(file, "r", encoding="utf-8") as f:
                l = f.readline()
                ref = "# - %s (pyTree) -"%functionName
                l = l[0:len(ref)]
                if l != ref:
                    print("Warning: bad header in file: %s."%file)
                    print("%s %s"%(ref, l))
                    ret += 1

    # check unitary test (array)
    file = os.path.join(path, functionName+'_t1.py')
    a = os.path.exists(file)
    if not a:
        if runMissing:
            print("Warning: file %s is missing."%(functionName+'_t1.py'))
            ret += 1
    else:
        # check header
        if runHeader:
            with open(file, "r", encoding="utf-8") as f:
                l = f.readline()
                ref = "# - %s (array) -"%functionName
                l = l[0:len(ref)]
                if l != ref:
                    print("Warning: bad header in file: %s."%file)
                    print("%s %s"%(ref, l))
                    ret += 1

    # check unitary test (pyTree)
    file = os.path.join(path, functionName+'PT_t1.py')
    a = os.path.exists(file)
    if not a:
        if runMissing:
            print("Warning: file %s is missing."%(functionName+'PT_t1.py'))
            ret += 1
    else:
        # check header
        if runHeader:
            with open(file, "r", encoding="utf-8") as f:
                l = f.readline()
                ref = "# - %s (pyTree) -"%functionName
                l = l[0:len(ref)]
                if l != ref:
                    print("Warning: bad header in file: %s."%file)
                    print("%s %s"%(ref, l))
                    ret += 1
    return ret

#==============================================================================
if __name__ == "__main__":
    modules = getModules()
    errorTot = 0
    for module in modules: # for each module
        errors = 0
        funcs = getFunctions(module)
        fileNames = funcs.keys()
        for file in fileNames: # for each python file
            functionNames = funcs[file]
            for functionName in functionNames:
<<<<<<< Updated upstream
                e = checkUnitaryTests(module, functionName,
                                      runMissing=False, runHeader=True, runTest=False)
=======
                e = checkUnitaryTests(module, functionName[0], functionName[1], 
                                      runMissing=False, runHeader=False, runTest=False, runDoc=True)
>>>>>>> Stashed changes
                errors += e
        if errors > 0:
            print("%s has %d python files."%(module, len(fileNames)))
            for file in fileNames:
                print("     %s/%s has %d functions."%(module, file, len(functionNames)))
            print("%s has %d errors."%(module, errors))
        errorTot += errors
    print("Total of %s errors."%errorTot)