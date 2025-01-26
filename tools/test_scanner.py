# Scan Cassiopee's MODULES and output in a log file all (non-internal) functions
# that are not listed in any tests of the module they were defined in.
# If a function is also available inplace, only one of these two functions needs
# to have a test.
#
# Usage:
#  - run: python test_scanner.py
#  - open the output file ./missingTests.txt
 
import os
import sys
import re

CASSIOPEE = os.path.join(os.environ.get("CASSIOPEE"), 'Cassiopee')

def findFileNames(rootDir, targets=[], patterns=[], excludeDirs=[]):
    """
    Search for all files with the specified name in any subdirectory

    :param rootDir: The root directory to start searching from
    :param fileNames: The name of the file to search for
    :return: A list of file paths matching any of the target file names
    """
    matches = []
    regexs = [re.compile(p) for p in patterns]
    
    for dirpath, _, filenames in os.walk(rootDir):
        if any(d in dirpath for d in excludeDirs): continue
        for target in targets:
            if target in filenames:
                matches.append(os.path.join(dirpath, target))
        for filename in filenames:
            for regex in regexs:
                if regex.match(filename):
                    matches.append(os.path.join(dirpath, filename))
    return matches
    
def findFunctionNames(fileNames, excludePrivate=True):
    def _find(fileName):
        names = []
        with open(fileName, 'r') as f:
            for line in f:
                if "def " in line:
                    name = line.strip()
                    if name[0] == '#': continue
                    name = name.split()[1].split('(')[0]
                    if excludePrivate and '__' in name: continue
                    names.append(name)
        return names
    return {fileName: _find(fileName) for fileName in fileNames}
    
def arePatternsFound(fileNames, patterns, include_inline=True):
    def _isPatternFound(fileNames, pattern):
        for fileName in fileNames:
            with open(fileName, 'r') as f:
                for line in f:
                    if pattern in line:
                        return True
                    elif include_inline and pattern[0] == "_" and pattern[1:] in line:
                        return True    
        return False
    return [_isPatternFound(fileNames, pattern) for pattern in patterns]


if __name__ == "__main__":
    # Read module names
    with open(os.path.join(CASSIOPEE, 'MODULES'), 'r') as f:
        moduleNames = f.readline().split("'")[1].split()
        
    foundOutfile = open('foundTests.txt', 'w')
    missingOutfile = open('missingTests.txt', 'w')
    
    # Loop over all modules
    for module in moduleNames:
        print(f">>> Module {module}")
        
        # Search for python source files
        rootDir = os.path.join(CASSIOPEE, module)
        fileNames = ["PyTree.py", f"{module}.py", "Internal.py", "Mpi.py"]
        pysrc = findFileNames(
            rootDir=rootDir,
            targets=fileNames,
            excludeDirs=['build', 'test', 'mystuff']
        )
        if not pysrc: continue
        #print("pysrc", pysrc)
        
        fnNames = findFunctionNames(pysrc)
        #print("fnNames", fnNames)
        
        # Search for test cases
        fileNames = [r".*_[tm][0-9]*\.py"]
        tests = findFileNames(
            rootDir=os.path.join(rootDir, 'test'),
            patterns=fileNames
        )
        if not tests: continue
        #print("tests", tests)
        
        foundFnNames, missingFnNames = {}, {}
        for src, fns in fnNames.items():
            testFound = arePatternsFound(tests, fns)
            foundFnNames[src] = [fn for i, fn in enumerate(fns) if testFound[i]]
            missingFnNames[src] = [fn for i, fn in enumerate(fns) if not testFound[i]]
        
        #print("fnNames", fnNames)

        # Print to file
        missingOutfile.write(f"Module {module}\n")
        for src, fns in missingFnNames.items():
            if not missingFnNames[src]: continue
            source = os.path.relpath(src, CASSIOPEE)
            missingOutfile.write(f"    - Source {source}\n")
            missingOutfile.write("\n".join(f"        + {fn}" for fn in fns))
            missingOutfile.write("\n\n")
        missingOutfile.write('-'*80 + "\n\n")
        
        foundOutfile.write(f"Module {module}\n")
        for src, fns in foundFnNames.items():
            if not foundFnNames[src]: continue
            source = os.path.relpath(src, CASSIOPEE)
            foundOutfile.write(f"    - Source {source}\n")
            foundOutfile.write("\n".join(f"        + {fn}" for fn in fns))
            foundOutfile.write("\n\n")
        foundOutfile.write('-'*80 + "\n\n")

    missingOutfile.close()
    foundOutfile.close()
