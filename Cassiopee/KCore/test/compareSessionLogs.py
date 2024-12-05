# Compare session logs printed by validCassiopee
# Usage: python compareSessionLogs.py --logs='log1 log2'
#        python compareSessionLogs.py --logs='log1 log2' --email
# Differences written in: compSession_DATE.txt where DATE is the current time
import os
import sys
from time import strptime, strftime

# Parse command-line arguments
def parseArgs():
    import argparse
    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--email", action="store_true",
                        help="Email results. Default: print in terminal")
    parser.add_argument("-l", "--logs", type=str, default='',
                        help="Single-quoted logs to compare.")
    parser.add_argument("-p", "--prod", type=str, default='',
                        help="Name of the production.")
    parser.add_argument("-u", "--update", action="store_true",
                        help="Update valid. log on stck. Default: no update")
    # Parse arguments
    return parser.parse_args()

# Read git info in a session log of validCassiopee
def readGitInfo(filename):
    if not os.access(filename, os.R_OK):
        raise Exception("Session log can't be read: {}".format(filename))

    gitInfo = dict.fromkeys(['Base', 'Git origin', 'Git branch', 'Commit hash'])
    with open (filename, 'r') as f:
        for line in f:
            if 'Base from' in line:
                gitInfo['Base'] = line.strip().split(' ')[-1]
            if 'Git origin' in line:
                gitInfo['Git origin'] = line.strip().split(' ')[-1]
            if 'Git branch' in line:
                info = line.strip().split(' ')
                gitInfo['Git branch'] = info[2][:-1]
                gitInfo['Commit hash'] = info[-1]
            if all(v is not None for v in gitInfo.values()):
                break
    return gitInfo

# Find a two session logs of validCassiopee for a given production
def findLogs(prodname):
    import os
    from glob import glob
    validDataFolder = "/stck/cassiope/git/Cassiopee/Cassiopee/ValidData_{}".format(prodname)
    if not os.access(validDataFolder, os.R_OK):
        raise Exception("Session logs can't be retrieved in {}".format(validDataFolder))

    logs = None
    refLogs = sorted(glob(os.path.join(validDataFolder, "REF-session-*.log")))
    sessionLogs = sorted(glob(os.path.join(validDataFolder, "session-*.log")))
    if refLogs: logs = [refLogs[-1], sessionLogs[-1]]
    elif len(sessionLogs) > 1: logs = sessionLogs[-2:]
    return logs

# Read a session log of validCassiopee
def readLog(filename):
    if not os.access(filename, os.R_OK):
        raise Exception("Session log can't be read: {}".format(filename))

    session = []
    with open (filename, 'r') as f:
        for line in f:
            if ' :' in line:
                test = line.strip().split(':')
                session.append([e.strip() for e in test])
    return session

# Get the name of the ValidData directory containing the session log
def getLogDirName(filename):
    if not os.access(filename, os.R_OK):
        raise Exception("Session log can't be read: {}".format(filename))
    return os.path.basename(os.path.dirname(filename))

# Get the name of the prod from a session log
def getProd(filename):
    return os.path.basename(os.path.dirname(filename))[10:]

# Return time of creation of a session log of validCassiopee in two different
# formats (email subject vs write to file)
def getTimeFromLog(filename):
    if not os.access(filename, os.R_OK):
        raise Exception("Session log can't be read: {}".format(filename))
    time = filename[-17:-4]
    time1 = strptime(time, "%y%m%d_%H%M%S")
    return strftime("%d/%m/%y at %T", time1), time.replace("_", "-")

# Check the status of two tests using the test strings of validCassiopee
def diffTest(test, ref, new):
    refStatus = ref[5]
    newStatus = new[5]
    if refStatus != newStatus:
        return test, ref, new
    else:
        return ''

# Return test execution time of the Base, the reference log and the new log
def getExecTime(test, ref, new):
    def _testStr2Time(t):
        t = t[:-1].split('m')
        if len(t) == 1: return 1e-12 # No data
        return float(t[0])*60. + float(t[1])
    baseExecTime = _testStr2Time(ref[1])
    refExecTime = _testStr2Time(ref[0])
    newExecTime = _testStr2Time(new[0])
    return baseExecTime, refExecTime, newExecTime

# Return test execution time difference in % between (ref, new) and the Base
def getDiffExecTime(test, ref, new):
    cutoff = 1.0 # cut-off time in sec under which a test is not considered
    baseExecTime, refExecTime, newExecTime = getExecTime(test, ref, new)
    if baseExecTime < cutoff: return 0., 0.
    diffRef = round((refExecTime-baseExecTime)/baseExecTime*100., 1)
    diffNew = round((newExecTime-baseExecTime)/baseExecTime*100., 1)
    return diffRef, diffNew

# Stringify test comparison
def stringify(test='', ref='', new=''):
    if not test:
        return test
    mod, test = test.split('/')
    test.split('.')[0]
    if not (ref or new):
        return "{:>15} | {:>42} |\n".format(mod, test)
    elif not isinstance(ref, list):
        return "{:>15} | {:>42} | {:>10} | {:>10} |\n".format(mod, test, ref, new)
    else:
        return "{:>15} | {:>42} | {:>10} | {:>10} |\n".format(mod, test, ref[5], new[5])

# Main
if __name__ == '__main__':
    script_args = parseArgs()
    if script_args.prod:
        script_args.logs = findLogs(script_args.prod)
        if not(isinstance(script_args.logs, list) and len(script_args.logs) == 2):
            raise Exception("Two session logs were not found for prod. {}".format(
                script_args.prod))
    else:
        script_args.logs = script_args.logs.split(' ')
        if len(script_args.logs) != 2:
            raise Exception("Two session logs must be provided using the flag -l "
                            "or --logs")

    # Read log files and git info
    refSession = readLog(script_args.logs[0])
    newSession = readLog(script_args.logs[1])
    gitInfo = readGitInfo(script_args.logs[1])

    # Draw a one-to-one correspondance between tests of each session
    # (module + testname)
    refDict = dict((test[0] + '/' + test[1], test[2:]) for test in refSession)
    newDict = dict((test[0] + '/' + test[1], test[2:]) for test in newSession)
    refSet = set(refDict)
    newSet = set(newDict)

    # Find common tests
    commonTests = sorted(refSet & newSet)
    # Find new tests (tests in newSession but not in refSession)
    newTests = sorted(newSet - refSet)
    # Find deleted tests (tests in refSession but not in newSession)
    deletedTests = sorted(refSet - newSet)
    # Find failed tests in newSession
    failedTests = sorted([k for k, v in newDict.items() if v[5] != 'OK'])

    # Write differences to terminal or send an email
    baseState = 'OK'
    compStr = ""
    header = "\n".join("{}: {}".format(k,v) for k,v in gitInfo.items())
    header += "\n\nREF = {}\n".format(script_args.logs[0])
    header += "NEW = {}\n\n\n".format(script_args.logs[1])
    header += "{} | {} | {} |\n{}\n".format("TESTS".center(60), "REF".center(10),
                                            "NEW".center(10), '*'*88)
    commonTestsHeader = "Tests that differ:\n{}\n".format('-'*17)
    for test in commonTests:
        compStr += stringify(*diffTest(test, refDict[test], newDict[test]))
    if len(compStr):
        compStr = commonTestsHeader + compStr
        baseState = 'FAILED'
    else: compStr = commonTestsHeader + "[none]\n"

    newTestsHeader = "\nNew tests:\n{}\n".format('-'*9)
    if len(newTests):
        compStr += newTestsHeader
        for test in newTests:
            compStr += stringify(test)
        if baseState == 'OK': baseState = 'ADDITIONS'
    else: compStr += newTestsHeader + "[none]\n"

    deletedTestsHeader = "\nDeleted tests:\n{}\n".format('-'*13)
    if len(deletedTests):
        compStr += deletedTestsHeader
        for test in deletedTests:
            compStr += stringify(test)
        if baseState == 'OK': baseState = 'DELETIONS'
        elif baseState == 'ADDITIONS': baseState += ' & DELETIONS'
    else: compStr += deletedTestsHeader + "[none]\n"

    failedTestsHeader = "\nReminder - Failed tests:\n{}\n".format('-'*23)
    if len(failedTests):
        compStr += failedTestsHeader
        for test in failedTests:
            compStr += stringify(test)
        baseState = 'FAILED'
    else: compStr += failedTestsHeader + "[none]\n"

    execTime = []
    for test in commonTests:
        execTime.append([test, *getDiffExecTime(test, refDict[test], newDict[test])])
    execTime.sort(key=lambda x: x[2])

    if script_args.email:
        threshold = 50.
        execTimeHeader = "\nExecution time - {:.0f}% threshold:\n{}\n".format(threshold, '-'*30)
        compStr += execTimeHeader
        compStr += "{} | {} | {} |\n{}\n".format(" "*60, "REF v Base".center(10),
                                                 "NEW v Base".center(10), '*'*88)
        cmpt = 0
        for test in execTime:
            if abs(test[2]) > threshold and abs(test[1]) < threshold:
                compStr += stringify(test[0], str(test[1])+'%', str(test[2])+'%')
                cmpt += 1
        if cmpt == 0: compStr += "[none]\n"

    baseStateMsg = ""
    prod = getProd(script_args.logs[1])
    tlog, tlog2 = getTimeFromLog(script_args.logs[1])
    messageSubject = "[validCassiopee - {}] {} - State: {}".format(prod, tlog, baseState)
    messageText = header + compStr + baseStateMsg
    if script_args.email:
        if baseStateMsg: baseStateMsg = '\n\n'+baseStateMsg
        try:
            from KCore.notify import notify
            notify(recipients=['vincent.casseau@onera.fr', 'christophe.benoit@onera.fr'],
                   messageSubject=messageSubject,
                   messageText=messageText)
        except ImportError:
            raise SystemError("Error: KCore is required to import notify.")
    else:
        print("{0}\n|{1:^86}|\n{0}\n{2}".format(88*'-', messageSubject, messageText))

    if script_args.update:
        # If the state of the Base is OK, set the new session log to be the
        # reference
        exitStatus = 0
        if (any(st in baseState for st in ['OK', 'ADDITIONS', 'DELETIONS']) and
                os.path.basename(script_args.logs[0]).startswith('REF-')):
            if os.access(script_args.logs[0], os.W_OK):
                import shutil
                os.remove(script_args.logs[0])
                newRef = os.path.join(os.path.dirname(script_args.logs[1]),
                                      'REF-' + os.path.basename(script_args.logs[1]))
                shutil.copyfile(script_args.logs[1], newRef)
            else: exitStatus = 2

        # Amend state of the base in logs/validation_status.txt
        logAllValids = "/stck/cassiope/git/logs/validation_status.txt"
        entry = "{} - {} - {}\n".format(prod, tlog2, baseState)
        if os.access(os.path.dirname(logAllValids), os.W_OK):
            with open(logAllValids, 'r') as f: contents = f.readlines()
            prodFound = False
            for i, line in enumerate(contents):
                if line.startswith(prod):
                    contents[i] = entry
                    prodFound = True
                    break

            if not prodFound: contents.append(entry)
            with open(logAllValids, 'w') as f: f.writelines(contents)

        sys.exit(exitStatus)
