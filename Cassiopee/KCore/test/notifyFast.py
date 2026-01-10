# Notify a user about:
#   - the installation status
#   - the checkout status
#   - the validation status
# of different prods, either in the terminal window (default) or by email
# (set the environment variable CASSIOPEE_EMAIL and run on localhost)
# Usage: python notifyFast.py --install
#        python notifyFast.py --install --email --recipients='a.b@onera.fr c.d@onera.fr'
#        python notifyFast.py --checkout
#        python notifyFast.py --valid
#        python notifyFast.py --valid --prod=<prod_name>
#        python notifyFast.py --valid --prod=<prod_name> --full
import os
import sys
from glob import glob
from time import strptime, strftime

try:
    import KCore.Dist as Dist
    from KCore.notify import notify
except ImportError:
    print("Error: KCore is required to execute notifyFast.py")
    sys.exit()


# Tests to ignore in non-debug mode
IGNORE_TESTS_NDBG = []
# Tests to ignore in debug mode
IGNORE_TESTS_DBG = []


# Parse command-line arguments
def parseArgs():
    import argparse
    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--checkout", action="store_true",
                        help="Show checkout status. Default: disabled")
    parser.add_argument("-e", "--email", action="store_true",
                        help="Email results. Default: print in terminal")
    parser.add_argument("-i", "--install", action="store_true",
                        help="Show installation status. Default: disabled")
    parser.add_argument("-f", "--full", action="store_true",
                        help="Show test case logs. Default: disabled")
    parser.add_argument("-l", "--logs", type=str, default='',
                        help="Single-quoted logs to compare.")
    parser.add_argument("-p", "--prod", type=str, default='',
                        help="Name of the production.")
    parser.add_argument("-r", "--recipients", type=str, default='',
                        help="Single-quoted space-separated list of recipients")
    parser.add_argument("-u", "--update", action="store_true",
                        help="Update valid. log on stck. Default: disabled")
    parser.add_argument("-v", "--valid", action="store_true",
                        help="Show validation status. Default: disabled")
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
def findLogs(prodname, findRef=True):
    validDataFolder = "/stck/cassiope/git/Cassiopee/Cassiopee/ValidData_{}".format(prodname)
    if not os.access(validDataFolder, os.R_OK):
        raise Exception("Session logs can't be retrieved in {}".format(validDataFolder))

    logs = None
    refLogs = []
    if findRef: refLogs = sorted(glob(os.path.join(validDataFolder, "REF-sessionFast-*.log")))
    sessionLogs = sorted(glob(os.path.join(validDataFolder, "sessionFast-*.log")))
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

# Return a list of tests to ignore for a given prod.
def tests2Ignore(prod):
    if '_DBG' in prod:
        return IGNORE_TESTS_DBG
    return IGNORE_TESTS_NDBG

# Get full test logs for a given list of tests
def getTestLogs(prodname, testList):
    testLog = ""
    modNames = [test.split('/')[0] for test in testList]
    testNames = [test.split('/')[1] for test in testList]
    validDataFolder = "/stck/cassiope/git/Cassiopee/Cassiopee/ValidData_{}".format(prodname)
    # Read the last purged logValiFast.dat file
    purgedLogs = sorted(glob(os.path.join(validDataFolder, "logValidFast_purged_*.dat")))
    if not purgedLogs or not os.access(purgedLogs[-1], os.R_OK): return testLog
    with open(purgedLogs[-1], 'r') as f: log = f.read()
    # Split log using "Running " as the delimiter
    failedTests = log.split("Running ")[1:]
    failedTestNames = [test.split(' ')[0] for test in failedTests]
    # Add logs of cases that failed
    for i, name in enumerate(testNames):
        try:
            pos = failedTestNames.index(name)
            testLog += "{}/{}{}\n\n".format(modNames[i], failedTests[pos], 88*'*')
        except ValueError: pass
    return testLog

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

# Check install status
def checkInstallStatus():
    log_entries = []
    with open('/stck/cassiope/git/logs/installation_status.txt', 'r') as f:
        for line in f:
            log_entries.append(line.strip().split(' - '))
    log_entries.sort(key=lambda x: x[3], reverse=True)

    # Get git info
    fastIncDir = '/stck/cassiope/git/Fast/Fast'
    gitOrigin = Dist.getGitOrigin(fastIncDir)
    gitInfo = "Git origin: {}".format(gitOrigin)

    baseState = 'OK'
    messageText = "Installation of Cassiopee, Fast and all "\
        "PModules:\n{}\n\n{}\n\n".format(48*'-', gitInfo)
    messageText += '{:^22} | {:^6} | {:^7} | {:^24} | {:^10}\n{}\n'.format(
        "PROD.", "BRANCH", "HASH", "DATE", "STATUS", 83*'-')
    for log_machine in log_entries:
        prod = log_machine[0]
        gitBranch = log_machine[1]
        gitHash = log_machine[2]
        date = strptime(log_machine[3], "%y%m%d-%H%M%S")
        date = strftime("%d/%m/%y at %T", date)
        status = log_machine[4]
        messageText += '  {:<20} | {:^6} | {:^7} | {:^24} | {:^10}\n'.format(
            prod, gitBranch, gitHash, date, status)
        if 'FAILED' in log_machine: baseState = 'FAILED'

    messageSubject = "[Install Fast] State: {}".format(baseState)
    if baseState == 'FAILED':
        messageText += '\n\nIf the prod. you wish to use is marked as FAILED, '\
            'please contact the maintainer:\nivan.mary@onera.fr'

    return messageSubject, messageText

# Check checkout status
def checkCheckoutStatus(sendEmail=False):
    log_entries = []
    with open('/stck/cassiope/git/logs/checkout_status.txt', 'r') as f:
        for line in f:
            log_entries.append(line.strip().split(' - '))
    log_entries.sort(key=lambda x: x[1], reverse=True)

    # Do not send a notification when everything is OK
    if not any('FAILED' in log_machine for log_machine in log_entries):
        if sendEmail: sys.exit(0)
        else:
            messageSubject = "[Checkout Fast] State: OK"
            messageText = ""
            return messageSubject, messageText

    # Get git info
    fastIncDir = '/stck/cassiope/git/Fast/Fast'
    gitOrigin = Dist.getGitOrigin(fastIncDir)
    gitBranch = Dist.getGitBranch(fastIncDir)
    gitHash = Dist.getGitHash(fastIncDir)[:7]
    gitInfo = "Git origin: {}\nGit branch: {}\nCommit hash: {}".format(
        gitOrigin, gitBranch, gitHash)

    messageSubject = "[Checkout Fast] State: FAILED"
    messageText = "Pulling updates for Cassiopee, Fast and all "\
        "PModules:\n{}\n\n{}\n\n".format(52*'-', gitInfo)
    messageText += '{:^20} | {:^15} | {:^30} | {:^10}\n{}\n'.format(
        "PROD.", "PCKG.", "DATE", "STATUS", 85*'-')
    for log_machine in log_entries:
        prod = log_machine[0]
        pckg = log_machine[1]
        date = strptime(log_machine[2], "%y%m%d-%H%M%S")
        date = strftime("%d/%m/%y at %T", date)
        status = log_machine[3]
        messageText += '{:^20} | {:^15} | {:^30} | {:^10}\n'.format(
            prod, pckg, date, status)

    messageText += '\n\nIf the prod. you wish to use is marked as FAILED, '\
        'please contact the maintainer:\nivan.mary@onera.fr'

    return messageSubject, messageText

# Check valid status
def checkValidStatus():
    log_entries = []
    with open('/stck/cassiope/git/logs/validationFast_status.txt', 'r') as f:
        for line in f:
            log_entries.append(line.strip().split(' - '))
    log_entries.sort(key=lambda x: x[3], reverse=True)

    # Get git info
    fastIncDir = '/stck/cassiope/git/Fast/Fast'
    gitOrigin = Dist.getGitOrigin(fastIncDir)
    gitInfo = "Git origin: {}".format(gitOrigin)

    vnvState = 'OK'
    messageText = "Non-regression testing of Fast:\n{}\n\n{}\n\n".format(30*'-', gitInfo)
    messageText += '{:^22} | {:^6} | {:^7} | {:^24} | {:^10}\n{}\n'.format(
        "PROD.", "BRANCH", "HASH", "DATE", "STATUS", 83*'-')
    for log_machine in log_entries:
        prod = log_machine[0]
        gitBranch = log_machine[1]
        gitHash = log_machine[2]
        date = strptime(log_machine[3], "%y%m%d-%H%M%S")
        date = strftime("%d/%m/%y at %T", date)
        status = log_machine[4]
        messageText += '  {:<20} | {:^6} | {:^7} | {:^24} | {:^10}\n'.format(
            prod, gitBranch, gitHash, date, status)
        if 'FAILED' in log_machine: vnvState = 'FAILED'

    messageSubject = "[V&V Fast] State: {}".format(vnvState)
    if vnvState == 'FAILED':
        messageText += '\n\nIf the prod. you wish to use is marked as FAILED, '\
            'please contact the maintainer:\nivan.mary@onera.fr\n'\
            'or list remaining issues with:\n'\
            'notifyFast --valid --prod=your_prod_name --full'

    return messageSubject, messageText

# Compare session logs
def compareSessionLogs(logFiles=[], showExecTimeDiffs=False,
                       showTestLogs=False, update=False):
    # Read log files and git info
    refSession = readLog(logFiles[0])
    newSession = readLog(logFiles[1])
    gitInfo = readGitInfo(logFiles[1])

    # Get prod name
    prod = getProd(logFiles[1])

    # Delete tests that were not run
    refSession = [test for test in refSession if test[-1] != '...']
    newSession = [test for test in newSession if test[-1] != '...']

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
    failedTests = [t for t in failedTests if t not in tests2Ignore(prod)]

    # Write differences to terminal or send an email
    baseState = 'OK'
    compStr = ""
    header = "\n".join("{}: {}".format(k,v) for k,v in gitInfo.items())
    header += "\n\nREF = {}\n".format(logFiles[0])
    header += "NEW = {}\n\n\n".format(logFiles[1])
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

    if showExecTimeDiffs:
        threshold = 50.
        execTimeHeader = "\nExecution time - {:.0f}% threshold:\n{}\n".format(
            threshold, '-'*30)
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
    tlog, tlog2 = getTimeFromLog(logFiles[1])
    messageSubject = "[validFast - {}] {} - State: {}".format(prod, tlog,
                                                              baseState)
    messageText = header + compStr + baseStateMsg

    if showTestLogs:
        testLogs = getTestLogs(prod, failedTests)
        if testLogs:
            messageText += f"\n\nFailed test logs:\n{'-'*16}\n{testLogs}"

    exitStatus = 0
    if update:
        # If the state of the Base is OK, set the new session log to be the
        # reference
        if (any(st in baseState for st in ['OK', 'ADDITIONS', 'DELETIONS']) and
                os.path.basename(logFiles[0]).startswith('REF-')):
            if os.access(logFiles[0], os.W_OK):
                import shutil
                os.remove(logFiles[0])
                newRef = os.path.join(os.path.dirname(logFiles[1]),
                                      'REF-' + os.path.basename(logFiles[1]))
                shutil.copyfile(logFiles[1], newRef)
            else: exitStatus = 2
        else: exitStatus = 1

        # Amend state of the base in logs/validationFast_status.txt
        logAllValids = "/stck/cassiope/git/logs/validationFast_status.txt"
        entry = "{} - {} - {} - {} - {}\n".format(prod, gitInfo['Git branch'],
                                                  gitInfo['Commit hash'], tlog2,
                                                  baseState)
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

    return messageSubject, messageText, exitStatus

# Main
if __name__ == '__main__':
    exitStatus = 0
    scriptArgs = parseArgs()
    recipients = scriptArgs.recipients.split(' ')
    if not recipients[0]: recipients = []
    if not (scriptArgs.install or scriptArgs.checkout or scriptArgs.valid):
        scriptArgs.install = True  # Show installation status by default

    if scriptArgs.install:
        messageSubject, messageText = checkInstallStatus()
    if scriptArgs.checkout:
        messageSubject, messageText = checkCheckoutStatus(sendEmail=scriptArgs.email)
    elif scriptArgs.valid:
        mode = "overview"
        if scriptArgs.prod:
            findRef = False if scriptArgs.logs == "latest" else True
            scriptArgs.logs = findLogs(scriptArgs.prod, findRef=findRef)
            if not(
                isinstance(scriptArgs.logs, list) and
                len(scriptArgs.logs) == 2
            ):
                raise Exception("Two session logs were not found for "
                                "prod. {}".format(scriptArgs.prod))
            mode = "compare"
        elif scriptArgs.logs:
            scriptArgs.logs = scriptArgs.logs.split(' ')
            if len(scriptArgs.logs) != 2:
                raise Exception("Two session logs must be provided using the "
                                "flag -l or --logs")
            mode = "compare"

        if mode == "overview":
            messageSubject, messageText = checkValidStatus()
        else:
            messageSubject, messageText, exitStatus = compareSessionLogs(
                logFiles=scriptArgs.logs,
                showExecTimeDiffs=scriptArgs.email,
                showTestLogs=scriptArgs.full,
                update=scriptArgs.update
            )

    if scriptArgs.email:
        notify(recipients=recipients,
               messageSubject=messageSubject,
               messageText=messageText)
    else:
        sep = 83*'-'
        print(f"{sep}\n|{messageSubject:^81}|\n{sep}\n{messageText}")
    sys.exit(exitStatus)
