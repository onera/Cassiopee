# Compare session logs printed by validCassiopee
# Usage: python compareSessionLogs.py log1 log2
# Differences written in: compSession_DATE.txt where DATE is the current time
import os
import sys
from time import localtime, strftime

if len(sys.argv) != 3:
  raise Exception("Missing one or both input session logs")
  
filename1 = sys.argv[1]
filename2 = sys.argv[2]

def readLog(filename):
  if not os.access(filename, os.R_OK):
    raise Exception("Session log can't be read: {}".format(filename))
    
  session = []
  with open (filename, 'r') as f:
    for line in f:
      if ':' in line: 
        test = line.strip().split(':')
        session.append([e.strip() for e in test])
  return session

def diffTest(test, ref, new):
  refStatus = ref[5]
  newStatus = new[5]
  if refStatus != newStatus:
    return test, ref, new
  else:
    return ''

def stringify(test='', ref='', new=''):
  if not test:
    return test
  mod, test = test.split('/')
  test.split('.')[0]
  if not (ref or new):
    return "{:>15} | {:>42} |\n".format(mod, test)
  else:
    return "{:>15} | {:>42} | {:>10} | {:>10} |\n".format(mod, test, ref[5], new[5])

if __name__ == '__main__':
  # Read log files
  refSession = readLog(filename1)
  newSession = readLog(filename2)
  
  # Draw a one-to-one correspondance between them (module & testname)
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
  
  # Write differences to file
  compStr = ""
  header = "{} | {} | {} |\n{}\n".format("TESTS".center(60), "REF".center(10),
                                         "NEW".center(10), '*'*88)
  commonTestsHeader = "Common tests that differ:\n{}\n".format('-'*24)
  for test in commonTests:
    compStr += stringify(*diffTest(test, refDict[test], newDict[test]))
  if len(compStr): compStr = commonTestsHeader + compStr
  
  if len(newTests): compStr += "\nNew tests:\n{}\n".format('-'*9)
  for test in newTests:
    compStr += stringify(test)
  
  if len(deletedTests): compStr += "\nDeleted tests:\n{}\n".format('-'*13)
  for test in deletedTests:
    compStr += stringify(test)
    
  if compStr:
    now = strftime("%y%m%d-%H%M%S", localtime())
    with open("compSession_{}.txt".format(now), 'w') as f:
      f.write(header + compStr)
