# Remplace le numero de version dans les modules
# usage: python version.py 3.5 3.6
import os, sys
CASSIOPEE = os.environ['CASSIOPEE']

# Recupere la liste des modules definie dans le shell MODULES
def getModuleList(path):
    f = open(path+'/MODULES')
    a = f.read()
    f.close()
    a = a.replace('FREEMODULES', '')
    a = a.replace('FULLMODULES', '')
    a = a.replace('OTHERS', '')
    a = a.replace('export', '')
    a = a.replace('=', '')
    a = a.replace("'", '')
    a = a.replace("$", '')
    a = a.split(' ')
    out = []
    for i in a:
        i = i.split('\n')
        for r in i:
            if r != ' ' and r != '': out.append(r)
    return out

def update(path, version1, version2):
    MODULES = getModuleList(path)
    for m in MODULES:
      print("Updating %s."%m)
      w = os.walk(path+'/'+m)
      for f in w:
          # exclude .svn, .git,  build
          dirName = f[0]
          if dirName.find('build') == -1 and dirName.find('.svn') == -1 and dirName.find('.git') == -1:
              files = f[2]
              for fi in files:
                  if fi == 'setup.py' or fi == 'setupLegacy.py' or fi == m+'.py' or fi == '__init__.py' or fi == 'version.py':
                      d = open(dirName+'/'+fi, "rb")
                      lines = d.read()
                      d.close()
                      lines = lines.replace(version1, version2)
                      d = open(dirName+'/'+fi, "wb")
                      d.write(lines)
                      d.close()

#==============================================================================
if __name__ == "__main__":

  if len(sys.argv) != 3: 
      print('Usage: version 3.0 3.1\nChange version 3.0 in 3.1.')
      sys.exit()
  version1 = sys.argv[1]
  version2 = sys.argv[2]
  
  print('Changing '+version1+' -> '+version2)

  version1 = bytes(version1, encoding='utf-8')
  version2 = bytes(version2, encoding='utf-8')
  
  print('in Cassiopee...')
  update(CASSIOPEE+'/Cassiopee/', version1, version2)
  
  
