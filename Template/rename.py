# Adapt the module name
import os, sys
if len(sys.argv) < 2:
    print('python rename.py <Module>')
    print('rename template module as <Module>.')
    import sys; sys.exit()

MODULE = sys.argv[1]
print('Renaming module as '+MODULE)

st1 = MODULE.capitalize()
st2 = MODULE.lower()
st3 = MODULE.upper()

w = os.walk('.')
for f in w:
    # exclude .svn, build
    dirName = f[0]
    if dirName.find('build') == -1 and dirName.find('.svn') == -1:
        files = f[2]
        for f in files:
            d = open(dirName+'/'+f)
            lines = d.read()
            d.close()
            lines = lines.replace('Template', st1)
            lines = lines.replace('template', st2)
            lines = lines.replace('TEMPLATE', st3)
            d = open(dirName+'/'+f, "w")
            d.write(lines)
            d.close()

# - Some file operations -
# Python file
os.rename('Template/Template.py', 'Template/'+st1+'.py')
# Module file
os.rename('Template/template.cpp', 'Template/'+st2+'.cpp')
# Module include
os.rename('Template/template.h', 'Template/'+st2+'.h')
# Source dir
os.rename('Template', st1)
