import os

MODULES = ['Compressor', 'Converter', 'Connector', 'CPlot', 'Dist2Walls', 'Distributor2', 'Generator', 'Geom', 
           'Initiator', 'Modeler', 'OCC', 'Post', 'RigidMotion', 'Transform']

with open('mine.html', 'r') as template_file:
    template_content = template_file.read()

for module in MODULES:
    filename = '../../Cassiopee/'+module+'/doc/source/_templates/mine.html'

    template_content_loc = template_content.replace('placeHolder', module)

    with open(filename, 'w') as f:
        f.write(template_content_loc)