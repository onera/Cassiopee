# - tkLuxRender -
# Interface avec Lux render
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.Internal as Internal
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Generator.PyTree as G
import KCore.Vector as Vector
import os

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Create Lxs file in rep
#==============================================================================
def createLxs(rep):
    
    # File open
    name = VARS[1].get()
    file = open(rep+'/'+name+'.lxs', 'w')
    file.write('# Lux Render main file - written by Cassiopee -\n')

    # Sampler
    file.write('Renderer "sampler"\n')
    file.write('\n')
    file.write('Sampler "metropolis"\n')
    file.write('    "float largemutationprob" [0.4]\n')
    file.write('\n')
    file.write('Accelerator "qbvh"\n')
    file.write('\n')

    # Integrators
    file.write('SurfaceIntegrator "bidirectional"\n')
    file.write('      "integer eyedepth" [48]\n')
    file.write('      "integer lightdepth" [48]\n')
    file.write('      "string lightpathstrategy" ["auto"]\n')
    file.write('      "string lightstrategy" ["auto"]\n')
    file.write('\n')
    file.write('VolumeIntegrator "multi"\n')
    file.write('      "float stepsize" [1.0]\n')
    file.write('\n')

    # Pixel filter
    file.write('PixelFilter "mitchell"\n')
    file.write('      "bool supersample" ["true"]\n')
    file.write('      "float B" [0.3333333]\n')
    file.write('      "float C" [0.3333333]\n')
    file.write('\n')

    # Camera
    eye = CPlot.getState('posEye')
    cam = CPlot.getState('posCam')
    dir = CPlot.getState('dirCam')

    file.write('LookAt '+str(cam[0])+' '+str(cam[1])+' '+str(cam[2])+' '+
               str(eye[0])+' '+str(eye[1])+' '+str(eye[2])+' '+
               str(dir[0])+' '+str(dir[1])+' '+str(dir[2])+'\n')
    file.write('Camera "perspective" "float fov" [50]\n')
    file.write('\n')

    # Film
    resolution = VARS[0].get()
    size = resolution.split('x'); xres = size[0]; yres = size[1]
    file.write('Film "fleximage"\n')
    file.write('       "integer xresolution" ['+xres+']\n')
    file.write('       "integer yresolution" ['+yres+']\n')
    file.write('\n')

    # World
    file.write('WorldBegin\n')
    file.write('Include "Scene/Materials.lxm"\n')
    file.write('Include "Scene/Geometry.lxo"\n')
    file.write('Include "Scene/Volumes.lxv"\n')
    file.write('\n')

    # Lights
    writeCassiopeeLamps(file)

    file.write('\n')
    file.write('Exterior "world"\n')
    file.write('WorldEnd\n')
    file.close()
    return

#==============================================================================
# Create geometry file
# Contient pour chaque zone:
# La reference a son fichier stl
# Son material (surface, interieur, exterieur)
#==============================================================================
def createGeo(rep):
    file = open(rep+'/Geometry.lxo', 'w')
    file.write('# Lux Render geometries - written by Cassiopee -\n')
    
    zones = Internal.getZones(CTK.t)
    c = 0
    for z in zones:
        material = 'Solid'; color = 'White'
        ri = Internal.getNodesFromName1(z, '.RenderInfo')
        if ri != []:
            # Material
            mt = Internal.getNodesFromName1(ri[0], 'Material')
            if mt != []: material = Internal.getValue(mt[0])
            # Shader parameters
            sp = Internal.getNodesFromName1(ri[0], 'ShaderParameters')
            if sp != []: intensity = 50.*sp[0][1][0]
            else: intensity = 50.
            
        if material == 'Light':
            file.write('AttributeBegin # light_'+str(c)+'\n')
            file.write('NamedMaterial "material'+str(c)+'"\n')
            file.write('LightGroup "'+str(c)+'"\n')
            file.write('AreaLightSource "area"\n')
            file.write('     "float importance" [1.0]\n')
            file.write('     "float power" ['+str(intensity)+']\n')
            file.write('      "float efficacy" [17.0]\n')
            file.write('      "texture L" ["lamp blackbody"]\n')
            file.write('      "integer nsamples" [1]\n')
            file.write('      "float gain" [0.1]\n')
            file.write('      Shape "stlmesh"\n')
            file.write('      "integer nsubdivlevels" 1\n')
            file.write('      "string filename" ["Scene/mesh_'+str(c)+'.stl"]\n')
            file.write('AttributeEnd # ""\n')
            file.write('\n')
        else:
            file.write('AttributeBegin # mesh_'+str(c)+'\n')
            file.write('NamedMaterial "material'+str(c)+'"\n')
            file.write('Exterior "world"\n')
            #file.write('Interior "air"\n')
            file.write('Shape "stlmesh"\n')
            # not already implemented in lux?
            #file.write('     "bool smooth" ["true"]\n')
            file.write('     "integer nsubdivlevels" 1\n')
            file.write('     "string filename" ["Scene/mesh_'+str(c)+'.stl"]\n')
            file.write('AttributeEnd # ""\n')
            file.write('\n')
        c += 1
    file.close()
    return

#==============================================================================
# write matte dans file
# Update dict qui retient si les textures sont deja ecrites
# Cree le matrial c
#==============================================================================
def writeMatte0(file, dict, c, colorR, colorG, colorB, scale):
    name = 'matte'
    if name not in dict: dict[name] = 0
    
    file.write('MakeNamedMaterial "material'+str(c)+'"\n')
    file.write('      "color Kd" ['+str(colorR)+' '+str(colorG)+' '+str(colorB)+']\n')
    file.write('      "float sigma" [8.000000000000000]\n')
    file.write('      "string type" ["matte"]\n')
    file.write('\n')

#==============================================================================
# write glass dans file
# Update dict qui retient si les textures sont deja ecrites
# Cree le matrial c
#==============================================================================
def writeGlass0(file, dict, c, colorR, colorG, colorB, scale):
    name = 'glass'
    if name not in dict: dict[name] = 0
    
    file.write('MakeNamedMaterial "material'+str(c)+'"\n')
    file.write('      "bool architectural" ["false"]\n')
    file.write('      "float cauchyb" [0.000000000000000]\n')
    file.write('      "float film" [0.000000000000000]\n')
    file.write('      "float filmindex" [1.333299994468689]\n')
    file.write('      "float index" [1.519000053405762]\n')
    #file.write('      "color Kr" [0.69999999 0.69999999 0.69999999]\n')
    file.write('      "color Kr" ['+str(colorR)+' '+str(colorG)+' '+str(colorB)+']\n')
    #file.write('      "color Kr" [0.69999999 0.69999999 0.69999999]\n')
    file.write('      "color Kt" [1.00000000 1.00000000 1.00000000]\n')
    file.write('      "string type" ["glass"]\n')
    file.write('\n')

#==============================================================================
# write chrome dans file
# Update dict qui retient si les textures sont deja ecrites
# Cree le material c
#==============================================================================
def writeChrome0(file, dict, c, colorR, colorG, colorB, scale):
    name = 'chrome'
    if name not in dict: dict[name] = 0
    
    file.write('MakeNamedMaterial "material'+str(c)+'"\n')   
    file.write('      "float film" [0.000000000000000]\n')
    file.write('      "float filmindex" [1.333299994468689]\n')
    #file.write('      "color Kr" [0.69999999 0.69999999 0.69999999]\n')
    file.write('      "color Kr" ['+str(colorR)+' '+str(colorG)+' '+str(colorB)+']\n')
    #file.write('      "color Ks" [0.04000000 0.04000000 0.04000000]\n')
    file.write('      "color Ks" ['+str(0.1*colorR)+' '+str(0.1*colorG)+' '+str(0.1*colorB)+']\n')
    file.write('      "color Ks" [0.04000000 0.04000000 0.04000000]\n')
    file.write('      "float uroughness" [0.100000001490116]\n')
    file.write('      "float vroughness" [0.100000001490116]\n')
    file.write('      "string type" ["shinymetal"]\n')
    file.write('\n')

#==============================================================================
# write chrome dans file
# Update dict qui retient si les textures sont deja ecrites
# Cree le material c
#==============================================================================
def writeMetal0(file, dict, c, colorR, colorG, colorB, scale):
    name = 'metal'
    if name not in dict: dict[name] = 0
    
    file.write('MakeNamedMaterial "material'+str(c)+'"\n')   
    file.write('	"bool multibounce" ["false"]\n')
    file.write('        "color Kd" ['+str(colorR)+' '+str(colorG)+' '+str(colorB)+']\n')
    #file.write('	"color Kd" [0.11725003 0.11927158 0.10916383]\n')
    file.write('	"color Ks" [0.23300919 0.23300919 0.23300919]\n')
    file.write('	"float index" [0.000000000000000]\n')
    file.write('	"float uroughness" [0.249940723180771]\n')
    file.write('	"float vroughness" [0.249940723180771]\n')
    file.write('	"float sigma" [0.000000000000000]\n')
    file.write('	"string type" ["glossy"]\n')
    file.write('\n')

#==============================================================================
# write marble dans file
# Update dict qui retient si les textures sont deja ecrites
# Cree le material c
#==============================================================================
def writeMarble0(file, dict, c, colorR, colorG, colorB, scale):
    name = 'marble'
    if name not in dict: dict[name] = 0
    
    file.write('Texture "Texture" "float" "blender_marble"\n')
    file.write('      "float bright" [1.000000000000000]\n')
    file.write('      "float contrast" [3.435120105743408]\n')
    file.write('      "string type" ["sharp"]\n')
    file.write('      "string noisetype" ["hard_noise"]\n')
    file.write('      "string noisebasis" ["voronoi_f2"]\n')
    file.write('      "string noisebasis2" ["tri"]\n')
    file.write('      "float noisesize" [1.572026491165161]\n')
    file.write('      "float turbulence" [17.708000183105469]\n')
    file.write('      "integer noisedepth" [1]\n')
    file.write('      "string coordinates" ["global"]\n')
    file.write('      "vector translate" [0.0 0.0 0.0]\n')
    file.write('      "vector rotate" [0.0 0.0 0.0]\n')
    file.write('      "vector scale" ['+str(scale)+' '+str(scale)+' '+str(scale)+']\n')
    file.write('\n')

    file.write('Texture "Texture.002" "color" "mix"\n')
    file.write('	"texture amount" ["Texture"]\n')
    file.write('	"color tex1" [0.74838847 0.74838847 0.74838847]\n')
    #file.write('    "color tex1" ['+str(colorR)+' '+str(colorG)+' '+str(colorB)+']\n')
    #file.write('	"color tex2" [0.52770847 0.52770847 0.52770847]\n')
    file.write('    "color tex2" ['+str(colorR)+' '+str(colorG)+' '+str(colorB)+']\n')
    file.write('\n')
        
    file.write('MakeNamedMaterial "material'+str(c)+'"\n')   
    file.write('      "bool multibounce" ["false"]\n')
    file.write('      "texture Kd" ["Texture.002"]\n')
    file.write('      "color Ks" [0.07657672 0.07657672 0.07657672]\n')
    file.write('      "float index" [0.0]\n')
    file.write('      "float uroughness" [0.015000002458692]\n')
    file.write('      "float vroughness" [0.015000002458692]\n')
    file.write('      "float sigma" [0.0]\n')
    file.write('      "string type" ["glossy"]\n')
    file.write('\n')
    
#==============================================================================
# write wood (alder) dans file
# Update dict qui retient si les textures sont deja ecrites
# Cree le matrial c
#==============================================================================
def writeWood0(file, dict, c, colorR, colorG, colorB, scale):
    name = 'wood alder'
    if name not in dict: dict[name] = 0
    file.write('Texture "wood alder part 1" "float" "blender_wood"\n')
    file.write('     "float bright" [1.000000000000000]\n')
    file.write('     "float contrast" [2.000000000000000]\n')
    file.write('     "string noisebasis" ["blender_original"]\n')
    file.write('     "string noisebasis2" ["sin"]\n')
    file.write('     "float noisesize" [0.250000000000000]\n')
    file.write('     "string noisetype" ["hard_noise"]\n')
    file.write('     "float turbulence" [2.000000000000000]\n')
    file.write('     "string type" ["ringnoise"]\n')
    file.write('     "string coordinates" ["local"]\n')
    file.write('     "vector translate" [0.0 0.0 0.0]\n')
    file.write('     "vector rotate" [0.0 0.0 0.0]\n')
    file.write('     "vector scale" ['+str(scale)+' '+str(scale)+' '+str(scale)+']\n')
    file.write('\n')

    file.write('Texture "wood alder part 2" "float" "blender_wood"\n')
    file.write('     "float bright" [1.000000000000000]\n')
    file.write('     "float contrast" [2.000000000000000]\n')
    file.write('     "string noisebasis" ["blender_original"]\n')
    file.write('     "string noisebasis2" ["saw"]\n')
    file.write('     "float noisesize" [0.250000000000000]\n')
    file.write('     "string noisetype" ["hard_noise"]\n')
    file.write('     "float turbulence" [2.000000000000000]\n')
    file.write('     "string type" ["ringnoise"]\n')
    file.write('     "string coordinates" ["local"]\n')
    file.write('     "vector translate" [0.0 0.0 0.0]\n')
    file.write('     "vector rotate" [0.0 0.0 0.0]\n')
    file.write('     "vector scale" ['+str(scale)+' '+str(scale)+' '+str(scale)+']\n')	
    file.write('\n')

    file.write('Texture "wood alder mix" "float" "mix"\n')
    file.write('     "float amount" [0.500000000000000]\n')
    file.write('     "texture tex1" ["wood alder part 1"]\n')
    file.write('     "texture tex2" ["wood alder part 2"]\n')
    file.write('\n')

    file.write('Texture "1b67fd720273dfb27261d" "float" "scale"\n')
    file.write('     "float tex1" [0.002000000094995]\n')
    file.write('     "texture tex2" ["wood alder mix"]\n')
    file.write('\n')
        
    file.write('Texture "wood alder part 3" "float" "blender_clouds"\n')
    file.write('     "float bright" [0.009999999776483]\n')
    file.write('     "float contrast" [1.200000047683716]\n')
    file.write('     "string noisetype" ["hard_noise"]\n')
    file.write('     "string noisebasis" ["voronoi_f3"]\n')
    file.write('     "float noisesize" [0.250000000000000]\n')
    file.write('     "integer noisedepth" [2]\n')
    file.write('     "string coordinates" ["local"]\n')
    file.write('     "vector translate" [0.0 0.0 0.0]\n')
    file.write('     "vector rotate" [0.0 0.0 0.0]\n')
    file.write('     "vector scale" ['+str(scale)+' '+str(scale)+' '+str(scale)+']\n')
    file.write('\n')
        
    file.write('Texture "665d14f3da96af422fdff" "float" "scale"\n')
    file.write('     "float tex1" [1.000000000000000]\n')
    file.write('     "texture tex2" ["wood alder part 3"]\n')

    file.write('Texture "wood alder diffuse 1" "color" "mix"\n')
    file.write('     "texture amount" ["wood alder mix"]\n')
    file.write('     "color tex1" [0.64313728 0.45594707 0.29407442]\n')
    file.write('     "color tex2" [0.36960801 0.25732201 0.16374999]\n')
    file.write('\n')
        
    file.write('Texture "wood alder diffuse 2" "color" "mix"\n')
    file.write('     "texture amount" ["665d14f3da96af422fdff"]\n')
    file.write('     "texture tex1" ["wood alder diffuse 1"]\n')
    file.write('     "color tex2" [0.17332031 0.09053276 0.04070793]\n')
    file.write('\n')

    file.write('MakeNamedMaterial "material'+str(c)+'"\n')
    file.write('      "texture bumpmap" ["1b67fd720273dfb27261d"]\n')
    file.write('      "bool multibounce" ["false"]\n')
    file.write('      "texture Kd" ["wood alder diffuse 2"]\n')
    file.write('      "color Ks" [0.04000000 0.04000000 0.04000000]\n')
    file.write('      "float index" [0.000000000000000]\n')
    file.write('      "float uroughness" [0.100000001490116]\n')
    file.write('      "float vroughness" [0.100000001490116]\n')
    file.write('      "float sigma" [0.000000000000000]\n')
    file.write('      "string type" ["glossy"]\n')
    file.write('\n')

#==============================================================================
def writeCassiopeeLamps(file):
    type = VARS[2].get()
    
    eye = CPlot.getState('posEye')
    cam = CPlot.getState('posCam')
    dir = CPlot.getState('dirCam')
    d = Vector.sub(eye, cam)
    n = Vector.cross(d, dir)
    dir = Vector.normalize(dir)
    n = Vector.normalize(n)
    norm = Vector.norm(d)
    d = Vector.normalize(d)
    n = Vector.sub(n, dir)
    n = Vector.add(n, d)
    n = Vector.mul(0.4*norm, n)
    pos = Vector.sub(eye, n)

    if type == 'Interior': pass        
##         # distant light
##         file.write('AttributeBegin\n')
##         file.write('LightGroup "default"\n')
##         file.write('Exterior "world"\n')
##         file.write('LightSource "distant"\n')
##         file.write('      "point from" ['+
##                    str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+
##                    '] "point to" ['+
##                    str(eye[0])+' '+str(eye[1])+' '+str(eye[2])+']\n')
##         file.write('      "color L" [5 5 5]\n')
##         file.write('AttributeEnd\n')

    else: # Exterior
        # sun/sky
        sundir = Vector.sub(pos, eye)
        sundir = Vector.normalize(sundir)
        file.write('AttributeBegin\n')
        file.write('LightGroup "default"\n')
        file.write('Exterior "world"\n')
        file.write('LightSource "sunsky"\n')
        file.write('      "vector sundir" ['+
                   str(sundir[0])+' '+str(sundir[1])+' '+str(sundir[2])+
                   ']\n')
        file.write('      "float turbidity" [2.0]\n')
        file.write('      "float gain" [0.005]\n')
        file.write('AttributeEnd\n')

    
#==============================================================================
# Create material file
#==============================================================================
def createMat(rep):
    # already created materials
    dict = {}
    
    file = open(rep+'/Materials.lxm', 'w')
    file.write('# Lux Render materials - written by Cassiopee -\n')

    # world (clear medium)
    file.write('MakeNamedVolume "world" "clear"\n')
    file.write('      "float fresnel" [1.00]\n')
    file.write('      "color absorption" [0.0 0.0 0.0]\n')
    file.write('\n')

    # air (avec des particules dedans)
    file.write('MakeNamedVolume "air" "homogeneous"\n')
    file.write('       "float fresnel" [1.0]\n')
    file.write('       "color g" [0.0 0.0 0.0]\n')
    file.write('       "color sigma_a" [0.0 0.0 0.0]\n')
    file.write('       "color sigma_s" [0.036 0.036 0.036]\n')
    file.write('\n')

    # lamp
    file.write('Texture "lamp blackbody" "color" "blackbody"\n')
    file.write('      "float temperature" [6500.000000000000000]\n')
    file.write('\n')
    file.write('MakeNamedMaterial "lamp"\n')
    file.write('     "string type" ["null"]\n')
    file.write('\n')

    # cree un materiau par zone
    zones = Internal.getZones(CTK.t)
    c = 0
    for z in zones:
        material = 'Solid'; color = 'White'; mode = 0; blending = 1 # default
        ri = Internal.getNodesFromName1(z, '.RenderInfo')
        if ri != []:
            # Material
            mt = Internal.getNodesFromName1(ri[0], 'Material')
            if mt != []: material = Internal.getValue(mt[0])
            # Color
            co = Internal.getNodesFromName1(ri[0], 'Color')
            if co != []: color = Internal.getValue(co[0])
            # Blending
            co = Internal.getNodesFromName1(ri[0], 'Blending')
            if co != []: blending = Internal.getValue(co[0])

        s = color.split(':')
        if len(s) == 2 and s[0] == 'Iso':
            vref = C.getVarNames(z)[0]
            for pos in xrange(len(vref)):
                if (vref[pos] == s[1]): break
            
            if pos == len(vref): color = 'White'; mode = 0
            else: color = 'Iso'; mode = pos+1; material = 'Iso'
        # traduction color
        if color[0] == '#':
            colorR = color[1:3]; colorG = color[3:5]; colorB = color[5:]
            colorR = int(colorR, 16); colorR = colorR / 255.
            colorG = int(colorG, 16); colorG = colorG / 255.
            colorB = int(colorB, 16); colorB = colorB / 255.
        elif color == 'White': colorR = 1; colorG = 1; colorB = 1
        elif color == 'Black': colorR = 0; colorG = 0; colorB = 0
        elif color == 'Grey': colorR = 0.69; colorG = 0.69; colorB = 0.69
        elif color == 'Blue': colorR = 0; colorG = 0; colorB = 1
        elif color == 'Red': colorR = 1; colorG = 0; colorB = 0
        elif color == 'Green': colorR = 0; colorG = 1; colorB = 0
        elif color == 'Yellow': colorR = 1; colorG = 1; colorB = 0
        elif color == 'Orange': colorR = 0.94; colorG = 0.737; colorB = 0.06
        elif color == 'Magenta': colorR = 1; colorG = 0; colorB = 1
        elif color == 'Brown': colorR = 0.588; colorG = 0.294; colorB = 0
        else: coloR = 1; colorG = 1; colorB = 1

        # Scale (utlise pour scaler les textures)
        bb = G.bbox(z)
        rx = bb[3]-bb[0]; ry = bb[4]-bb[1]; rz = bb[5]-bb[2]
        scale = 0.5 * min(rx, ry, rz)
            
        if material == 'Solid':
            writeMatte0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Glass':
            writeGlass0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Chrome':
            writeChrome0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Metal':
            writeMetal0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'XRay':
            writeMatte0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Wood':
            writeWood0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Marble':
            writeMarble0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Granite':
            writeMatte0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Smoke':
            writeMatte0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Brick':
            writeMatte0(file, dict, c, colorR, colorG, colorB, scale)
        elif material == 'Light':
            writeGlass0(file, dict, c, colorR, colorG, colorB, scale)
        else:
            writeMatte0(file, dict, c, colorR, colorG, colorB, scale)
        c += 1
    file.close()
    return

#==============================================================================
# Create volume file
#==============================================================================
def createVol(rep):
    file = open(rep+'/Volumes.lxv', 'w')
    file.write('# Lux Render volumes - written by Cassiopee -\n')
    file.close()
    return

#==============================================================================
# Create stl files in rep
#==============================================================================
def createStl(rep):
    zones = Internal.getZones(CTK.t)
    c = 0
    for z in zones:
        z = C.convertArray2Tetra(z)
        
        # Get Material/Color
        material = 'Solid'; color = 'White'; mode = 0  # default
        ri = Internal.getNodesFromName1(z, '.RenderInfo')
        if ri != []:
            mt = Internal.getNodesFromName1(ri[0], 'Material')
            if mt != []: material = Internal.getValue(mt[0])
            co = Internal.getNodesFromName1(ri[0], 'Color')
            if co != []: color = Internal.getValue(co[0])
        s = color.split(':')
        if len(s) == 2 and s[0] == 'Iso':
            vref = C.getVarNames(z)[0]
            for pos in xrange(len(vref)):
                if vref[pos] == s[1]: break
            
            if pos == len(vref): color = 'White'; mode = 0
            else: color = 'Iso'; mode = pos+1; material = 'Iso'

        # write file
        C.convertPyTree2File(z, rep+'/mesh_'+str(c)+'.stl')
        c += 1
    
    return
    
#==============================================================================
# Create all Lux Render files
#==============================================================================
def createFiles():
    if CTK.t == []: return
    rep = VARS[1].get()
    dir = os.path.dirname(CTK.FILE)
    rep = os.path.join(dir, rep)
    a = os.access(rep, os.F_OK)
    if a == False: os.mkdir(rep)
    sceneRep = rep+'/'+'Scene'
    a = os.access(sceneRep, os.F_OK)
    if a == False: os.mkdir(sceneRep)
    createStl(sceneRep)
    createGeo(sceneRep)
    createMat(sceneRep)
    createVol(sceneRep)
    createLxs(rep)
    CTK.TXT.insert('START', 'Luxrender files written.\n')
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkLuxRender', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Export to LuxRender.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkLuxRender')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Image size
    V = TK.StringVar(win); V.set('800x600'); VARS.append(V)
    if 'tkLuxRenderSize' in CTK.PREFS: 
        V.set(CTK.PREFS['tkLuxRenderSize'])
    # -1- Rep name
    V = TK.StringVar(win); V.set('LuxRender'); VARS.append(V)
    if 'tkLuxRenderOutput' in CTK.PREFS: 
        V.set(CTK.PREFS['tkLuxRenderOutput'])
    # -2- Interior / exterior
    V = TK.StringVar(win); V.set('Exterior'); VARS.append(V)   
    if 'tkLuxRenderType' in CTK.PREFS: 
        V.set(CTK.PREFS['tkLuxRenderType'])

    # - Type of scene
    B = TTK.OptionMenu(Frame, VARS[2], 'Exterior', 'Interior')
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of scene.')

    # - Rep -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White')
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Directory name for output.')

    # - Image size -
    B = TTK.OptionMenu(Frame, VARS[0], '320x200', '800x600', '1024x768',
                       '1600x1200')
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Output image resolution.')

    # - Export scene -
    B = TTK.Button(Frame, text="Export", command=createFiles)
    B.grid(row=2, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create all lux render files.')
    
#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW)

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkLuxRenderSize'] = VARS[0].get()
    CTK.PREFS['tkLuxRenderOutput'] = VARS[1].get()
    CTK.PREFS['tkLuxRenderType'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('800x600')
    VARS[1].set('LuxRender')
    VARS[2].set('Exterior')
    CTK.PREFS['tkLuxRenderSize'] = VARS[0].get()
    CTK.PREFS['tkLuxRenderOutput'] = VARS[1].get()
    CTK.PREFS['tkLuxRenderType'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)
    
#==============================================================================
if (__name__ == "__main__"):
    import sys
    if (len(sys.argv) == 2):
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkLuxRender '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
