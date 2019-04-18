# - tkPovRay -
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
import subprocess

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Save PovRay file
#==============================================================================
def savePovFile():
    if CTK.t == []: return

    # Sauvegarde toutes les zones au format pov
    rep = VARS[2].get()
    dir = os.path.dirname(CTK.FILE)
    rep = os.path.join(dir, rep)
    os.chdir(rep)
    zones = Internal.getZones(CTK.t)
    c = 0; files = []
    colors = []; materials = []; blendings = []; shader1s = []
    scales = []; centers = []
    for z in zones:
        z = C.convertArray2Tetra(z)

        # Scale (utlise pour scaler les textures)
        bb = G.bbox(z)
        rx = bb[3]-bb[0]; ry = bb[4]-bb[1]; rz = bb[5]-bb[2]
        scale = min(rx, ry, rz)
        scales.append(scale)
        centers.append([bb[0]+0.5*rx, bb[1]+0.5*ry, bb[2]+0.5*rz])
        
        # Material/Color/Blending
        material = 'Solid'; color = 'White'; mode = 0;
        blending = 1; shader1 = 1.
        ri = Internal.getNodesFromName1(z, '.RenderInfo')
        if (ri != []):
            # Material
            mt = Internal.getNodesFromName1(ri[0], 'Material')
            if (mt != []): material = Internal.getValue(mt[0])
            # Color
            co = Internal.getNodesFromName1(ri[0], 'Color')
            if (co != []): color = Internal.getValue(co[0])
            # Blending
            co = Internal.getNodesFromName1(ri[0], 'Blending')
            if (co != []): blending = Internal.getValue(co[0])
            # Shader parameter 1
            co = Internal.getNodesFromName1(ri[0], 'ShaderParameters')
            if (co != []): shader1 = co[0][1][0]
            else: shader1 = 1.
        s = color.split(':')
        if (len(s) == 2 and s[0] == 'Iso'): # couleur = iso field
            vref = C.getVarNames(z)[0]
            for pos in xrange(len(vref)):
                if (vref[pos] == s[1]): break
            
            if (pos == len(vref)): color = 'White'; mode = 0
            else: color = 'Iso'; mode = pos+1
        # traduction color si #FFFFFF
        if (color[0] == '#'):
            colorR = color[1:3]; colorG = color[3:5]; colorB = color[5:]
            colorR = int(colorR, 16); colorR = colorR / 255.
            colorG = int(colorG, 16); colorG = colorG / 255.
            colorB = int(colorB, 16); colorB = colorB / 255.
            color = 'rgbf<'+str(colorR)+','+str(colorG)+','+str(colorB)+'>'
        colors.append(color); materials.append(material)
        blendings.append(blending); shader1s.append(shader1)
        
        nt = C.newPyTree(['Base'])
        nt[2][1][2].append(z)
        try:
            if (mode == 0):
                C.convertPyTree2File(nt, 'mesh_'+str(c)+'.pov')
            else:
                C.convertPyTree2File(nt, 'mesh_'+str(c)+'.pov',
                                     colormap=mode) # avec iso
            files.append('mesh_'+str(c)+'.pov')
            c += 1
        except: pass

    # Cam position
    eye = CPlot.getState('posEye')
    cam = CPlot.getState('posCam')
    dir = CPlot.getState('dirCam')

    # Ecriture du fichier PovRay
    file = open('scene.pov', 'w')
    file.write('// POV-Ray version 3.6 scenery file written by *Cassiopee*\n')
    file.write('// Please render this file with :\n')
    file.write('// povray -W800 -H600 +a0.3 +SP16 scene.pov +P\n')
    file.write('#version 3.6;\n')
    file.write('#include "colors.inc"\n')
    file.write('#include "textures.inc"\n')
    file.write('#include "woods.inc"\n')
    file.write('#include "stones.inc"\n')
    
    # Brushed metal texture
    file.write('#declare Brushed_Depth = 10; // Bump size\n')
    file.write('#declare Brushed_Pigment = pigment {colour rgb 0.73} \n')
    file.write('#declare Brushed_Finish = finish {ambient 0 diffuse 0.95 specular 0.96 roughness 0.0005 phong 0.43 phong_size 25 brilliance 3.15 reflection 0.33 metallic metallic on }\n')
    file.write('// The brushed metal texture.\n')
    file.write('#declare Brushed_Texture = texture { average texture_map { [ pigment {Brushed_Pigment} normal {wood +Brushed_Depth ramp_wave rotate 90*x scale 50} finish {Brushed_Finish} ] [pigment {Brushed_Pigment} normal {wood -Brushed_Depth ramp_wave rotate 90*x scale 50} finish {Brushed_Finish} ] } }\n')

    # XRay texture
    file.write('#declare XRayTexture2 = texture { pigment { slope{'+'<'
               +str(cam[0])+' , '+str(cam[1])+' , '+str(cam[2])+'> - <'
               +str(eye[0])+' , '+str(eye[1])+' , '+str(eye[2])+'>}\n')
    file.write('pigment_map {[0 color rgbt 2*<1,1,1,0.1>] [0.75 color rgbt <0.1,0.6,2,1>*1] [1    color rgbt <1,1,1,1>] } } finish {ambient 3} }\n')

    # - Radiosity -
    #file.write('global_settings { assumed_gamma 1 radiosity { \n')
    #file.write('pretrace_start 0.08 \n')
    #file.write('pretrace_end   0.02 \n')
    #file.write('count 50 \n')
    #file.write('error_bound 0.5 \n')
    #file.write('recursion_limit 1 } } \n')
    
    # - Camera -
    file.write('#declare Cam0 = camera {angle 50 \n')
    file.write('#location  <'+str(cam[0])+' , '+str(cam[1])+' , '+
               str(cam[2])+'> \n')
    file.write('#look_at <'+str(eye[0])+' , '+str(eye[1])+' , '+
               str(eye[2])+'>\n')
    file.write('#direction <-1,0,0>\n')
    file.write('#sky <'+str(dir[0])+','+str(dir[1])+','+str(dir[2])+'> }\n')
    # focal blur (experimental)
    #file.write('#focal_point <'+str(eye[0])+' , '+str(eye[1])+' , '+
    #           str(eye[2])+'>\n')
    #file.write('#aperture 0.4\n')
    #file.write('#blur_samples 20 }\n')
    file.write('camera{Cam0}\n')

    # - Lumieres: point -
    #d = Vector.sub(eye, cam)
    #n = Vector.cross(d, dir)
    dir = Vector.normalize(dir)
    #n = Vector.normalize(n)
    #norm = Vector.norm(d)
    #d = Vector.normalize(d)
    #n = Vector.sub(n, dir)
    #n = Vector.add(n, d)
    #n = Vector.mul(0.4*norm, n)
    #pos = Vector.sub(eye, n)

    dir = Vector.mul(0.1, dir)
    pos = Vector.add(cam, dir)
    
    c = 0; light = 0
    for f in files:
        material = materials[c]
        if (material == 'Light'):
            xc = centers[c]; color = colors[c]
            light = 1
            intensity = shader1s[c]*5.
            file.write('light_source{<'+str(xc[0])+' , '+str(xc[1])+' , '+
                       str(xc[2])+'> color '+color+'*'+str(intensity)+'}\n')
        c += 1

    if (light == 0): # pas de lumiere dans l'arbre, on met celle par defaut
        file.write('light_source{<'+str(pos[0])+' , '+str(pos[1])+' , '+
                   str(pos[2])+'> color White*4}\n')
        
    # - Background -
    bckgrd = VARS[0].get()
    if (bckgrd == 'Blue sky'):
        # Ciel bleu
        file.write('sky_sphere { pigment { gradient <0,0,1> turbulence 0\n')
        file.write('       color_map { [0.00 rgb <0.6,0.7,1.0>]\n')
        file.write('                   [0.35 rgb <0.1,0.2,0.8>]\n')
        file.write('                   [0.65 rgb <0.1,0.2,0.8>]\n')
        file.write('                   [1.00 rgb <0.6,0.7,1.0>]\n')
        file.write('                 }\n')
        file.write('       scale 2\n')     
        file.write('     } // end of pigment\n')
        file.write('  } //end of skysphere\n')
    elif (bckgrd == 'Cloudy sky'): # on pourrait faire beaucoup mieux
        file.write('sky_sphere { \n')
        file.write('pigment{ bozo turbulence 0.76\n')
        file.write('       color_map { [0.5 rgb <0.20, 0.20, 1.0>]\n')
        file.write('                   [0.6 rgb <1,1,1>]\n')
        file.write('                   [1.0 rgb <0.5,0.5,0.5>]\n')
        file.write('                 }\n')
        file.write('     } // end of pigment\n')
        file.write('  } //end of skysphere\n')

        # Avec les macros de pov, pas satisfaisant
        #file.write('#include "skies.inc" \n')
        #file.write("object{ O_Cloud1 rotate 90*x}\n")
        #file.write("sphere { <0,0,0>, 100 \n")
        #file.write("texture {T_Cloud3} scale 100 }\n")
        
    elif (bckgrd == 'Starfield'):
        file.write('#include "stars.inc"\n')
        file.write('sphere { <0,0,0>, 1\n')
        file.write('texture { Starfield1 }\n')
        file.write('scale 10000\n')
        file.write('  } //end of sphere\n')
        
    elif (bckgrd == 'White'):
        file.write('sky_sphere { pigment { White } }\n')

    # Objets
    c = 0
    for f in files:
        color = colors[c]; material = materials[c]
        blend = str(1.-blendings[c])
        
        if (material == 'Solid' or material == 'None'): # OK
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            if (color != 'Iso'):
                file.write('texture{pigment{color '+color+' filter '+
                           blend+'}\n')
            else:
                file.write('texture{pigment{filter '+
                           blend+'}\n')
            file.write('finish {ambient 0.1 diffuse 0.1 reflection 0.05 phong 0.5 }}\n')
            file.write('}\n')
        elif (material == 'Flat'):
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            if (color != 'Iso'):
                file.write('texture{pigment{color '+color+' filter '+
                           blend+'}\n')
            else:
                file.write('texture{pigment{filter '+
                           blend+'}\n')
            file.write('finish {ambient 0.1 diffuse 0.1 reflection 0.0 phong 0.0 }}\n')
            file.write('}\n')
        elif (material == 'Glass'): # OK
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            if (color != 'Iso'):
                file.write('texture{pigment{color '+color+' filter '+
                           str(-0.1*blendings[c] +1.)+' }\n')
            else:
                 file.write('texture{pigment{filter '+
                            str(-0.1*blendings[c] +1.)+' }\n')
            file.write('finish {reflection 0.2 phong 0.7 }}\n')
            file.write('interior { ior 1.3 }\n')
            file.write('}\n')
        elif (material == 'Chrome'): # OK
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            if (color != 'Iso'):
                file.write('texture{pigment{color '+color+' filter '+
                           blend+'}\n')
            else:
                file.write('texture{pigment{filter '+
                           blend+'}\n')
            file.write('finish {ambient 0.25 brilliance 4 diffuse 0.5 reflection 0.4 specular 0.2 metallic roughness 1/80 }}\n')
            file.write('}\n')
        elif (material == 'Metal'):
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            if (color == 'White'):
                file.write('texture {pigment{color rgb 0.73}\n')
            elif (color != 'Iso'):
                file.write('texture {pigment{color '+color+' filter '+
                           blend+'}\n')
            else:
                file.write('texture {pigment{filter '+
                           blend+'}\n')
            file.write('finish {ambient 0 diffuse 0.95 specular 0.26 roughness 0.0005 phong 0.33 phong_size 2 brilliance 3.15 reflection 0.33 metallic metallic on } }\n')
            file.write('}\n')
        elif (material == 'XRay'):
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            file.write('texture { XRayTexture2 } }\n')
        elif (material == 'Wood'): # OK
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            if (color == 'White'):
                file.write('texture{T_Wood3 \n')
            elif (color == 'Black'):
                file.write('texture{T_Wood2 \n')
            elif (color == 'Blue'):
                file.write('texture{T_Wood31 \n')
            elif (color == 'Red'):
                file.write('texture{T_Wood6 \n')
            elif (color == 'Green'):
                file.write('texture{T_Wood32 \n')
            elif (color == 'Yellow'):
                file.write('texture{T_Wood35 \n')
            elif (color == 'Orange'):
                file.write('texture{T_Wood7 \n')
            elif (color == 'Magenta'):
                file.write('texture{T_Wood4 \n')
            else: file.write('texture{T_Wood32 \n')
            file.write("scale "+str(scales[c])+"\n");
            file.write('finish {ambient 0.7 brilliance 0.2 diffuse 0.1 reflection 0.01 specular 0.1 roughness 1/20 }}\n')
            file.write('}\n')
        elif (material == 'Marble' or material == 'Granite'):
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            if (color == 'White'):
                #file.write('texture{T_Grnt20 \n')
                file.write('texture{White_Marble \n')
            elif (color == 'Black'):
                file.write('texture{T_Grnt15 \n')
            elif (color == 'Blue'):
                file.write('texture{T_Wood6 \n')
            elif (color == 'Red'):
                file.write('texture{T_Wood28 \n')
            elif (color == 'Green'):
                file.write('texture{T_Grnt21 \n')
            elif (color == 'Yellow'):
                file.write('texture{T_Wood1 \n')
            elif (color == 'Orange'):
                file.write('texture{T_Wood13 \n')
            elif (color == 'Magenta'):
                file.write('texture{T_Grnt14 \n')
            else: file.write('texture{T_Grnt20 \n')
            file.write('scale '+str(scales[c])+'\n');
            file.write('finish {ambient 0.3 brilliance 0.3 diffuse 0.1 reflection 0.1 specular 0.1 }}\n')
            file.write('}\n')
        elif (material == 'Smoke'):
            file.write('#include "'+f+'"\n')
            file.write('object {mesh_'+str(c)+'\n')
            file.write('texture{pigment{color '+color+' filter 1.}\n')
            #file.write('texture{pigment { rgbt 1 }\n') 
            file.write('finish {ambient 0.1 diffuse 0.1  }}\n')
            file.write('hollow\n')
            file.write('interior{ //---------------------\n')
            file.write('media{ method 2 \n')
            file.write('emission 0. \n')
            file.write('scattering{ 1, // Type \n')
            file.write('<1,1,1>*0.2 // color of scattering haze \n')
            file.write('extinction  1. \n')
            file.write('// how fast the scattering media absorbs light \n')
            file.write('// useful if the media absorbs too much light \n')
            file.write('} // end scattering \n')
            file.write('density{ bozo \n')
            file.write('turbulence 8.0 \n')
            file.write('        color_map { \n')
            file.write('        [0.00 rgb 0] \n')
            file.write('        [0.05 rgb 0] \n')
            file.write('        [0.20 rgb 0.2] \n')
            file.write('        [0.30 rgb 0.6] \n')
            file.write('        [0.40 rgb 1] \n')
            file.write('        [1.00 rgb 1] \n')
            file.write('       } // end color_map \n')
            file.write('scale '+str(scales[c])+'\n');
            file.write('} // end of density  \n')
            file.write('samples 1,1   // 3,3 for adaptive sampling \n')
            file.write('intervals 10   // increase up to 15 \n')
            file.write('} // end of media --------------------------- \n')
            file.write('} // end of interior \n')

            file.write('}\n')
        c += 1

    file.close()
    os.chdir('..')
    
#==============================================================================
# Render scene using PovRay
#==============================================================================
def render():
    if (CTK.t == []): return
    rep = VARS[2].get()
    dir = os.path.dirname(CTK.FILE)
    rep = os.path.join(dir, rep)
    a = os.access(rep, os.F_OK)
    if (a == False): os.mkdir(rep)
    savePovFile()
    size = VARS[1].get()
    size = size.split('x')
    proc = subprocess.Popen('cd '+rep+'; povray -W'+size[0]+' -H'+size[1]+' +a0.3 +SP16 scene.pov +P', stdout=subprocess.PIPE, shell=True)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkPovRay', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Export to povRay ray tracer.\nCtrl+c to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkPovRay')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- background -
    V = TK.StringVar(win); V.set('Black'); VARS.append(V)
    if 'tkPovRayBackground' in CTK.PREFS: 
        V.set(CTK.PREFS['tkPovRayBackground'])
    # -1- Image size
    V = TK.StringVar(win); V.set('800x600'); VARS.append(V)
    if 'tkPovRaySize' in CTK.PREFS: 
        V.set(CTK.PREFS['tkPovRaySize'])
    # -2- Dir name (file.pov et file.png)
    V = TK.StringVar(win); V.set('PovRay'); VARS.append(V)    
    if 'tkPovRayOutput' in CTK.PREFS: 
        V.set(CTK.PREFS['tkPovRayOutput'])

    # - File -
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White')
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Directory name for output (.pov and .png).')

    # - Background selection -
    B = TTK.OptionMenu(Frame, VARS[0], 'Black', 'White', 'Blue sky', 
                       'Cloudy sky', 'Starfield' )
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Output background.')

    # - Image size -
    B = TTK.OptionMenu(Frame, VARS[1], '320x200', '800x600',
                       '1024x768', '1600x1200')
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Image resolution.')

    # - Render scene -
    B = TTK.Button(Frame, text="Render scene", command=render)
    B.grid(row=2, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Render pov scene.')
    
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
    CTK.PREFS['tkPovRayBackground'] = VARS[0].get()
    CTK.PREFS['tkPovRaySize'] = VARS[1].get()
    CTK.PREFS['tkPovRayOutput'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('Black')
    VARS[1].set('800x600')
    VARS[2].set('PovRay')
    CTK.PREFS['tkPovRayBackground'] = VARS[0].get()
    CTK.PREFS['tkPovRaySize'] = VARS[1].get()
    CTK.PREFS['tkPovRayOutput'] = VARS[2].get()
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
    (win, menu, file, tools) = CTK.minimal('tkPovRay '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
