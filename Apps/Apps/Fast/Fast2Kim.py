# Converter Fast probe file to KIM files
# orig: par Gabriel Reboul (ONERA)
import os
import struct
import numpy

import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Converter.Distributed as Distributed
import Converter.PyTree as C
import Post.PyTree as P

#=========================================
# Data
#=========================================
INT_SIZE = 4
DOUBLE_SIZE = 8

DEFAULT_SOUND_VELOCITY = 343.204
DEFAULT_DENSITY = 1.20432
DEFAULT_PRESSURE = 101325
DEFAULT_MACH = 0.
DEFAULT_BETAM = 0.
DEFAULT_NBMIC = 1
DEFAULT_FORM = 'solid'
DEFAULT_CORFLUX = '.False.'
DEFAULT_FORMFILE = 'fmt'
DEFAULT_TITLE = 'FAST4KIM'
DEFAULT_ID = 0
DEFAULT_OVERN = 0
DEFAULT_OMEGA = 0

#=======================================
class KIM:
                                   
    def __init__(self, config):
        '''Object constructor'''
        # Data input
        self.config = config
        self.kimDir = './KIM'
        self.kimEntreeDir = self.kimDir+os.sep+'ENTREES'         
        self.kimAeroDir = self.kimEntreeDir+os.sep+'AERO' 
        
        # Automatically executed functions
        self.createDirs()
        self.update_config_w_default_values()

    def update_config_w_default_values(self):
        '''Update the config dict file using default values'''
        if "density" not in self.config:
            self.config["density"] = DEFAULT_DENSITY
        if "sound_velocity" not in self.config:
            self.config["sound_velocity"] = DEFAULT_SOUND_VELOCITY
        if "pressure_infinity" not in self.config:
            self.config["pressure_infinity"]  = DEFAULT_PRESSURE
        if "mach_inf" not in self.config:
            self.config["mach_inf"] = DEFAULT_MACH
        if "betam" not in self.config:
            self.config["betam"] = DEFAULT_BETAM
        if "number_mic" not in self.config:
            self.config["number_mic"] = DEFAULT_NBMIC
        if "form"  not in self.config:
            self.config["form"] = DEFAULT_FORM 
        if "corflux" not in self.config:
            self.config["corflux"] = DEFAULT_CORFLUX
        if "title" not in self.config:
            self.config["title"] = DEFAULT_TITLE 
        if "dom_id" not in self.config:
            self.config["dom_id"] = DEFAULT_ID
        if "overn" not in self.config:
            self.config["overn"] = DEFAULT_OVERN            
        if "omega" not in self.config:
            self.config["omega"] = DEFAULT_OMEGA
            
    def createDirs(self):
        '''Create the KIM directories'''
        try: os.makedirs(self.kimDir)
        except: print ('KIM directory already exists.')
        try: os.makedirs(self.kimEntreeDir)
        except: print ('KIM ENTREES directory already exists.')
        try: os.makedirs(self.kimAeroDir)
        except: print ('KIM AERO directory already exists.')

    def write_fortran_block(self,mesh_file, vector, type_size, vector_type='d'):
        '''Write a vector in a binary file, with header and footer to comply with Fortran block format

        Parameters
        ----------
        mesh_file : file
        File object to write in
        x : nd.array
        Numpy array to be written in the the file
        type_size : int
        Size of the vector element
        '''
        tmp = type_size * len(vector)
        mesh_file.write(struct.pack('i', tmp))
        mesh_file.write(vector.astype(vector_type).tostring())
        mesh_file.write(struct.pack('i', tmp))

    def createKimDonn(self):
        '''Create the KIM donn.in'''
        # Creation of donn.in dictionnary
        donnD = {}
        donnD['titre']      = self.config['title']
        donnD['kirchhoff']  = '.False.'
        donnD['fwh_fixe']   = '.False.'
        if self.config['form'] == 'solid':
            donnD['fwh_soli']   = '.True.'
            donnD['fwh_mobi']   = '.False.'
        else:
            donnD['fwh_soli']   = '.False.'
            donnD['fwh_mobi']   = '.True.'          
        donnD['integrat']   = 0
        donnD['reflexion']  = '.False.'
        donnD['mail_def']   = '.True.'
        donnD['mail_tou']   = '.False.'
        donnD['cham_ins']   = '.True.'
        donnD['cham_tou']   = '.False.'
        donnD['cham_pro']   = '.True.'
        donnD['freefreq']   = '.False.'
        donnD['periodic']   = 0
        donnD['cor_flux']   = self.config['corflux']
        donnD['vit_moy']    = '.False.'
        donnD['ndom']       = self.config['number_dom']
        donnD['ntrimax']    = 16
        if self.config['form'] == 'porous':
            donnD['nvar']   = 5
        else:
            donnD['nvar']   = 1        
        donnD['ntemiss']    = self.config['number_dt']
        donnD['ntrecep']    = 2*self.config['number_dt']
        donnD['nper']       = 3
        donnD['fich_aero']  = '.True.'
        if self.config["form_file"] == 'bin':
            donnD['form_mail']  = '"unformatted"'
            donnD['form_cham']  = '"unformatted"'
        else:
            donnD['form_mail']  = '"formatted"'
            donnD['form_cham']  = '"formatted"'
        donnD['litendian']  = '.True.'    
        donnD['krepmaill']  = 1
        donnD['krepvites']  = 1
        donnD['ndupang']    = 1
        donnD['nduptem']    = 1
        donnD['lref']       = 1
        donnD['alpha']      = 0
        donnD['beta']       = 0
        donnD['alphar']     = 0
        donnD['betar']      = 0
        donnD['omega']      = self.config['omega']
        donnD['presinf']    = self.config['pressure_infinity']
        donnD['ainf']       = self.config['sound_velocity']
        donnD['machinf']    = self.config['mach_inf']
        donnD['duremis']    = self.config['number_dt']*self.config['dt']
        donnD['durecep']    = 2*self.config['number_dt']*self.config['dt']
        donnD['anacontri']  = '.False.'
        donnD['strucanac']  = 111
        donnD['dref']       = 1
        donnD['alpham']     = 0
        donnD['betam']      = self.config['betam']
        donnD['nmic_cart']  = self.config['number_mic']
        donnD['jmax_mic']   = -1
        donnD['trac_mic']   = '.True.'
        donnD['Mco']        = 1
        donnD['anime']      = '.True.'
        donnD['ktdebe']     = 1
        donnD['ktfine']     = min(self.config['number_dt'],100)
        donnD['ktpase']     = 2
        donnD['animd']      = '.False.'
        donnD['ktdebr']     = 0
        donnD['ktfinr']     = 0
        donnD['ktpasr']     = 0
        donnD['animr']      = '.False.'
        donnD['animt']      = '.False.'
        donnD['anitec']     = '.True.'
        donnD['fondam']     = 1
        donnD['inddeb']     = 1
        donnD['indfin']     = self.config['number_dt']
        donnD['nharli']     = 498
        donnD['nbharm']     = 498
        donnD['tiersdoc']   = '.False.'
        donnD['nivdb']      = '.True.'
        donnD['hardeb']     = 1
        donnD['harfin']     = 498
        donnD['nivdba']     = '.True.'
        donnD['mapdb']      = '.True.'
        donnD['mapdba']     = '.True.'    
        donnD['mapcomplex'] = '.False.'

        # Creation of donn.in lines
        donnLines = '''!
&CALCUL
g_titre      =  %(titre)s
g_kirchhoff  =  %(kirchhoff)s
g_fwh_fixe   =  %(fwh_fixe)s
g_fwh_mobi   =  %(fwh_mobi)s
g_fwh_soli   =  %(fwh_soli)s
g_integrat   =  %(integrat)d
g_reflexion  =  %(reflexion)s
g_mail_def   =  %(mail_def)s
g_mail_tou   =  %(mail_tou)s
g_cham_pro   =  %(cham_pro)s
g_cham_ins   =  %(cham_ins)s
g_cham_tou   =  %(cham_tou)s
g_freefreq   =  %(freefreq)s
g_periodic   =  %(periodic)d
g_cor_flux   =  %(cor_flux)s
g_vit_moy    =  %(vit_moy)s
/
!
&DIMENS
g_ndom       =  %(ndom)d 
g_ntrimax    =  %(ntrimax)d 
g_nvar       =  %(nvar)d
g_ntemiss    =  %(ntemiss)d 
g_ntrecep    =  %(ntrecep)d 
g_nper       =  %(nper)d
/
!
&PARAMS
g_fich_aero  =  .True.
g_form_mail  =  %(form_mail)s
g_form_cham  =  %(form_cham)s
g_litendian  =  %(litendian)s
g_krepmaill  =  %(krepmaill)d
g_krepvites  =  %(krepvites)d
g_ndupang    =  %(ndupang)d
g_nduptem    =  %(nduptem)d
/
!
&CONSTS
g_lref       =  %(lref)f
g_alpha      =  %(alpha)f
g_beta       =  %(beta)f
g_alphar     =  %(alphar)f
g_betar      =  %(betar)f
g_omega      =  %(omega)f
g_presinf    =  %(presinf)f
g_ainf       =  %(ainf)f
g_machinf    =  %(machinf)f
g_duremis    =  %(duremis)f
g_durecep    =  %(durecep)f
/
!
&PROCED
g_anacontri  =  %(anacontri)s
g_strucanac  =  %(strucanac)s
/
!
&MICROS
g_dref       =  %(dref)f
g_alpham     =  %(alpham)f
g_betam      =  %(betam)f
g_nmic_cart  =  %(nmic_cart)d
g_jmax_mic   =  %(jmax_mic)d
/
!
&CONTRO
g_trac_mic   =  %(trac_mic)s
g_mco        =  %(Mco)d
/
!
&ANIMAT
g_anime      =  %(anime)s
g_ktdebe     =  %(ktdebe)d
g_ktfine     =  %(ktfine)d
g_ktpase     =  %(ktpase)d
g_animd      =  %(animd)s
g_ktdebr     =  %(ktdebr)d
g_ktfinr     =  %(ktfinr)d
g_ktpasr     =  %(ktpasr)d 
g_animr      =  %(animr)s
g_animt      =  %(animt)s
g_anitec     =  %(anitec)s
/
!
&ANALYS
g_fondam     =  %(fondam)d
g_nbharm     =  %(nbharm)d
g_tiersdoc   =  %(tiersdoc)s
g_nivdb      =  %(nivdb)s
g_nharli     =  %(nharli)d
g_hardeb     =  %(hardeb)d
g_harfin     =  %(harfin)d
g_nivdba     =  %(nivdba)s
g_mapdb      =  %(mapdb)s
g_mapdba     =  %(mapdba)s
g_mapcomplex =  %(mapcomplex)s
/
!\n'''%(donnD)
        donn = open(self.kimEntreeDir+os.sep+'donn.in','w')
        donn.write(donnLines)
        donn.close()

    def createKimDimdomHead(self):
        '''Create the KIM dimdom.in head lines'''
        # Creation of dimdom.in dictionnary
        dimdomD = {}
        dimdomD['Version']    = 6
        dimdomD['Strucaero']  = 317
        dimdomD['NbDecimal']  = 8  
        dimdomD['Ndom']       = self.config['number_dom']
        if self.config['form'] == 'porous':
            dimdomD['Nvar']   = 8 # x,y,z,ro,rou,rov,row,p
        else:
            dimdomD['Nvar']   = 4 # x,y,z,p
        dimdomD['NbIns']      = self.config['number_dt'] # nbre d'instants
        dimdomD['NumVarMail'] = '1 2 3' # variables du maillage
        if self.config['form'] == 'solid':
            dimdomD['NumVarCham'] = '4' # pression est le 4eme champ
        else:
            dimdomD['NumVarCham'] = '4 5 6 7 8'
        dimdomD['TransX']     = 0
        dimdomD['TransY']     = 0
        dimdomD['TransZ']     = 0
        dimdomD['Perm']       = 0
        dimdomD['Rot']        = 0
        dimdomD['VarRef']     ='%f %f %f'%(self.config["density"],self.config["sound_velocity"],self.config["pressure_infinity"])
        if self.config['form_file'] == 'bin':
            dimdomD['Saut']   = '0 0 0'
        else:
            if self.config['form'] == 'solid': dimdomD['Saut'] = '10 0 0'
            else: dimdomD['Saut'] = '14 0 0'
        # Creation of dimdom lines
        dimdomLines = ''' Version du fichier dimdom
%(Version)d
 Structure des donnees aero. et nombre de decimales
%(Strucaero)d %(NbDecimal)d
 Nombre de domaines
%(Ndom)d
 Nombre total de variables des fichiers aero.
%(Nvar)d
 Nombre d"instants
%(NbIns)d
 Numeros des variables x, y et z
%(NumVarMail)s
 Numeros des variables aero. a copier
%(NumVarCham)s
 Coordonnees x, y et z d'un pt pour trans de l'origine
%(TransX)f %(TransY)f %(TransZ)f
 Nombre de perm xyz -> yzx et angle de rot / x
%(Perm)d %(Rot)f
References rho0, a0, p0
%(VarRef)s
 Sauts de lignes en ascii
%(Saut)s
idom nuda imax jmax kmax imi ima jmi  ...\n'''%(dimdomD)
#%(Dim)s\n'''%(dimdomD)
        if self.config['dom_id'] == 0:
            dimdom = open(self.kimEntreeDir+os.sep+'dimdom.in','w')
        else:
            dimdom = open(self.kimEntreeDir+os.sep+'dimdom'+str(self.config['dom_id'])+'.in','w')
        dimdom.write(dimdomLines)
        dimdom.close()

    def createKimAero(self):
        '''Create the KIM aero files'''

        # Lit "walls.cgns" pour avoir les dimensions uniquement
        tw = Cmpi.convertFile2SkeletonTree('stress.cgns')
        zonesw = Internal.getZones(tw)

        # write dim dom entetes des domaines
        dimdom = open(self.kimEntreeDir+os.sep+'dimdom.in','a')
        for c, z in enumerate(zonesw):
            nB_Pts = C.getNPts(z)
            nb_Elts = C.getNCells(z)
            dimz = Internal.getZoneDim(z)
            indom = c+1
            dim = [1,1,1]
            #dim[0] = 1
            #dim[1] = -nB_Pts
            #dim[2] = nb_Elts
            dim[0] = dimz[1]-1
            dim[1] = dimz[2]-1
            dim[2] = 1
            nuikim = 1
            invnorm = 1
            dimdom.write(str(indom)+' '+str(indom)+' '+str(dim[0])+' '+str(dim[1])+' '+str(dim[2])+'  1  '+ str(dim[0])+'  0  0  1  '+str(dim[2])+' '+str(nuikim)+' '+str(invnorm)+'\n')
        dimdom.close()

        # Lit le fichier des probes
        t = Cmpi.convertFile2SkeletonTree(self.config['filename'])
        zones = Internal.getZones(t)
        
        # Pour chaque domaine, load tous les instants
        checkZones = [] # Pour verification
        for iDom, zone in enumerate(zones):
            zonew = zonesw[iDom]
            dimz = Internal.getZoneDim(zonew)
            print(zone[0], zonew[0], flush=True)

            nodes = Internal.getNodesFromName(zone, 'FlowSolution#*')
            ncont = len(nodes) # nbre de containers pleins
            print('ncont=%d'%ncont); ncont = 1

            itglob = 1
            for nc in range(ncont):
                paths = ['CGNSTree/Base/%s/GridCoordinates#%d'%(zone[0],nc)]
                paths += ['CGNSTree/Base/%s/FlowSolution#%d'%(zone[0],nc)]
                nodes2 = Distributed.readNodesFromPaths(self.config['filename'], paths)
                px = Internal.getNodeFromName(nodes2[0], 'CoordinateX')[1]
                nrec = px.shape[0] # 50
                print('nrec=%d'%nrec); nrec = 1

                for it in range(0, nrec):
                    px = Internal.getNodeFromName(nodes2[0], 'CoordinateX')[1][it,:]
                    py = Internal.getNodeFromName(nodes2[0], 'CoordinateY')[1][it,:]
                    pz = Internal.getNodeFromName(nodes2[0], 'CoordinateZ')[1][it,:]
                    p =  Internal.getNodeFromName(nodes2[1], 'Pressure')[1][it,:]
            #        kimAeroFile = self.kimAeroDir+os.sep+'aerod%03dt%04d'%(iDom+self.config['dom_id'],itglob); itglob += 1
            #        fich = open(kimAeroFile, 'wb')                    
            #        self.write_fortran_block(fich, px, DOUBLE_SIZE) #X
            #        self.write_fortran_block(fich, py, DOUBLE_SIZE) #Y
            #        self.write_fortran_block(fich, pz, DOUBLE_SIZE) #Z
            #        self.write_fortran_block(fich, p, DOUBLE_SIZE) #P
            #        #for i in range(C.getNCells(data)):
            #        #    element = numpy.array([data[0][2][0][i],data[0][2][1][i],data[0][2][2][i]])
            #        #    self.write_fortran_block(fich,element, INT_SIZE, vector_type='i')#Connectivity
            #        fich.close()
                    
                    # Rebuild zone for check
                    cz = Internal.newZone(zone[0], zsize=[[dimz[1]-1,dimz[1]-2,0],[dimz[2]-1,dimz[2]-2,0]], ztype='Structured')
                    gc = Internal.newGridCoordinates(parent=cz)
                    ox = Internal.newDataArray('CoordinateX', value=px, parent=gc)
                    oy = Internal.newDataArray('CoordinateY', value=py, parent=gc)
                    oz = Internal.newDataArray('CoordinateZ', value=pz, parent=gc)                    
                    fs = Internal.newFlowSolution('FlowSolution', 'Vertex', parent=cz)
                    op = Internal.newDataArray('Pressure', value=p, parent=fs)
                    checkZones.append(cz)

        Internal.printTree(checkZones)
        C.convertPyTree2File(checkZones, 'out%d.cgns'%nc)