import shutil
import os
import pathlib

##NOTE:: Should only be used for SATOR
directory_cbenoit = '/stck/benoit/Cassiopee/Dist/bin'
directory_local   = os.environ['CASSIOPEE']
elsaprod          = os.environ['ELSAPROD']
directory_cbenoit += '/'+elsaprod
directory_local   += '/Dist/bin/'+elsaprod

file_path = pathlib.Path(directory_cbenoit)
if file_path.exists():
    print(elsaprod+' exists.')
    onlyfiles = [f for f in os.listdir(directory_cbenoit+'/lib/') if os.path.isfile(os.path.join(directory_cbenoit+'/lib/', f))]
    file_list_lib=[]
    for file_tmp in onlyfiles:
        if file_tmp.split('.')[1] =='so':
            file_list_lib.append(file_tmp)
            print('lib: ',file_tmp)
    del onlyfiles


    file_list_include=['png.h',
                       'pngconf.h',
                       'pnglibconf.h']
    list_folders     =['GL','X11']


    ##Copy lib
    for file_tmp in file_list_lib:
        shutil.copyfile(directory_cbenoit+'/lib/'+file_tmp,directory_local+'/lib/'+file_tmp)

    ##Copy include
    file2include  = directory_local+'/include/'
    file_path_tmp = pathlib.Path(file2include)
    if file_path_tmp.exists():
        for f_tmp in list_folders:
            if pathlib.Path(file2include+f_tmp).exists():
                shutil.rmtree(file2include+f_tmp)

            if pathlib.Path(directory_cbenoit+'/include/'+f_tmp).exists():
                print('include:',f_tmp)
                shutil.copytree(directory_cbenoit+'/include/'+f_tmp ,file2include+f_tmp)

        for file_tmp in file_list_include:
            if os.path.exists(directory_cbenoit+'/include/'+file_tmp):
                print('include:',file_tmp)
                shutil.copyfile(directory_cbenoit+'/include/'+file_tmp,file2include+file_tmp)
    else:
        shutil.copytree(directory_cbenoit+'/include/',directory_local+'/include/')
else:
    print(elsaprod+'=does not exist')
