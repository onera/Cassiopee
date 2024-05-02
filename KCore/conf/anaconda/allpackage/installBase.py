import socket
from os import environ as env


installDict = {
    socket.gethostname(): [
        'anaconda',
        env['FC'],  # f77compiler
        env['FC'],  # f90compiler
        env['CXX'],  # Cppcompiler
        env['CXXFLAGS'].split(),  # CppAdditionalOptions
        env['FFLAGS'].split(),  # f77AdditionalOptions
        False,  # useOMP
        False,  # static
        False,  # CPlotOffScreen
        env['CPPPATH'].split(),  # additionalIncludePaths
        #[],  # additionalLibs
        env['LIBPATH'].split(), # additionalLibs
        env['LIBPATH'].split(),  # additionalLibPaths
    ],
}
