import glob,os,platform

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

pulsarDbLib = libEnv.StaticLibrary('pulsarDb', listFiles(['src/*.cxx']))

progEnv.Tool('pulsarDbLib')
gtephcompBin = progEnv.Program('gtephcomp', listFiles(['src/gtephcomp/*.cxx']))
gtpulsardbBin = progEnv.Program('gtpulsardb', listFiles(['src/gtpulsardb/*.cxx']))

progEnv.Tool('registerObjects', package = 'pulsarDb', libraries = [pulsarDbLib], binaries = [gtephcompBin, gtpulsardbBin], includes = listFiles(['pulsarDb/*.h']), pfiles = listFiles(['pfiles/*.par']))