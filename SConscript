#$Id: SConscript,v 1.3 2008/02/26 05:27:32 glastrm Exp $
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('pulsarDbLib', depsOnly = 1)
pulsarDbLib = libEnv.StaticLibrary('pulsarDb', listFiles(['src/*.cxx']))

progEnv.Tool('pulsarDbLib')
gtephemBin = progEnv.Program('gtephem', listFiles(['src/gtephem/*.cxx']))
gtpulsardbBin = progEnv.Program('gtpulsardb', listFiles(['src/gtpulsardb/*.cxx']))
test_pulsarDbBin = progEnv.Program('test_pulsarDb', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', package = 'pulsarDb', libraries = [pulsarDbLib], binaries = [gtephemBin, gtpulsardbBin], testApps = [test_pulsarDbBin], includes = listFiles(['pulsarDb/*.h']),
             pfiles = listFiles(['pfiles/*.par']), data = listFiles(['data/*'], recursive = True))
