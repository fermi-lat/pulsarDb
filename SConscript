# -*- python -*-
# $Id: SConscript,v 1.12 2008/10/25 02:30:51 glastrm Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>
# Version: pulsarDb-08-01-01
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
