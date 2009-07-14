# -*- python -*-
# $Id: SConscript,v 1.18 2009/06/02 23:31:01 glastrm Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>
# Version: pulsarDb-08-04-00
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

progEnv.Tool('registerTargets', package = 'pulsarDb',
             staticLibraryCxts = [[pulsarDbLib, libEnv]],
             binaryCxts = [[gtephemBin, progEnv], [gtpulsardbBin, progEnv]],
             testAppCxts = [[test_pulsarDbBin, progEnv]],
             includes = listFiles(['pulsarDb/*.h']),
             pfiles = listFiles(['pfiles/*.par']),
             data = listFiles(['data/*'], recursive = True))
