# -*- python -*-
# $Id: SConscript,v 1.28 2010/03/18 01:01:34 hirayama Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>
# Version: pulsarDb-08-05-06
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

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
