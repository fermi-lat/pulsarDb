# -*- python -*-
# $Id: SConscript,v 1.24 2010/02/18 01:00:38 jrb Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>
# Version: pulsarDb-08-05-02
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
