/**
    \mainpage pulsarDb package

    \author  Masaharu Hirayama hirayama@jca.umbc.edu
             James Peachey peachey@milkyway.gsfc.nasa.gov

    \section synopsis Synopsis 
    The purpose of this package is to provide an interface to the standard pulsar ephemeris
    database. It contains a library and two utilities, gtpulsardb and gtephcomp. The
    gtpulsardb utility is used to create, filter, and/or combine pulsar ephemerides database files.
    The gtephcomp utility is used to present the user with the best ephemeris available in a pulsar
    database for a given pulsar and instant of time.

    \subsection gtpulsardb_prerequisites gtpulsardb Prerequisites
    One or more pulsar ephemerides database files in either GLAST D4 FITS format, and/or
    a text format which will be described below.

    \subsection gtpulsardb_parameters gtpulsardb Parameters
 
\verbatim
psrdbfile [file name]
    Name of input file containing ephemerides. Multiple files
    may be combined by listing them in a text file, one per line,
    and supplying the list file name preceded by an @ sign.

outfile [file name]
    Name of output file, which will be in GLAST D4 FITS format.

filter = NONE [string]
    Type of filtering to be performed. Valid choices are
    NAME, TIME, or NONE. If filter is NAME, the pulsar name
    (see psrname parameter) will be used to select only ephemerides
    for the named pulsar. If filter is TIME, the parameters tstart and
    tstop will be used to filter the ephemerides based on that
    time range. If filter is NONE, no filtering will be performed.

psrname [string]
    The name of the pulsar, used to select only ephemerides
    valid for a particular pulsar. This only has effect
    when the filter parameter is NAME.

tstart [double]
    Time used for the beginning of the interval used for
    time filtering. This only has effect when the filter
    parameter is TIME.

tstop [double]
    Time used for the end of the interval used for
    time filtering. This only has effect when the filter
    parameter is TIME.

\endverbatim

    \subsection gtephcomp_prerequisites gtephcomp Prerequisites
    A single pulsar ephemerides database file in GLAST D4 FITS format.

    \subsection gtephcomp_parameters gtephcomp Parameters
 
\verbatim
psrname [string]
    The name of the pulsar, used to select only ephemerides
    valid for a particular pulsar. This only has effect
    when the filter parameter is NAME.

epoch [double]
    The time for which an ephemeris will be selected, if any is
    available in the input file. The interpretation of this number
    is determined by the timeformat and timesys parameters.

timeformat = MJD [string]
    String describing the representation used for the epoch.
    Valid choices are MJD and GLAST (MET).

timesys = TDB [string]
    String describing the time system used for the epoch.
    Valid choices are TDB and TT.

(psrdbfile = DEFAULT) [file name]
    Name of pulsar ephemerides database file, in GLAST D4
    FITS format. If psrdbfile is DEFAULT, the canonical pulsar
    database file (master_pulsardb.fits), which is distributed
    with the pulsar tools, will be used.

(strict = no) [bool]
    If strict is yes, only spin ephemerides whose stated
    range of validity contains the epoch will be selected.
    If strict is no, the ephemeris closest to the epoch
    will be selected, regardless of its stated range of
    validity.

\endverbatim

*/
