/**
    \mainpage pulsarDb package

    \author  Masaharu Hirayama hirayama@jca.umbc.edu
             James Peachey James.Peachey-1@nasa.gov

    \section synopsis Synopsis 
    The purpose of this package is to provide an interface to the standard pulsar ephemeris
    database. It contains a library and two utilities, gtpulsardb and gtephem. The
    gtpulsardb utility is used to create, filter, and/or combine pulsar ephemerides database files.
    The gtephem utility is used to present the user with the best ephemeris available in a pulsar
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

filter = NONE [enumerated string (NAME|TIME|SOLAREPH|NONE)]
    Type of filtering to be performed. If filter is NAME, the pulsar name
    (see psrname parameter) will be used to select only ephemerides
    for the named pulsar. If filter is TIME, the parameters tstart and
    tstop will be used to filter the ephemerides based on that
    time range. If filter is SOLAREPH, only ephemerides that match the solareph
    parameter will be selected. If filter is NONE, no filtering will be performed.

psrname = ANY [string]
    Name of the pulsar, used to select only ephemerides valid for a
    particular pulsar. This only has effect when the filter parameter
    is NAME.

tstart = 0. [double]
    Time used for the beginning of the interval used for
    time filtering. This only has effect when the filter
    parameter is TIME.

tstop = 1.e5 [double]
    Time used for the end of the interval used for
    time filtering. This only has effect when the filter
    parameter is TIME.

solareph = JPL DE405 [enumerated string (JPL DE200|JPL DE405)]
    Solar system ephemeris used to filter the ephemerides.
    This only has effect when the filter parameter is SOLAREPH.

\endverbatim

    \subsection gtephem_prerequisites gtephem Prerequisites
    A single pulsar ephemerides database file in GLAST D4 FITS format.

    \subsection gtephem_parameters gtephem Parameters
 
\verbatim
psrdbfile [file name]
    Name of pulsar ephemerides database file, in GLAST D4
    FITS format.

psrname = ANY [string]
    Name of the pulsar, used to select only ephemerides valid for a
    particular pulsar.

reftime = 0. [string]
    Reference time for which an ephemeris will be selected, if any is
    available in the input file. The interpretation of this number is
    determined by the timeformat and timesys parameters.

timeformat = MJD [string]
    String describing the representation used for the reference time.
    Valid choices are MJD and GLAST (MET).

timesys = TDB [string]
    String describing the time system used for the reference time.
    Valid choices are TAI, TDB, TT and UTC.

(strict = no) [bool]
    If strict is yes, only spin ephemerides whose stated
    range of validity contains the epoch will be selected.
    If strict is no, the ephemeris closest to the epoch
    will be selected, regardless of its stated range of
    validity.

(solareph = JPL DE405) [enumerated string (JPL DE200|JPL DE405)]
    Solar system ephemeris used to select the ephemerides. This
    only has effect when the matchsolareph parameter is either
    ALL or PSRDB.

(matchsolareph = ALL) [enumerated string (NONE|PSRDB|ALL)]
    String that controls whether to use the name of the solar system
    ephemeris given by the solareph parameter to select ephemerides
    in the pulsar database. If matchsolareph is ALL or PSRDB, the string
    given by the solareph parameter is used to select the ephemerides.
    If matchsolareph is NONE, no selection will be performed.

(leapsecfile = DEFAULT) [file name]
    Name of the file containing the name of the leap second table, in
    OGIP-compliant leap second table format. If leapsecfile is the
    string DEFAULT, the default leap-second file (leapsec.fits), which
    is distributed with the extFiles package, will be used.

\endverbatim

    \section known-issues Known Issues
    None.

    \section resolved-issues Resolved Issues
\verbatim
Problem: Binary demodulation does not work on Windows. Some kind of
floating-precision problem prevents convergence. 

Resolution: Improved the handling of precision time expressions so as
not to rely on implementation-specific properties of the long double
type.
\endverbatim
*/
