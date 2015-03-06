# Extension of Pulsar Ephemerides Database to store spin ephemerides of a simple period history model
XTENSION     = 'BINTABLE'                  / binary table extension
BITPIX       = 8                           / 8-bit bytes
NAXIS        = 2                           / 2-dimensional binary table
NAXIS1       =                             / width of table in bytes
NAXIS2       =                             / number of rows in table
PCOUNT       =                             / size of special data area
GCOUNT       = 1                           / one data group (required keyword)
TFIELDS      =                             / number of fields in each row
CHECKSUM     = ''                          / checksum for entire HDU
DATASUM      = ''                          / checksum for data table
TELESCOP     = 'GLAST'                     / name of telescope generating data
INSTRUME     = 'LAT'                       / name of instrument generating data
EQUINOX      = 2000.0                      / equinox for ra and dec
RADECSYS     = 'FK5'                       / world coord. system for this file (FK5 or FK4)
DATE         =                             / file creation date (YYYY-MM-DDThh:mm:ss UT)
EXTNAME      = 'SPIN_PARAMETERS'           / name of this binary table extension
EPHSTYLE     = 'PER'                       / name of pulsar ephemeris model
PDBTGEN      = 2                           / second generation table of pulsar ephemerides database
TTYPE1       = 'PSRNAME'                   / pulsar name in PSR Jxxxx+xx[xx[aa]] format whenever available, or in any format otherwise
TFORM1       = '32A'                       / data format of field: character
TTYPE2       = 'RA'                        / right ascension (J2000) of the pulsar position
TFORM2       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT2       = 'deg'                       / physical unit of field
TLMIN2       = 0.0                         / minimum value
TLMAX2       = 360.0                       / maximum value
TTYPE3       = 'DEC'                       / declination (J2000) of the pulsar position
TFORM3       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT3       = 'deg'                       / physical unit of field
TLMIN3       = -90.0                       / minimum value
TLMAX3       = 90.0                        / maximum value
TTYPE4       = 'VALID_SINCE'               / first date for valid timing parameters in MJD
TFORM4       = 'J'                         / data format of field: 4-byte signed INTEGER
TUNIT4       = 'd'                         / physical unit of field
TTYPE5       = 'VALID_UNTIL'               / last date for valid timing parameters in MJD
TFORM5       = 'J'                         / data format of field: 4-byte signed INTEGER
TUNIT5       = 'd'                         / physical unit of field
TTYPE6       = 'EPOCH_INT'                 / integer part of the epoch for P0, P1, and P2 in MJD (TDB)
TFORM6       = 'J'                         / data format of field: 4-byte signed INTEGER
TUNIT6       = 'd'                         / physical unit of field
TTYPE7       = 'EPOCH_FRAC'                / fractional part of the epoch for P0, P1, and P2 in MJD (TDB)
TFORM7       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT7       = 'd'                         / physical unit of field
TLMIN7       = 0.0                         / minimum value
TLMAX7       = 1.0                         / maximum value
TTYPE8       = 'TOAGEO_INT'                / integer part of infinite-frequency geocentric pulse arrival time in MJD (TT)
TFORM8       = 'J'                         / data format of field: 4-byte signed INTEGER
TUNIT8       = 'd'                         / physical unit of field
TTYPE9       = 'TOAGEO_FRAC'               / fractional part of infinite-frequency geocentric pulse arrival time in MJD (TT)
TFORM9       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT9       = 'd'                         / physical unit of field
TLMIN9       = 0.0                         / minimum value
TLMAX9       = 1.1                         / maximum value
TTYPE10      = 'TOABARY_INT'               / integer part of infinite-frequency, binary-demodulated, barycentric pulse arrival time in MJD (TDB)
TFORM10      = 'J'                         / data format of field: 4-byte signed INTEGER
TUNIT10      = 'd'                         / physical unit of field
TTYPE11      = 'TOABARY_FRAC'              / fractional part of infinite-frequency, binary-demodulated, barycentric pulse arrival time in MJD (TDB)
TFORM11      = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT11      = 'd'                         / physical unit of field
TLMIN11      = 0.0                         / minimum value
TLMAX11      = 1.0                         / maximum value
TTYPE12      = 'P0'                        / pulsar rotation period
TFORM12      = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT12      = 's'                         / physical unit of field
TTYPE13      = 'P1'                        / first time derivative of pulsar period
TFORM13      = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT13      = ''                          / physical unit of field: dimensionless
TTYPE14      = 'P2'                        / second time derivative of pulsar period
TFORM14      = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT14      = 's**(-1)'                   / physical unit of field
TTYPE15      = 'RMS'                       / root-mean-square radio timing residual in milli-periods
TFORM15      = 'E'                         / data format of field: 4-byte REAL
TUNIT15      = ''                          / physical unit of field: dimensionless
TLMIN15      = 0.0                         / minimum value
TLMAX15      = 100000.0                    / maximum value
TTYPE16      = 'OBSERVER_CODE'             / source of timing information
TFORM16      = '4A'                        / data format of field: character
TTYPE17      = 'BINARY_FLAG'               / true for binary pulsars, false for single pulsars
TFORM17      = 'L'                         / data format of field: logical
TTYPE18      = 'SOLAR_SYSTEM_EPHEMERIS'    / name of solar system ephemeris used for barycentric corrections ("JPL DE200" or "JPL DE405")
TFORM18      = '32A'                       / data format of field: character
END
