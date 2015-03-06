# Extension of Pulsar Ephemerides Database to store orbital ephemerides of the simplified DD model
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
EXTNAME      = 'ORBITAL_PARAMETERS'        / name of this binary table extension
EPHSTYLE     = 'ELL1'                      / name of binary orbital model
PDBTGEN      = 2                           / second generation table of pulsar ephemerides database
TTYPE1       = 'PSRNAME'                   / pulsar name in PSR Jxxxx+xx[xx[aa]] format whenever available, or in any format otherwise
TFORM1       = '32A'                       / data format of field: character
TTYPE2       = 'PB'                        / orbital period
TFORM2       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT2       = 's'                         / physical unit of field
TTYPE3       = 'PBDOT'                     / first time derivative of PB (orbital period)
TFORM3       = 'D'                         / data format of field: 4-byte REAL
TUNIT3       = ''                          / physical unit of field: dimensionless
TTYPE4       = 'A1'                        / projected semi-major axis in light seconds
TFORM4       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT4       = 'lt-s'                      / physical unit of field
TTYPE5       = 'XDOT'                      / first time derivative of A1 (projected semi-major axis)
TFORM5       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT5       = 'lt-s / s'                  / physical unit of field
TTYPE6       = 'EPS1'                      / eccentricity multiplied by the sine of the periastron longitude
TFORM6       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT6       = ''                          / physical unit of field: dimensionless
TTYPE7       = 'EPS1DOT'                   / first time derivative of EPS1
TFORM7       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT7       = 's**(-1)'                   / physical unit of field
TTYPE8       = 'EPS2'                      / eccentricity multiplied by the cosine of the periastron longitude
TFORM8       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT8       = ''                          / physical unit of field: dimensionless
TTYPE9       = 'EPS2DOT'                   / first time derivative of EPS2
TFORM9       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT9       = 's**(-1)'                   / physical unit of field
TTYPE10      = 'TASC'                      / time of ascending node in MJD (TDB)
TFORM10      = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT10      = 'd'                         / physical unit of field
TTYPE11      = 'GAMMA'                     / time-dilation and gravitational redshift parameter
TFORM11      = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT11      = 's'                         / physical unit of field
TTYPE12      = 'SHAPIRO_R'                 / range parameter of Shapiro delay in binary system
TFORM12      = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT12      = 'us'                        / physical unit of field
TTYPE13      = 'SHAPIRO_S'                 / shape parameter of Shapiro delay in binary system
TFORM13      = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT13      = ''                          / physical unit of field: dimensionless
TTYPE14      = 'OBSERVER_CODE'             / source of orbital parameter information
TFORM14      = '4A'                        / data format of field: character
TTYPE15      = 'SOLAR_SYSTEM_EPHEMERIS'    / name of solar system ephemeris used for barycentric corrections ("JPL DE200" or "JPL DE405")
TFORM15      = '32A'                       / data format of field: character
END
