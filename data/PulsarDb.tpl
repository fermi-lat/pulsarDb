# Definition of Pulsar Ephemerides Database in FITS format (D4)
\include PulsarDb_primary.tpl
\include PulsarDb_spin_freq.tpl
\include PulsarDb_spin_per.tpl
\include PulsarDb_spin_hp.tpl
\include PulsarDb_orbital_dd.tpl
\include PulsarDb_orbital_bt.tpl
\include PulsarDb_orbital_ell1.tpl
\include PulsarDb_orbital_mss.tpl

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
EXTNAME      = 'REMARKS'                   / name of this binary table extension
TTYPE1       = 'PSRNAME'                   / pulsar name in PSR Jxxxx+xx[xx[aa]] format whenever available, or in any format otherwise
TFORM1       = '32A'                       / data format of field: character
TTYPE2       = 'EFFECTIVE_SINCE'           / time in MJD (TDB) since when remark on ephemeri(de)s is effective
TFORM2       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT2       = 'd'                         / physical unit of field
TTYPE3       = 'EFFECTIVE_UNTIL'           / time in MJD (TDB) until when remark on ephemeri(de)s is effective
TFORM3       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT3       = 'd'                         / physical unit of field
TTYPE4       = 'DESCRIPTION'               / description of remark on ephemeri(de)s
TFORM4       = '128A'                      / data format of field: character
TTYPE5       = 'OBSERVER_CODE'             / source of timing information
TFORM5       = '4A'                        / data format of field: character
END

XTENSION    = 'BINTABLE'                  / binary table extension
BITPIX      = 8                           / 8-bit bytes
NAXIS       = 2                           / 2-dimensional binary table
NAXIS1      =                             / width of table in bytes
NAXIS2      =                             / number of rows in table
PCOUNT      =                             / size of special data area
GCOUNT      = 1                           / one data group (required keyword)
TFIELDS     =                             / number of fields in each row
CHECKSUM    = ''                          / checksum for entire HDU
DATASUM     = ''                          / checksum for data table
TELESCOP    = 'GLAST'                     / name of telescope generating data
INSTRUME    = 'LAT'                       / name of instrument generating data
EQUINOX     = 2000.0                      / equinox for ra and dec
RADECSYS    = 'FK5'                       / world coord. system for this file (FK5 or FK4)
DATE        =                             / file creation date (YYYY-MM-DDThh:mm:ss UT)
EXTNAME     = 'OBSERVERS'                 / name of this binary table extension
TTYPE1      = 'OBSERVER_CODE'             / observer code
TFORM1      = '4A'                        / data format of field: character
TTYPE2      = 'OBSERVATORY'               / name of observatory or place of observatory
TFORM2      = '128A'                      / data format of field: character
TTYPE3      = 'CONTACT_PERSON'            / name of contact person
TFORM3      = '128A'                      / data format of field: character
TTYPE4      = 'REFERENCE'                 / reference for publications
TFORM4      = '1024A'                     / data format of field: character
END

XTENSION    = 'BINTABLE'                  / binary table extension
BITPIX      = 8                           / 8-bit bytes
NAXIS       = 2                           / 2-dimensional binary table
NAXIS1      =                             / width of table in bytes
NAXIS2      =                             / number of rows in table
PCOUNT      =                             / size of special data area
GCOUNT      = 1                           / one data group (required keyword)
TFIELDS     =                             / number of fields in each row
CHECKSUM    = ''                          / checksum for entire HDU
DATASUM     = ''                          / checksum for data table
TELESCOP    = 'GLAST'                     / name of telescope generating data
INSTRUME    = 'LAT'                       / name of instrument generating data
EQUINOX     = 2000.0                      / equinox for ra and dec
RADECSYS    = 'FK5'                       / world coord. system for this file (FK5 or FK4)
DATE        =                             / file creation date (YYYY-MM-DDThh:mm:ss UT)
EXTNAME     = 'ALTERNATIVE_NAMES'         / name of this binary table extension
TTYPE1      = 'ALTNAME'                   / alternative name for a pulsar
TFORM1      = '32A'                       / data format of field: character
TTYPE2      = 'PSRNAME'                   / pulsar name that appears in other extension(s)
TFORM2      = '32A'                       / data format of field: character
END
