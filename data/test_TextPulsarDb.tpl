# FITS template to test the text parser of PulsarDb class.
SIMPLE      = T                                   / file does conform to FITS standard
BITPIX      = 8                                   / number of bits per data pixel
NAXIS       = 0                                   / number of data axes
EXTEND      = T                                   / FITS dataset may contain extensions
CHECKSUM    = ''                                  / checksum for entire HDU
DATASUM     = ''                                  / checksum for data table
DATE        =                                     / file creation date (YYYY-MM-DDThh:mm:ss UT)
FILENAME    =                                     / name of this file
ORIGIN      = 'LISOC'                             / name of organization making file
AUTHOR      =                                     / name of person responsible for file generation
CREATOR     =                                     / software and version creating file
VERSION     = 1                                   / release version of the file
END

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
DATE         =                             / file creation date (YYYY-MM-DDThh:mm:ss UT)
EXTNAME      = 'SPIN_PARAMETERS'           / name of this binary table extension
EPHSTYLE     = 'TEST'                      / name of pulsar ephemeris model
PDBTGEN      = 1                           / first generation table of pulsar ephemerides database
TTYPE1       = 'PSRNAME'                   / name of field
TFORM1       = '32A'                       / data format of field
TTYPE2       = '1J_COLUMN'                 / name of field
TFORM2       = '1J'                        / data format of field
TTYPE3       = '1D_COLUMN'                 / name of field
TFORM3       = '1D'                        / data format of field
TTYPE4       = '3J_COLUMN'                 / name of field
TFORM4       = '3J'                        / data format of field
TTYPE5       = '3D_COLUMN'                 / name of field
TFORM5       = '3D'                        / data format of field
TTYPE6       = 'PJ_COLUMN'                 / name of field
TFORM6       = 'PJ'                        / data format of field
TTYPE7       = 'PD_COLUMN'                 / name of field
TFORM7       = 'PD'                        / data format of field
END
