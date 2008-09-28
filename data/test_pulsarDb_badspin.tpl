# FITS template to test PulsarDb class.
# This file represents a template in a wrong format for the following reason.
# --- EPHSTYLE header keyword is missing in SPIN_PARAMETERS extension.
SIMPLE      = T                                   / file does conform to FITS standard
BITPIX      = 8                                   / number of bits per data pixel
NAXIS       = 0                                   / number of data axes
EXTEND      = T                                   / FITS dataset may contain extensions
CHECKSUM    =                                     / checksum for entire HDU
DATASUM     =                                     / checksum for data table
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
CHECKSUM     =                             / checksum for entire HDU
DATASUM      =                             / checksum for data table
DATE         =                             / file creation date (YYYY-MM-DDThh:mm:ss UT)
EXTNAME      = 'SPIN_PARAMETERS'           / name of this binary table extension
TTYPE1       = 'STRING_VALUE'              / name of field
TFORM1       = '16A'                       / data format of field
END
