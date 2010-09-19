/** \file FormattedEph.h
    \brief Interface for FormattedEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_FormattedEph_h
#define pulsarDb_FormattedEph_h

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "st_stream/Stream.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

namespace tip {
  class Header;
}

namespace pulsarDb {

  /** \class ParameterFormatter
      \brief Class that formats parameter value listing in a text output. This class is designed to be used for the shift
      operator (<<) for PulsarEph and OrbitalEph classes.
  */
  template <typename DataType>
  struct ParameterFormatter {
    /** \brief Construct a ParameterFormatter object.
        \param param_name Name of parameter to be printed in a formatted text.
        \param param_obj Object to be printed in a formatted text.
        \param param_unit Character string that represents the physical unit for the parameter value.
        \param separator Character string to be used as a separator between a parameter name and a parameter value in a formatted text.
    */
    ParameterFormatter(const std::string & param_name, const DataType & param_obj, const std::string & param_unit,
      const std::string & separator): m_name(param_name), m_obj(&param_obj), m_unit(param_unit), m_separator(separator) {}

    /** \brief Write a formatted text to an output stream.
        \param os Output stream to write a formatted text to.
    */
    inline st_stream::OStream & write(st_stream::OStream & os) const {
      os.prefix().width(16); os << m_name << m_separator << *m_obj;
      if (!m_unit.empty()) os << " " << m_unit;
      return os;
    }

    std::string m_name;
    const DataType * m_obj;
    std::string m_unit;
    std::string m_separator;
  };

  template <typename DataType>
  inline st_stream::OStream & operator <<(st_stream::OStream & os, const ParameterFormatter<DataType> & fmt) {
    return fmt.write(os);
  }

  /** \class FormattedEph
      \brief Class that provides virtual methods and helper functions for reading and writing pulsar ephemerides,
             both spin and orbital.
  */
  class FormattedEph {
    public:
      /// \brief Output text expression of this FormattedEph to a given output stream.
      virtual st_stream::OStream & write(st_stream::OStream & os) const = 0;

    protected:
      /** \class CellReadError
          \brief Class that indicates an error in reading a FITS cell.
      */
      class CellReadError : public std::runtime_error {
        public:
          /** \brief Construct a CellReadError object by calling a base class constructor.
              \param what_arg Character string that describes the error.
          */
          explicit CellReadError(const std::string & what_arg): std::runtime_error(what_arg) {}
      };

      /** \brief Return a ParameterFormatter object to be used to format a text output of a given parameter.
          \param param_name Name of parameter to appear in a formatted text output.
          \param param_obj Object that holds the parameter value to output. The object must support a shift operator (<<)
                 for st_stream::OStream.
          \param param_unit Character string that represents the physical unit for the parameter value.
      */
      template <typename DataType>
      inline ParameterFormatter<DataType> format(const std::string & param_name, const DataType & param_obj,
        const std::string & param_unit, const std::string & separator = " = ") const {
        return ParameterFormatter<DataType>(param_name, param_obj, param_unit, separator);
      }

      /** \brief Helper method to get a value from a cell, returning it through an output argument (data_value).
                 This method throws an exception if the cell contains a special value like a NULL and a Not-A-Number.
                 It also throws an exception for an error in reading a FITS cell, such as a named column does not exist
                 in the FITS table that the given cell belongs to. If it does not throw an exception, it returns a logical
                 false as a return value of the method, indicating no error has occurred in reading out the cell.
          \param record Record of tip::Table that contains the cell to read its value from.
          \param field_name Name of the field for the cell.
          \param data_value Output argument to which the value read from the cell is set.
      */
      template <typename DataType>
      bool read(const tip::Table::ConstRecord & record, const std::string & field_name, DataType & data_value) const;

      /** \brief Helper method to get a value from a cell, returning it through an output argument (data_value).
                 If the cell contains a NULL or a Not-A-Number, or if a named cell does not exist in the given
                 FITS record (record), then the given default value (default_value) is set to the output argument,
                 and a logical true is returned to indicate that. This method throws an exception for other errors
                 in reading a FITS cell. In all other cases, it returns a logical false to indicate no error has
                 occurred in reading out the cell.
          \param record Record of tip::Table that contains the cell to read its value from.
          \param field_name Name of the field for the cell.
          \param data_value Output argument to which the value read from the cell is set.
          \param default_value Default value to set to the output argument (data_value) when the named column does
                 not exist, the content is a NULL or a Not-A-Number.
      */
      template <typename DataType>
      bool read(const tip::Table::ConstRecord & record, const std::string & field_name, DataType & data_value,
        const DataType & default_value) const;

      /** \brief Helper method to get a data array from a cell, returning it through an output argument (data_array).
                 If the cell contains a NULL or a Not-A-Number, then the given default value (default_value) is set to
                 the corresponding element of the output argument, the corresponding element of an output Boolean array
                 (is_default) is set, and a logical OR of the Boolean array is returned. If a named cell does not exist
                 in the given FITS record (record), then this method removes all elements from the output argument and
                 the output Boolean array (which make them zero-length arrays), and returns a logical true to indicate
                 the error (i.e., no such column exists). This method throws an exception for other errors in reading
                 a FITS cell. In all other cases, this method returns a logical false, indicating no error has occurred
                 in reading out the cell.
          \param record Record of tip::Table that contains the cell to read its value from.
          \param field_name Name of the field for the cell.
          \param data_array Output argument to which the values read from the cell are set.
          \param default_value Default value to set to an element of the output argument (data_array) when the cell
                 contains a NULL or a Not-A-Number.
          \param is_default Output Boolean array which indicates that the given default value is set to the corresponding
                 element of the output argument (data_array).
      */
      template <typename DataType>
      bool read(const tip::Table::ConstRecord & record, const std::string & field_name, std::vector<DataType> & data_array,
        const DataType & default_value, std::vector<bool> & is_default) const;

      // TODO: Split the above method into the following two methods.
      //template <typename DataType>
      //bool read(const tip::Table::ConstRecord & record, const std::string & field_name, std::vector<DataType> & data_array,
      //  const DataType & default_value) const;
      //template <typename DataType>
      //bool read(const tip::Table::ConstRecord & record, const std::string & field_name,
      //  std::vector<std::pair<DataType, bool> > & data_array) const;

      /** \brief Helper method that returns the fractional part of a given value, making sure that
                 the return value is in the range of [0, 1).
          \param phase_value Phase value whose fractional part is to be returned.
          \param phase_offset Phase value to be added to the computed pulse or orbital phase.
      */
      double trimPhaseValue(double phase_value, double phase_offset = 0.) const;

    private:
      /** \brief Helper method to check whether a given templated data is an INDEF.
          \param data_value Data to be examined.
      */
      template <typename DataType>
      inline bool isNan(DataType /* data_value */) const { return false; }
      inline bool isNan(float data_value) const { return isNotANumber(data_value); }
      inline bool isNan(double data_value) const { return isNotANumber(data_value); }

      /** \brief Helper method to check whether a given floating-point number is Not-A-Number.
          \param data_value Data to be examined.
      */
      template <typename DataType>
      bool isNotANumber(DataType data_value) const;

      /** \brief Helper method to check whether elements of a given data array are Not-A-Number.
                 The results are stored in an output array of Boolean values (nan_array), and this
                 method returns a logical OR of all the Boolean values.
          \param data_array Array of data to examine.
          \param nan_array Output array to which a logical true is set if a corresponding data element
                 is Not-A-Number.
      */
      template <typename DataType> bool hasNan(DataType /* data_array */, std::vector<bool> & /* nan_array */) const;
      template <typename DataType> bool hasNan(std::vector<DataType> data_array, std::vector<bool> & nan_array) const;
  };

  template <typename DataType>
  bool FormattedEph::read(const tip::Table::ConstRecord & record, const std::string & field_name, DataType & data_value) const {
    // Get the field index and check whether it is a scalar column or a vector.
    tip::FieldIndex_t field_index = 0;
    try {
      field_index = record.getExtensionData()->getFieldIndex(field_name);
    } catch (const tip::TipException & x) {
      // Re-throw the error as CellReadError, so that the caller can distinguish it from other exceptions.
      throw CellReadError(x.what());
    }
    bool is_scalar = record.getExtensionData()->getColumn(field_index)->isScalar();

    // Get the cell.
    const tip::TableCell & cell = record[field_name];

    // Check whether this cell contains a defined value.
    if (is_scalar) {
      if (cell.isNull()) throw CellReadError("Field \"" + field_name + "\" is undefined");
    } else {
      std::vector<bool> null_array;
      bool has_null = record.getExtensionData()->getColumn(field_index)->getNull(record.getIndex(), null_array);
      if (has_null) throw CellReadError("Field \"" + field_name + "\" has an undefined element");
      // TODO: Replace the above with the below once tip::TableCell::getNull is implemented.
      //if (cell.getNull(null_array)) throw CellReadError("Field \"" + field_name + "\" has an undefined element");
    }

    // Try to get the cell content.
    cell.get(data_value);

    // Throw an exception if the cell contains an Not-A-Number.
    if (is_scalar) {
      if (isNan(data_value)) throw CellReadError("Value of field \"" + field_name + "\" is not a number");
    } else {
      std::vector<bool> nan_array;
      if (hasNan(data_value, nan_array)) throw CellReadError("Value(s) of field \"" + field_name + "\" are not a number");
    }

    // Always return a logical false.
    return false;
  }

  template <typename DataType>
  bool FormattedEph::read(const tip::Table::ConstRecord & record, const std::string & field_name, DataType & data_value,
    const DataType & default_value) const {
    // Prepare a return value.
    bool is_default = false;

    // Try to read out the cell.
    try {
      read(record, field_name, data_value);
    } catch (const CellReadError &) {
      data_value = default_value;
      is_default = true;
    }

    // Return the flag.
    return is_default;
  }

  template <typename DataType>
  bool FormattedEph::read(const tip::Table::ConstRecord & record, const std::string & field_name, std::vector<DataType> & data_array,
    const DataType & default_value, std::vector<bool> & is_default) const {
    // Get the field index.
    bool no_such_field = false;
    tip::FieldIndex_t field_index = 0;
    try {
      field_index = record.getExtensionData()->getFieldIndex(field_name);
    } catch (const tip::TipException & x) {
      no_such_field = true;
    }
    if (no_such_field) {
      data_array.clear();
      is_default.clear();
      return true;
    }

    // Throw an exception if it is a scalar column.
    if (record.getExtensionData()->getColumn(field_index)->isScalar()) {
      throw CellReadError("Attempted to read an array of data from a scalar column");
    }

    // Get the cell.
    const tip::TableCell & cell = record[field_name];

    // Try to get the cell content, and prepare for the output flags.
    cell.get(data_array);
    typename std::vector<DataType>::size_type array_size = data_array.size();
    is_default.resize(array_size);
    if (is_default.size() != array_size) throw std::runtime_error("Failed to allocate memory space for default flags");

    // Check whether this cell contains a defined value.
    std::vector<bool> null_array;
    record.getExtensionData()->getColumn(field_index)->getNull(record.getIndex(), null_array);
    // TODO: Replace the above with the below once tip::TableCell::getNull is implemented.
    //cell.getNull(null_array));
    if (null_array.size() != array_size) throw std::runtime_error("Failed to allocate memory space for NULL flags");

    // Check whether the cell contains an Not-A-Number.
    std::vector<bool> nan_array;
    hasNan(data_array, nan_array);
    if (nan_array.size() != array_size) throw std::runtime_error("Failed to allocate memory space for Not-A-Number flags");

    // Take a logical OR of NULL and Not-A-Number, and set the default value to a corresponding element.
    bool return_value = false;
    for (std::vector<bool>::size_type ii = 0; ii < array_size; ++ii) {
      is_default[ii] = null_array[ii] || nan_array[ii];
      if (is_default[ii]) {
        data_array[ii] = default_value;
        return_value = true;
      }
    }

    // Return the logical OR of all NULL-or-Not-A-Number flags.
    return return_value;
  }

  inline double FormattedEph::trimPhaseValue(double phase_value, double phase_offset) const {
    double int_part; // ignored, needed for modf.
    double phase = std::modf(phase_value + phase_offset, &int_part);
    if (phase < 0.) ++phase;
    return phase;
  }

  template <typename DataType>
  inline bool FormattedEph::isNotANumber(DataType data_value) const {
#ifdef WIN32
    return 0 != _isnan(data_value);
#else
    return 0 != std::isnan(data_value);
#endif
  }

  template <typename DataType>
  bool FormattedEph::hasNan(DataType /* data_array */, std::vector<bool> & /* nan_array */) const {
    throw std::runtime_error("Unsupported array type is used for reading a vector column");
  }

  template <typename DataType>
  bool FormattedEph::hasNan(std::vector<DataType> data_array, std::vector<bool> & nan_array) const {
    // Allocate memory space for Boolean flags.
    typename std::vector<DataType>::size_type array_size = data_array.size();
    nan_array.resize(array_size);
    if (nan_array.size() != array_size) throw std::runtime_error("Failed to allocate memory space for Not-A-Number flags");

    // Check the array contents.
    bool has_nan = false;
    for (typename std::vector<DataType>::size_type ii = 0; ii < array_size; ++ii) {
      nan_array[ii] = isNan(data_array[ii]);
      has_nan |= nan_array[ii];
    }

    // Return the logical OR of the results.
    return has_nan;
  }

  inline st_stream::OStream & operator <<(st_stream::OStream & os, const FormattedEph & eph) { return eph.write(os); }

}

#endif
