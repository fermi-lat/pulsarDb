/** \file FormattedEph.h
    \brief Interface for FormattedEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_FormattedEph_h
#define pulsarDb_FormattedEph_h

#include <cmath>
#include <stdexcept>
#include <string>

#include "st_stream/Stream.h"
#include "tip/Table.h"

namespace tip {
  class Header;
}

namespace pulsarDb {

  /** \class ParameterFormatter
      \brief Class that formats parameter value listing in a text output. This class is desgined to be used for the shift
      operator (<<) for PulsarEph and OrbitalEph classes.
  */
  template <typename DataType>
  struct ParameterFormatter {
    ParameterFormatter(const std::string & param_name, const DataType & param_obj): m_name(param_name), m_obj(&param_obj) {}

    inline st_stream::OStream & write(st_stream::OStream & os) const {
      os.prefix().width(14); os << m_name + " = " << *m_obj;
      return os;
    }

    std::string m_name;
    const DataType * m_obj;
  };

  template <typename DataType>
  inline st_stream::OStream & operator <<(st_stream::OStream & os, const ParameterFormatter<DataType> & fmt) {
    return fmt.write(os);
  }

  /** \class FormattedEph
      \brief Class that provides helper functions for reading and writing pulsar ephemerides, both spin and orbital.
  */
  class FormattedEph {
    public:
      /** \brief Return a ParameterFormatter object to be used to format a text output of a given parameter.
          \param param_name Name of parameter to appear in a formatted text output.
          \param param_obj Object that holds the parameter value to output. The object must support a shift operator (<<)
                 for st_stream::OStream.
      */
      template <typename DataType>
      inline ParameterFormatter<DataType> format(const std::string & param_name, const DataType & param_obj) const {
        return ParameterFormatter<DataType>(param_name, param_obj);
      }

      /** \brief Helper method to get a value from a cell, returning it as the temmplated type, without handling the
          case where the value is null or a special value like Not-A-Number (throws an exception for such cases).
          \param record The record of tip::Table that contains the cell whose value to get.
          \param field_name The name of the field for the cell.
          \param data_value Variable to store the value.
      */
      template <typename DataType>
      void read(const tip::Table::ConstRecord & record, const std::string & field_name, DataType & data_value) {
        // Get the cell.
        const tip::TableCell & cell = record[field_name];

        // Try to get the cell content.
        if (cell.isNull()) throw std::runtime_error("tip::TableCell for field \"" + field_name + "\" is null");
        else cell.get(data_value);

        // Check the cell content for Not-A-Number type of data.
        if (IsNotAValue(data_value)) throw std::runtime_error("tip::TableCell for field \"" + field_name + "\" contains Not-A-Value");
      }

      /** \brief Helper method to get a value from a cell, returning it as the temmplated type, and handling the
          case where the cell doesn't exist, the value is null or a special value like Not-A-Number.
          \param record The record of tip::Table that contains the cell whose value to get.
          \param field_name The name of the field for the cell.
          \param data_value Variable to store the value.
          \param failed_value The value to set to data_value when exception occurs.
      */
      template <typename DataType>
      void read(const tip::Table::ConstRecord & record, const std::string & field_name, DataType & data_value,
        const DataType & failed_value) {
        try {
          read(record, field_name, data_value);
        } catch (...) {
          data_value = failed_value;
        }
      }

    private:
      /** \brief Helper method to check whether a given templated data is a valid value.
          \param x The data to be examined.
      */
      template <typename DataType>
      inline bool IsNotAValue(DataType /* x */) { return false; }

      /** \brief Specialized helper method to check whether a given double variable is Not-A-Number.
          \param x The data to be examined.
      */
      inline bool IsNotAValue(double x) {
#ifdef WIN32
        return 0 != _isnan(x);
#else
        return 0 != std::isnan(x);
#endif
      }
  };
}

#endif
