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
    /** \brief Construct a ParameterFormatter object.
        \param param_name Name of parameter to be printed in a formatted text.
        \param param_obj Object to be printed in a formatted text.
        \param separator Character string to be used as a separator between a parameter name and a parameter value in a formatted text.
    */
    ParameterFormatter(const std::string & param_name, const DataType & param_obj, const std::string & separator):
      m_name(param_name), m_obj(&param_obj), m_separator(separator) {}

    /** \brief Write a formatted text to an output stream.
        \param os Output stream to write a formatted text to.
    */
    inline st_stream::OStream & write(st_stream::OStream & os) const {
      os.prefix().width(16); os << m_name << m_separator << *m_obj;
      return os;
    }

    std::string m_name;
    const DataType * m_obj;
    std::string m_separator;
  };

  template <typename DataType>
  inline st_stream::OStream & operator <<(st_stream::OStream & os, const ParameterFormatter<DataType> & fmt) {
    return fmt.write(os);
  }

  /** \class FormattedEph
      \brief Class that provides helper functions for reading and writing pulsar ephemerides, both spin and orbital.
  */
  class FormattedEph {
    protected:
      /** \brief Return a ParameterFormatter object to be used to format a text output of a given parameter.
          \param param_name Name of parameter to appear in a formatted text output.
          \param param_obj Object that holds the parameter value to output. The object must support a shift operator (<<)
                 for st_stream::OStream.
      */
      template <typename DataType>
      inline ParameterFormatter<DataType> format(const std::string & param_name, const DataType & param_obj,
        const std::string & separator = " = ") const {
        return ParameterFormatter<DataType>(param_name, param_obj, separator);
      }

      /** \brief Helper method to get a value from a cell, returning it as the temmplated type, without handling the
          case where the value is null or a special value like Not-A-Number (throws an exception for such cases).
          \param record The record of tip::Table that contains the cell whose value to get.
          \param field_name The name of the field for the cell.
          \param data_value Variable to store the value.
      */
      template <typename DataType>
      void read(const tip::Table::ConstRecord & record, const std::string & field_name, DataType & data_value) const {
        // Get the cell.
        const tip::TableCell & cell = record[field_name];

        // Check whether this cell has a defined value.
        if (cell.isNull()) throw std::runtime_error("Field \"" + field_name + "\" is undefined");

        // Try to get the cell content.
        cell.get(data_value);

        // Throw an exception if the cell content is INDEF or not.
        if (isNan(data_value)) throw std::runtime_error("Value of field \"" + field_name + "\" is not a number");
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
        const DataType & failed_value) const {
        try {
          read(record, field_name, data_value);
        } catch (...) {
          data_value = failed_value;
        }
      }

      /** \brief Helper method that returns the fractional part of a given value, making sure that
          the return value is in the range of [0, 1).
          param phase_value Phase value whose fractional part is to be returned.
          \param phase_offset Phase value to be added to the computed pulse or orbital phase.
      */
      double trimPhaseValue(double phase_value, double phase_offset = 0.) const {
        double int_part; // ignored, needed for modf.
        double phase = std::modf(phase_value + phase_offset, &int_part);
        if (phase < 0.) ++phase;
        return phase;
      }

    private:
      /** \brief Helper method to check whether a given templated data is an INDEF.
          \param x The data to be examined.
      */
      template <typename DataType>
      inline bool isNan(DataType /* data_value */) const { return false; }

      inline bool isNan(float data_value) const { return isNotANumber(data_value); }
      inline bool isNan(double data_value) const { return isNotANumber(data_value); }

      /** \brief Helper method to check whether a given double variable is Not-A-Number.
          \param x The data to be examined.
      */
      template <typename T>
      inline bool isNotANumber(T x) const {
#ifdef WIN32
        return 0 != _isnan(x);
#else
        return 0 != std::isnan(x);
#endif
      }
  };
}

#endif
