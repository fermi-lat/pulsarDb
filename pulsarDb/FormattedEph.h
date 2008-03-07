/** \file FormattedEph.h
    \brief Interface for FormattedEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_FormattedEph_h
#define pulsarDb_FormattedEph_h

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
      ParameterFormatter<DataType> format(const std::string & param_name, const DataType & param_obj) const {
        return ParameterFormatter<DataType>(param_name, param_obj);
      }

      /** \brief Helper method to get a value from a cell, returning it as a double, and handling the
          case where the value is null.
          \param cell The cell whose value to get.
      */
      inline double get(const tip::TableCell & cell) const {
        double value = 0.;
        get(cell, value);
        return value;
      }

      /** \brief Helper method to get a value from a cell, returning it as the temmplated type, and handling the
          case where the value is null.
          \param cell The cell whose value to get.
          \param value Variable to store the value.
      */
      template <typename T>
      inline void get(const tip::TableCell & cell, T & value) const {
        // WARNING: This will break for a string column.
        if (cell.isNull()) value = 0;
        else cell.get(value);
      }
  };
}

#endif
