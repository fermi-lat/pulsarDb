/** \file EphStatus.h
    \brief Interface for EphStatus class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_EphStatus_h
#define pulsarDb_EphStatus_h

#include <limits>
#include <list>
#include <sstream>
#include <string>

#include "timeSystem/AbsoluteTime.h"

namespace timeSystem {
  template <typename TimeRepType>
  class TimeFormat;
}

namespace pulsarDb {

  enum EphStatusCodeType { Unavailable, Extrapolated, Remarked };

  /** \class EphStatus
      \brief Class representing status on pulsar ephemerides for a specific time interval, such as availability of ephemeris
             during the time interval, and a special remark on a glitch or excessive timing noise whose occurrence during the
             time interval was reported.
  */
  class EphStatus {
    public:
      /** \brief Create a ephemeris status object with given parameters.
          \param effective_since Beginning of time period during which this status is considered effective.
          \param effective_until End of time period during which this status is considered effective.
          \param description Description of the status on pulsar ephemerides.
      */
      EphStatus(const timeSystem::AbsoluteTime & effective_since, const timeSystem::AbsoluteTime & effective_until,
         const EphStatusCodeType & status_code, const std::string & description);

      /** \brief Return a start time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the start time itself.
      */
      const timeSystem::AbsoluteTime & getEffectiveSince() const;

      /** \brief Return an end time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the end time itself.
      */
      const timeSystem::AbsoluteTime & getEffectiveUntil() const;

      /// \brief Return a description of this status on ephemerides.
      const EphStatusCodeType & getStatusCode() const;

      /// \brief Return a description of this status on ephemerides.
      const std::string & getDescription() const;

      /** \brief Report this ephemeris status as a character string.
          \param time_system_name The time system to be used to represent end times of effectiveness window.
          \param time_format The time format to be used to represent end times of effectiveness window.
          \param precision Numerical precision to be used to represent end times of effectiveness window.
      */
      template <typename TimeRepType>
      std::string report(const std::string & time_system_name, const timeSystem::TimeFormat<TimeRepType> & time_format,
        std::streamsize precision = std::numeric_limits<double>::digits10) const;

      /** \brief Return true if this ephemeris status is effective in a certain part of a time interval between given times,
                 and return false otherwise.
          \param at1 One end of the time interval to be examined.
          \param at1 The other end of the time interval to be examined.
      */
      bool effectiveBetween(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const;

    private:
      timeSystem::AbsoluteTime m_since;
      timeSystem::AbsoluteTime m_until;
      EphStatusCodeType m_code;
      std::string m_description;
  };


  template <typename TimeRepType>
  std::string EphStatus::report(const std::string & time_system_name, const timeSystem::TimeFormat<TimeRepType> & time_format,
    std::streamsize precision) const {
    std::string report_string;

    // Create the start of a sentense.
    if (Unavailable == m_code) report_string += "No ephemeris available";
    else if (Extrapolated == m_code) report_string += "Ephemeris to be extrapolated";
    else if (Remarked == m_code) report_string += "Remarked";
    else report_string += "Unknown status reported";

    // Append a status description, if appropriate.
    if (Unavailable == m_code || Extrapolated == m_code) {
      if (!m_description.empty()) report_string += " (" + m_description + ")";
    } else {
      report_string += " \"" + m_description + "\"";
    }

    // Append a time part of the report.
    report_string += " since " + m_since.represent(time_system_name, time_format, precision) + " until " +
      m_until.represent(time_system_name, time_format, precision);

    // Return the report.
    return report_string;
  }

  typedef std::list<EphStatus> EphStatusCont;
}

#endif
