/** \file PdotCanceler.h
    \brief Interface for PdotCanceler class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PdotCanceler_h
#define pulsarDb_PdotCanceler_h

#include <string>
#include <vector>

#include "timeSystem/AbsoluteTime.h"

namespace timeSystem {
  class TimeSystem;
}

namespace pulsarDb {

  class PulsarEph;

  /** \class PdotCanceler
      \brief Class representing a single pulsar ephemeris. Warning: f0, f1, f2 depend on the time system.
      While AbsoluteTime objects are interchangeable, these other values are not!
  */
  class PdotCanceler {
    public:
      /** \brief Construct a PdotCanceler object.
          \param time_system_name Name of time system in which it performs pdot cancellation.
          \param time_origin Origin of time for the pdot cancellation.
          \param fdot_ratio Array of ratios of frequency derivatives over frequency at the origin of time.
                 Its first element is ratio of f1 over f0, the second is ratio of f2 over f0, and so on,
                 where f0 is the frequency at the origin of time, f1 the first time-derivative of frequency
                 at the origin of time, f2 the second, and so on.
      */
      PdotCanceler(const std::string & time_system_name, const timeSystem::AbsoluteTime & time_origin,
        const std::vector<double> & fdot_ratio);

      /** \brief Construct a PdotCanceler object.
          \param pulsar_eph Pulsar's spin ephemeris from which frequency and its time-derivatives are computed.
          \param max_derivative Maximum order of time derivatives of frequency to consider. Setting 0 (zero) to
                 this parameter implies no pdot cancellation will be performed. Setting 2 (two) to this parameter,
                 the first and the second time-derivatives of frequency will be computed and used to perform pdot
                 cancellations.
      */
      PdotCanceler(const timeSystem::AbsoluteTime & time_origin, const PulsarEph & pulsar_eph, int max_derivative);

      virtual ~PdotCanceler() {}

      /** \brief Correct an absolute time to account for pdot cancellation.
          \param abs_time Time to correct.
      */
      virtual void cancelPdot(timeSystem::AbsoluteTime & abs_time) const;

    protected:
      const timeSystem::TimeSystem * m_time_system;
      timeSystem::AbsoluteTime m_time_origin;
      std::vector<double>  m_fdot_ratio;
  };

}

#endif
