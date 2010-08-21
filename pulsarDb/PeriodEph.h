/** \file PeriodEph.h
    \brief Interface for PeriodEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PeriodEph_h
#define pulsarDb_PeriodEph_h

#include <string>
#include <utility>

#include "pulsarDb/EphStatus.h"
#include "pulsarDb/PulsarEph.h"

#include "tip/Table.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/TimeInterval.h"
#include "timeSystem/TimeSystem.h"

namespace st_stream {
  class OStream;
}

namespace tip {
  class Header;
}

namespace pulsarDb {

  /** \class PeriodEph
      \brief Class representing a single pulsar ephemeris expressed with three period coefficients.
             Note: p0, p1, p2 depend on the time system, while AbsoluteTime objects don't.
  */
  class PeriodEph : public PulsarEph {
    public:
      /** \brief Create a pulsar ephemeris object with the given parameters.
          \param time_system_name Name of time system to interpret period parameters, such as "TDB" or "UTC".
          \param valid_since Beginning of time period during which this ephemeris is considered valid.
          \param valid_until End of time period during which this ephemeris is considered valid.
          \param epoch Reference epoch of frequency parameters (p0, p1, and p2).
          \param ra Right Ascension of the pulsar in degrees.
          \param dec Declination of the pulsar in degrees.
          \param phi0 Pulse phase at the given epoch (dimensionless).
          \param p0 Pulse period in seconds at the given epoch.
          \param p1 First time derivative of pulse period at the given epoch (dimensionless).
          \param p2 Second time derivative of pulse period in the units of s^(-1) at the given epoch.
      */
      PeriodEph(const std::string & time_system_name, const timeSystem::AbsoluteTime & valid_since,
        const timeSystem::AbsoluteTime & valid_until, const timeSystem::AbsoluteTime & epoch, double ra, double dec,
        double phi0, double p0, double p1, double p2): m_system(&timeSystem::TimeSystem::getSystem(time_system_name)),
        m_since(valid_since), m_until(valid_until), m_epoch(epoch), m_ra(ra), m_dec(dec), m_phi0(phi0), m_p0(p0), m_p1(p1), m_p2(p2) {}

      /** \brief Create a pulsar ephemeris object with the parameters stored in tip record.
          \param record FITS row from which all parameters for an ephemeris being created are to be read.
          \param header FITS header to read other information if necessary (not used).
      */
      PeriodEph(const tip::Table::ConstRecord & record, const tip::Header & header);

      /// \brief Destruct this PeriodEph object.
      virtual ~PeriodEph() {}

      /// \brief Return a time system to be used to interpret return values of calcFrequency method.
      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      /** \brief Return a start time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the start time itself.
      */
      virtual const timeSystem::AbsoluteTime & getValidSince() const { return m_since; }

      /** \brief Return an end time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the end time itself.
      */
      virtual const timeSystem::AbsoluteTime & getValidUntil() const { return m_until; }

      /// \brief Return a reference epoch of this ephemeris.
      virtual const timeSystem::AbsoluteTime & getEpoch() const { return m_epoch; }

      /// \brief Return the container of ephemeris remarks.
      virtual const EphStatusCont & getRemark() const { static const EphStatusCont s_remark_cont; return s_remark_cont; }

      /// \brief Create a copy of this object, and return a pointer to it.
      virtual PulsarEph * clone() const { return new PeriodEph(*this); }

      /** \brief Compute a pulse phase of a given time, and return it.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time for which a pulse phase is to be computed.
          \param phase_offset Value to be added to a pulse phase. This value is added to the computed pulse phase
                 before truncated to a value in range [0, 1).
      */
      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      /** \brief Compute a pulse frequency or its time derivative at a given time in the time system given to
                 this object upon its construction. Call getSystem method to obtain the time system to interpret
                 the return value of this method. The unit of frequency and its derivatives must be derived from
                 the time unit of seconds, i.e., Hz (sE-1), Hz/s (sE-2), and so on.
                 Note: validity of the ephemeris (valid since and valid until) are not checked. 
          \param ev_time Absolute time at which a pulse frequency is to be computed.
          \param derivative_order Order of frequency derivative to be computed. Set zero (0) to this argument
                 to compute a pulse frequency, one (1) the first time derivative of pulse frequency, and so on.
      */
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const;

      /** \brief Compute the pulsar position at a given time, and return it.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time at which the pulsar position is to be computed.
      */
      virtual timeSystem::SourcePosition calcPosition(const timeSystem::AbsoluteTime & ev_time) const;

    protected:
      /// \brief Output text expression of subclass-specific parameters of this object to a given output stream.
      virtual void writeModelParameter(st_stream::OStream & os) const;

    private:
      /** \brief Compute the number of elapsed seconds at a given absolute time since the reference epoch of this ephemeris,
                 measured in the time system in which this ephemeris is defined, and return it.
          \param at Absolute time at which the number of elapsed seconds is to be computed.
      */
      inline double calcElapsedSecond(const timeSystem::AbsoluteTime & at) const {
        return (at - m_epoch).computeDuration(m_system->getName(), "Sec");
      }

      const timeSystem::TimeSystem * m_system;
      timeSystem::AbsoluteTime m_since;
      timeSystem::AbsoluteTime m_until;
      timeSystem::AbsoluteTime m_epoch;
      double m_ra;
      double m_dec;
      double m_phi0;
      double m_p0;
      double m_p1;
      double m_p2;
  };

}

#endif
