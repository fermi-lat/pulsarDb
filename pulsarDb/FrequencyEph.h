/** \file FrequencyEph.h
    \brief Interface for FrequencyEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_FrequencyEph_h
#define pulsarDb_FrequencyEph_h

#include <string>
#include <utility>

#include "pulsarDb/PulsarEph.h"

#include "tip/Table.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/TimeSystem.h"

namespace tip {
  class Header;
}

namespace st_stream {
  class OStream;
}

namespace pulsarDb {

  /** \class FrequencyEph
      \brief Class representing a single pulsar ephemeris expressed with three frequency coefficients.
             Note: f0, f1, f2 depend on the time system, while AbsoluteTime objects don't.
  */
  class FrequencyEph : public PulsarEph {
    public:
      /** \brief Create a pulsar ephemeris object with the given parameters.
          \param valid_since Beginning of time period during which this ephemeris is considered valid.
          \param valid_until End of time period during which this ephemeris is considered valid.
          \param epoch Reference epoch of frequency parameters (f0, f1, and f2).
          \param ra Right Ascension of the pulsar.
          \param dec Declination of the pulsar.
          \param phi0 Pulse phase at the given epoch.
          \param f0 Pulse frequency at the given epoch.
          \param f1 First time derivative of pulse frequency at the given epoch.
          \param f2 Second time derivative of pulse frequency at the given epoch.
      */
      FrequencyEph(const std::string & time_system_name, const timeSystem::AbsoluteTime & valid_since,
        const timeSystem::AbsoluteTime & valid_until, const timeSystem::AbsoluteTime & epoch, double ra, double dec,
        double phi0, double f0, double f1, double f2): m_system(&timeSystem::TimeSystem::getSystem(time_system_name)),
        m_since(valid_since), m_until(valid_until), m_epoch(epoch), m_ra(ra), m_dec(dec), m_phi0(phi0), m_f0(f0), m_f1(f1), m_f2(f2) {}

      /** \brief Create a pulsar ephemeris object with the parameters stored in tip record.
          \param time_system_name Name of time system to interpret frequency parameters, such as "TDB" or "UTC".
          \param record Record that stores all parameters for an ephemeris being created.
      */
      FrequencyEph(const tip::Table::ConstRecord & record, const tip::Header & header);

      /// \brief Destruct this FrequencyEph object.
      virtual ~FrequencyEph() {}

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

      /// \brief Create a copy of this object, and return a pointer to it.
      virtual PulsarEph * clone() const { return new FrequencyEph(*this); }

      /** \brief Compute a pulse phase of a given time, and return it.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time for which a pulse phase is to be computed.
          \param phase_offset Value to be added to a pulse phase. This value is added to the computed pulse phase
                 before truncated to a value in range [0, 1).
      */
      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      /** \brief Compute a pulse frequency at a given time, and return it.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time at which a pulse frequency is to be computed.
          \param derivative_order Order of frequency derivative to be computed. Set zero (0) to this argument
                 to compute a pulse frequency, one (1) the first time derivative of pulse frequency, and so on.
      */
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const;

      /** \brief Compute a sky position at a given time, and return it. The returned value is a pair of Right Ascension
                 and Declination of the position.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Absolute time at which a sky position is to be computed.
      */
      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const;

    protected:
      /// \brief Output text expression of subclass-specific parameters of this object to a given output stream.
      virtual void writeModelParameter(st_stream::OStream & os) const;

    private:
      /** \brief Compute the number of elapsed seconds between given absolute times in the time system in which this ephemeris
                 is defined, and return it.
          \param at1 Absolute time at which the number of elapsed seconds is to be computed.
          \param at2 Absolute time that gives a start time to count the number of elapsed seconds.
      */
      virtual double calcElapsedSecond(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const;

      /** \brief Compute the number of elapsed seconds at a given absolute time in the time system in which this ephemeris
                 is defined, and return it.
          \param at Absolute time at which the number of elapsed seconds is to be computed.
      */
      virtual double calcElapsedSecond(const timeSystem::AbsoluteTime & at) const { return calcElapsedSecond(at, m_epoch); }

      const timeSystem::TimeSystem * m_system;
      timeSystem::AbsoluteTime m_since;
      timeSystem::AbsoluteTime m_until;
      timeSystem::AbsoluteTime m_epoch;
      double m_ra;
      double m_dec;
      double m_phi0;
      double m_f0;
      double m_f1;
      double m_f2;
  };
}

#endif
