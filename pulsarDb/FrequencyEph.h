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

      virtual ~FrequencyEph() {}

      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      virtual const timeSystem::AbsoluteTime & getValidSince() const { return m_since; }

      virtual const timeSystem::AbsoluteTime & getValidUntil() const { return m_until; }

      virtual const timeSystem::AbsoluteTime & getEpoch() const { return m_epoch; }

      virtual PulsarEph * clone() const { return new FrequencyEph(*this); }

      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const;

      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const;

      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const;

    protected:
      virtual void writeModelParameter(st_stream::OStream & os) const;

    private:
      virtual double calcElapsedSecond(const timeSystem::AbsoluteTime & at1, const timeSystem::AbsoluteTime & at2) const;
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
