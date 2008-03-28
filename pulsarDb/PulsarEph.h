/** \file PulsarEph.h
    \brief Interface for PulsarEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarEph_h
#define pulsarDb_PulsarEph_h

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "pulsarDb/FormattedEph.h"

#include "st_stream/Stream.h"
#include "tip/Table.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/Duration.h"
#include "timeSystem/TimeSystem.h"

namespace tip {
  class Header;
}

namespace pulsarDb {

  /** \class PulsarEph
      \brief Base class representing a single pulsar ephemeris.
  */
  class PulsarEph: public FormattedEph {
    public:
      virtual ~PulsarEph() {}

      /// \brief Return a time system to be used to interpret return values of calcFrequency method.
      virtual const timeSystem::TimeSystem & getSystem() const = 0;

      /** \brief Return a start time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the start time itself.
      */
      virtual const timeSystem::AbsoluteTime & getValidSince() const = 0;

      /** \brief Return an end time of a time interval, during which this ephemeris is considered valid.
                 Note: The ephemeris is also considered valid on the end time itself.
      */
      virtual const timeSystem::AbsoluteTime & getValidUntil() const = 0;

      /// \brief Return a reference epoch of this ephemeris.
      virtual const timeSystem::AbsoluteTime & getEpoch() const = 0;

      /** \brief Compute the spin phase of the given time.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Time of the event.
          \param phase_offset Phase value to be added to the computed pulse phase.
      */
      virtual double calcPulsePhase(const timeSystem::AbsoluteTime & ev_time, double phase_offset = 0.) const = 0;

      /** \brief Compute the pulse frequency at a given time in the time system given to this object
                 upon its construction. Call getSytem method to obtain the time system to interpret
                 the return value of this method. The unit of frequency and its derivatives must be
                 based on seconds, i.e., Hz (sE-1), Hz/s (sE-2), and so on.
                 Note: validity of the ephemeris (valid since and valid until) are not checked. 
          \param ev_time Time of the event.
          \param derivative_order Order of time derivative of frequency to compute. Set 0 (zero) to obtain
                 frequency, 1 (one) the first time derivative of frequency, and so on.
      */
      virtual double calcFrequency(const timeSystem::AbsoluteTime & ev_time, int derivative_order = 0) const = 0;

      /** \brief Compute the Right Ascension and Declination at a given time.
                 Note: validity of the ephemeris (valid since and valid until) are not checked.
          \param ev_time Time of the event.
      */
      virtual std::pair<double, double> calcSkyPosition(const timeSystem::AbsoluteTime & ev_time) const = 0;

      /// \brief Create a copy of this object.
      virtual PulsarEph * clone() const = 0;

      /// \brief Output text expression of this PulsarEph to a given output stream.
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    protected:
      PulsarEph() {};

      /// \brief Output text expression of subclass-specific parameters of this PulsarEph to a given output stream.
      virtual void writeModelParameter(st_stream::OStream & os) const = 0;
  };

  inline st_stream::OStream & operator <<(st_stream::OStream & os, const PulsarEph & eph) { return eph.write(os); }

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

  /** \class PeriodEph
      \brief Class representing a single pulsar ephemeris expressed with three period coefficients.
             Note: p0, p1, p2 depend on the time system, while AbsoluteTime objects don't.
  */
  class PeriodEph : public PulsarEph {
    public:
      /** \brief Create a pulsar ephemeris object with the given parameters.
          \param valid_since Beginning of time period during which this ephemeris is considered valid.
          \param valid_until End of time period during which this ephemeris is considered valid.
          \param epoch Reference epoch of frequency parameters (f0, f1, and f2).
          \param ra Right Ascension of the pulsar.
          \param dec Declination of the pulsar.
          \param phi0 Pulse phase at the given epoch.
          \param f0 Pulse period at the given epoch.
          \param f1 First time derivative of pulse period at the given epoch.
          \param f2 Second time derivative of pulse period at the given epoch.
      */
      PeriodEph(const std::string & time_system_name, const timeSystem::AbsoluteTime & valid_since,
        const timeSystem::AbsoluteTime & valid_until, const timeSystem::AbsoluteTime & epoch, double ra, double dec,
        double phi0, double p0, double p1, double p2): m_system(&timeSystem::TimeSystem::getSystem(time_system_name)),
        m_since(valid_since), m_until(valid_until), m_epoch(epoch), m_ra(ra), m_dec(dec), m_phi0(phi0), m_p0(p0), m_p1(p1), m_p2(p2) {}

      virtual ~PeriodEph() {}

      virtual const timeSystem::TimeSystem & getSystem() const { return *m_system; }

      virtual const timeSystem::AbsoluteTime & getValidSince() const { return m_since; }

      virtual const timeSystem::AbsoluteTime & getValidUntil() const { return m_until; }

      virtual const timeSystem::AbsoluteTime & getEpoch() const { return m_epoch; }

      virtual PulsarEph * clone() const { return new PeriodEph(*this); }

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
      double m_p0;
      double m_p1;
      double m_p2;
  };

  typedef std::vector<PulsarEph *> PulsarEphCont;
}

#endif
